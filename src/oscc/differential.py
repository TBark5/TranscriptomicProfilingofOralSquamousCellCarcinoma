"""Paired negative-binomial inference with convergence auditing."""
import warnings
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from .validation import validate_counts, design_matrix


def bh_adjust(pvalues):
    p = pd.Series(pvalues, dtype=float)
    finite = p.notna()
    if ((p[finite] < 0) | (p[finite] > 1) | ~np.isfinite(p[finite])).any():
        raise ValueError("P-values outside [0,1]")
    result = pd.Series(np.nan, index=p.index)
    if finite.any():
        result.loc[finite] = multipletests(p[finite], method="fdr_bh")[1]
    return result


def classify(result, alpha=.05, lfc=1.):
    if not 0 < alpha < 1 or lfc < 0:
        raise ValueError("Require 0 < alpha < 1 and lfc >= 0")
    significant = (result.padj < alpha) & (result.log2FoldChange.abs() > lfc)
    return pd.Series(np.where(significant, np.where(result.log2FoldChange > 0, "up", "down"), "not_significant"), index=result.index)


def fit_paired(counts, metadata):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    counts, metadata = validate_counts(counts, metadata)
    design_matrix(metadata)
    metadata = metadata.copy()
    metadata["patient_id"] = metadata.patient_id.astype(str)
    metadata["condition"] = pd.Categorical(metadata.condition, categories=["normal", "tumor"])
    inference = DefaultInference(n_cpus=1)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        dds = DeseqDataSet(counts=counts.T, metadata=metadata, design="~patient_id + condition",
                           refit_cooks=False, inference=inference, quiet=True)
        actual_design = dds.obsm["design_matrix"]
        rank = int(np.linalg.matrix_rank(actual_design))
        if rank != actual_design.shape[1] or rank >= len(metadata):
            raise ValueError("PyDESeq2 design is rank deficient or saturated")
        dds.deseq2()
        stats = DeseqStats(dds, contrast=["condition", "tumor", "normal"],
                          independent_filter=False, cooks_filter=True, inference=inference, quiet=True)
        stats.summary()
    result = stats.results_df.copy()
    # Optimizer failure is not evidence of differential expression.
    flags = ["_genewise_converged", "_MAP_converged", "_LFC_converged"]
    missing = set(flags) - set(dds.var.columns)
    if missing:
        raise RuntimeError(f"Cannot audit convergence in this PyDESeq2 version: {missing}")
    audit = dds.var[flags + ["dispersions"]].copy()
    audit["model_converged"] = audit[flags].fillna(False).all(axis=1)
    result["pvalue_model"] = result.pvalue
    result["model_converged"] = audit.model_converged
    result.loc[~result.model_converged, "pvalue"] = np.nan
    result["padj"] = bh_adjust(result.pvalue)
    result["test_status"] = np.where(~result.model_converged, "optimizer_not_converged",
                                    np.where(result.pvalue.isna(), "cooks_or_undefined", "tested"))
    normalized = pd.DataFrame(dds.layers["normed_counts"].T, index=counts.index, columns=counts.columns)
    size_factors = pd.Series(dds.obs["size_factors"], index=counts.columns, name="size_factor")
    if not np.isfinite(normalized.to_numpy()).all() or not (size_factors > 0).all():
        raise RuntimeError("Invalid fitted normalization")
    diagnostics = {"design": "~patient_id + condition", "contrast": "tumor / normal", "design_rank": rank,
                   "residual_df": len(metadata) - rank, "failed_convergence": int((~audit.model_converged).sum()),
                   "missing_pvalue": int(result.pvalue.isna().sum()),
                   "warnings": sorted(set(str(w.message) for w in captured)),
                   "normalization": "DESeq2 median-of-ratios size factors", "cooks_filter": True,
                   "outlier_replacement": False, "independent_filter": False,
                   "bh_family": "all finite, convergence-qualified, Cook-filtered p-values"}
    return result, normalized, size_factors, audit, diagnostics
