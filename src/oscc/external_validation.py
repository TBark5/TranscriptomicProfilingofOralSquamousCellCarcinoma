"""Independent paired cohort analysis; original discovery bundle is read-only."""
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .acquisition import acquire, sha256
from .dashboard import load_results
from .datasets import get_dataset
from .differential import fit_paired, classify
from .enrichment import run_enrichment, HALLMARK_SHA
from .geo_validation import read_geo_metadata, read_geo_counts
from .preprocessing import filter_counts, log_expression
from .visualization import make_figures


def compare_genes(discovery, validation, raw_symbols, retained_symbols, alpha=.05, lfc=1., validation_alpha=.05):
    """Preserve every discovery row and explicitly distinguish unavailable tests."""
    columns = ["log2FoldChange", "lfcSE", "pvalue", "padj", "test_status"]
    left = discovery[columns].add_prefix("discovery_")
    right = validation[columns].add_prefix("validation_")
    compared = left.join(right, how="left")
    compared["discovery_candidate"] = classify(discovery, alpha, lfc) != "not_significant"
    compared["mapping_status"] = "not_in_validation_source"
    compared.loc[compared.index.isin(raw_symbols), "mapping_status"] = "filtered_low_expression"
    compared.loc[compared.index.isin(retained_symbols), "mapping_status"] = "matched_exact_symbol"
    compared["testable"] = (compared.mapping_status.eq("matched_exact_symbol") &
                            np.isfinite(compared.discovery_pvalue) & np.isfinite(compared.validation_pvalue))
    compared["same_direction"] = compared.testable & (
        np.sign(compared.discovery_log2FoldChange) == np.sign(compared.validation_log2FoldChange))
    compared["replicated"] = (compared.discovery_candidate & compared.same_direction &
                              (compared.validation_padj < validation_alpha))
    compared["replicated_large_effect"] = compared.replicated & (compared.validation_log2FoldChange.abs() > lfc)
    compared["validation_outcome"] = "not_a_discovery_candidate"
    candidate = compared.discovery_candidate
    compared.loc[candidate, "validation_outcome"] = "not_testable"
    compared.loc[candidate & compared.testable, "validation_outcome"] = "not_significant"
    compared.loc[candidate & compared.testable & (compared.validation_padj < validation_alpha) & ~compared.same_direction,
                 "validation_outcome"] = "significant_opposite_direction"
    compared.loc[compared.replicated, "validation_outcome"] = "replicated"
    compared["lfc_difference"] = compared.validation_log2FoldChange - compared.discovery_log2FoldChange
    compared.index.name = "gene"
    return compared


def compare_pathways(discovery, validation, alpha=.05):
    keys = ["library", "Term"]
    columns = keys + ["NES", "pvalue", "padj"]
    left = discovery[columns].rename(columns={c: f"discovery_{c}" for c in columns if c not in keys})
    right = validation[columns].rename(columns={c: f"validation_{c}" for c in columns if c not in keys})
    table = left.merge(right, on=keys, how="outer", validate="one_to_one")
    table["testable"] = table.discovery_padj.notna() & table.validation_padj.notna()
    table["same_direction"] = table.testable & (np.sign(table.discovery_NES) == np.sign(table.validation_NES))
    table["discovery_significant"] = table.discovery_padj < alpha
    table["replicated"] = table.discovery_significant & table.same_direction & (table.validation_padj < alpha)
    return table


def preservation_hashes(root):
    root = Path(root)
    paths = [root / "REPORT.md"]
    for directory in ["results", "data/processed"]:
        paths.extend(p for p in (root / directory).rglob("*") if p.is_file() and
                     not any(part in {"GSE20116", "GSE184616"} for part in p.relative_to(root).parts))
    return {p.relative_to(root).as_posix(): sha256(p) for p in paths if p.is_file()}


def comparison_figures(comparison, pathways, figures):
    import matplotlib.pyplot as plt
    from .visualization import save as save_figure

    # Use all mutually testable genes, not a validation-selected subset.
    common = comparison.loc[comparison.testable]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(common.discovery_log2FoldChange, common.validation_log2FoldChange,
               s=8, alpha=.18, color="#627A88", rasterized=True)
    candidates = common.loc[common.discovery_candidate]
    ax.scatter(candidates.discovery_log2FoldChange, candidates.validation_log2FoldChange,
               s=10, alpha=.4, color="#A7443C", label="Discovery candidates")
    ax.axhline(0, color="grey", lw=.7); ax.axvline(0, color="grey", lw=.7)
    ax.set(xlabel="GSE20116 log2 fold change", ylabel="GSE184616 log2 fold change",
           title="Independent cohort effect-size comparison")
    ax.legend(frameon=False)
    save_figure(fig, figures, "effect_size_comparison")

    # These six were highlighted in the original README before validation.
    featured = comparison.reindex(["PTHLH", "LAMC2", "COL4A6", "TMPRSS11B", "PTGFR", "PYGM"])
    fig, ax = plt.subplots(figsize=(8, 4.7))
    for offset, prefix, color, label in [(-.13, "discovery", "#627A88", "GSE20116"),
                                        (.13, "validation", "#A7443C", "GSE184616")]:
        ax.errorbar(featured[f"{prefix}_log2FoldChange"], np.arange(len(featured)) + offset,
                    xerr=1.96 * featured[f"{prefix}_lfcSE"], fmt="o", capsize=3, color=color, label=label)
    ax.set_yticks(np.arange(len(featured)), featured.index)
    ax.axvline(0, color="grey", lw=.8)
    ax.set(xlabel="Tumor / normal log2 fold change (approximate 95% Wald interval)",
           title="Candidates highlighted before external validation")
    ax.legend(frameon=False)
    save_figure(fig, figures, "candidate_effects")

    common_pathways = pathways.loc[pathways.testable]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(common_pathways.discovery_NES, common_pathways.validation_NES,
               c=np.where(common_pathways.replicated, "#A7443C", "#627A88"), s=35)
    for row in common_pathways.itertuples():
        if row.Term in ["Myogenesis", "E2F Targets"]:
            ax.annotate(row.Term, (row.discovery_NES, row.validation_NES), xytext=(5, 5), textcoords="offset points")
    ax.axhline(0, color="grey", lw=.7); ax.axvline(0, color="grey", lw=.7)
    ax.set(xlabel="GSE20116 Hallmark NES", ylabel="GSE184616 Hallmark NES",
           title="Hallmark replication (red: both FDR < 0.05, same direction)")
    save_figure(fig, figures, "pathway_replication")


def run_validation(root, offline=False, seed=42, permutations=1000):
    from .validation_report import write_validation_report
    if permutations < 100:
        raise ValueError("At least 100 GSEA permutations required")
    root = Path(root).resolve()
    cohort = get_dataset("GSE184616")
    output = cohort.output(root)
    tables, figures = output / "tables", output / "figures"
    for directory in [tables, figures]:
        directory.mkdir(parents=True, exist_ok=True)
    original_hashes = preservation_hashes(root)
    manifest_path = output / "manifest.json"
    report_path = root / "VALIDATION_REPORT.md"
    # Prevent a stale successful report surviving a failed subsequent attempt.
    report_path.write_text("# Independent validation\n\nStatus: running. Inspect results/GSE184616/manifest.json.\n", encoding="utf-8")
    manifest = {"status": "running", "dataset": cohort.accession, "role": cohort.role, "data_kind": "real",
                "started_utc": datetime.now(timezone.utc).isoformat(),
                "original_snapshot_sha256": original_hashes,
                "config": {"alpha": .05, "lfc": 1., "min_count": 10, "min_samples": cohort.min_samples,
                           "seed": seed, "permutations": permutations, "libraries": ["MSigDB_Hallmark_2020"],
                           "mapping": "exact source gene symbol; no alias guesses or count pooling",
                           "replication": "discovery candidate; same LFC sign; validation genome-wide BH < 0.05"}}
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    try:
        np.random.seed(seed)
        discovery = load_results(root, "GSE20116", legacy=True)
        manifest["discovery_input_sha256"] = {"manifest": sha256(root / "results/manifest.json"),
                    "differential_expression": sha256(root / "results/tables/differential_expression.csv")}
        discovery_config = discovery["manifest"]["config"]
        frozen = discovery["differential_expression"].loc[
            classify(discovery["differential_expression"], discovery_config["alpha"], discovery_config["lfc"]) != "not_significant"]
        frozen.to_csv(tables / "frozen_discovery_candidates.csv", lineterminator="\n")
        manifest["candidate_selection"] = {"count": len(frozen), "alpha": discovery_config["alpha"],
            "lfc": discovery_config["lfc"], "sha256": sha256(tables / "frozen_discovery_candidates.csv"),
            "frozen_utc": datetime.now(timezone.utc).isoformat()}
        print("Verifying GSE184616 source files and matched sample metadata...", flush=True)
        manifest["sources"] = acquire(root, offline, cohort.accession)
        metadata = read_geo_metadata(cohort.raw(root) / "GSE184616_family.soft.gz")
        raw_counts = read_geo_counts(cohort.raw(root) / "GSE184616_unnormalisedGeneCounts.txt.gz", metadata)
        if set(metadata.index) & set(discovery["metadata"].index):
            raise ValueError("Discovery and validation share GEO sample accessions")
        counts, audit = filter_counts(raw_counts, min_count=10, min_samples=cohort.min_samples)
        metadata.to_csv(tables / "metadata.csv", lineterminator="\n")
        audit.to_csv(tables / "filter_audit.csv", index_label="gene", lineterminator="\n")
        counts.to_csv(tables / "counts.csv", lineterminator="\n")
        manifest.update(source_genes=len(raw_counts), samples=len(metadata), pairs=int(metadata.patient_id.nunique()))
        print(f"Verified {len(raw_counts):,} integer count rows; fitting {len(counts):,} genes, 15 pairs...", flush=True)
        result, normalized, sizes, convergence, diagnostics = fit_paired(counts, metadata)
        result.index.name = "gene"
        result["direction"] = classify(result)
        log_expr = log_expression(normalized)
        for name, frame in [("differential_expression", result), ("normalized_counts", normalized),
                            ("log_expression", log_expr), ("convergence", convergence)]:
            frame.to_csv(tables / f"{name}.csv", lineterminator="\n")
        sizes.to_csv(tables / "size_factors.csv", lineterminator="\n")
        pd.DataFrame({"retained_count_sum": counts.sum(), "size_factor": sizes}).to_csv(tables / "library_sizes.csv", lineterminator="\n")
        manifest["diagnostics"] = diagnostics
        (tables / "model_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
        comparison = compare_genes(discovery["differential_expression"], result, raw_counts.index, counts.index,
                                   discovery_config["alpha"], discovery_config["lfc"])
        comparison.to_csv(tables / "gene_concordance.csv", lineterminator="\n")
        candidates = comparison.loc[comparison.discovery_candidate]
        candidates.to_csv(tables / "candidate_concordance.csv", lineterminator="\n")
        print("Computing independent Hallmark enrichment and replication...", flush=True)
        discovery_hallmark = [e for e in discovery["manifest"]["enrichment"] if e["library"] == "MSigDB_Hallmark_2020"]
        if not discovery_hallmark or discovery_hallmark[0].get("source", {}).get("sha256") != HALLMARK_SHA:
            raise ValueError("Discovery Hallmark definition differs from validation")
        ora, gsea, status, curve = run_enrichment(result, root, tables, ["MSigDB_Hallmark_2020"],
                                                 seed=seed, permutations=permutations, offline=offline)
        if not status or any(e["status"] != "complete" for e in status) or gsea.empty:
            raise RuntimeError(f"Pathway validation did not complete: {status}")
        manifest["enrichment"] = status
        pathways = compare_pathways(discovery["gsea"].loc[discovery["gsea"].library == "MSigDB_Hallmark_2020"], gsea)
        pathways.to_csv(tables / "pathway_replication.csv", index=False, lineterminator="\n")
        pca, explained, correlations = make_figures(counts, log_expr, metadata, result, gsea, curve,
                                                    figures, dataset=cohort.accession)
        pca.to_csv(tables / "pca.csv", lineterminator="\n")
        correlations.to_csv(tables / "sample_correlations.csv", lineterminator="\n")
        comparison_figures(comparison, pathways, figures)
        with (figures / "CAPTIONS.md").open("a", encoding="utf-8") as handle:
            handle.write("\n\neffect_size_comparison: All mutually testable exact-symbol genes; highlighted candidates were selected only in discovery.\n\ncandidate_effects: Six genes highlighted in the original README, shown regardless of external outcome; unshrunk approximate 95% Wald intervals.\n\npathway_replication: Matching Hallmark terms, same NES sign and FDR < 0.05 in both cohorts highlighted. Gene backgrounds differ.\n")
        common = comparison.loc[comparison.testable]
        rho = float(spearmanr(common.discovery_log2FoldChange, common.validation_log2FoldChange).statistic)
        manifest["qc"] = {"pca_variance": explained, "library_min": int(counts.sum().min()), "library_max": int(counts.sum().max())}
        manifest["summary"] = {"genes": len(result), "tested": int(result.pvalue.notna().sum()),
            "up": int((result.direction == "up").sum()), "down": int((result.direction == "down").sum()),
            "discovery_candidates": len(candidates), "testable_candidates": int(candidates.testable.sum()),
            "concordant_candidates": int(candidates.same_direction.sum()), "replicated_candidates": int(candidates.replicated.sum()),
            "replicated_large_effect": int(candidates.replicated_large_effect.sum()),
            "common_tested_genes": len(common), "lfc_spearman": rho,
            "discovery_significant_pathways": int(pathways.discovery_significant.sum()),
            "testable_discovery_pathways": int((pathways.discovery_significant & pathways.testable).sum()),
            "replicated_pathways": int(pathways.replicated.sum())}
        manifest["versions"] = {name: importlib.metadata.version(name) for name in
                                ["numpy", "pandas", "scipy", "pydeseq2", "gseapy", "matplotlib"]}
        manifest["python"] = platform.python_version()
        if preservation_hashes(root) != original_hashes:
            raise RuntimeError("Original discovery snapshot changed during validation")
        manifest["original_snapshot_unchanged"] = True
        write_validation_report(root, manifest, comparison, pathways, metadata)
        paths = [p for p in output.rglob("*") if p.is_file() and p != manifest_path] + [report_path]
        manifest["output_sha256"] = {p.relative_to(root).as_posix(): sha256(p) for p in paths}
        manifest["status"] = "complete"
        manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
        print(json.dumps(manifest["summary"], indent=2), flush=True)
    except Exception as exc:
        manifest["status"], manifest["error"] = "failed", str(exc)
        report_path.write_text(f"# Independent validation\n\nStatus: **failed**.\n\nBlocker: {exc}\n\n"
                               "Partial outputs must not be interpreted as completed validation. See results/GSE184616/manifest.json.\n", encoding="utf-8")
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
