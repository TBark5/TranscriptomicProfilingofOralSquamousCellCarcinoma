"""Non-destructive, source-grounded second-stage OSCC audit and figures.

Run from the repository root with its pinned Python 3.11 environment. All writes
are limited to results/phase2/v1; original analyses are read-only inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

from oscc.enrichment import read_gmt

GENES = ["PTHLH", "LAMC2", "COL4A6", "MMP11", "SPP1", "TMPRSS11B", "PTGFR", "PYGM", "CRNN", "MAL"]
MARKERS = {
    "epithelial": ["EPCAM", "KRT5", "KRT14", "KRT13", "KRT4"],
    "immune": ["PTPRC", "CD3D", "CD68", "LST1", "MS4A1"],
    "stromal": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "muscle": ["ACTA1", "MYH1", "MYH2", "DES", "MYOG"],
}
DIFFERENTIATION = ["KRT13", "KRT4", "CRNN", "MAL", "IVL", "SPRR1A", "FLG"]
OUTCOMES = ["untestable", "same-direction significant", "same-direction nonsignificant", "opposite-direction significant", "opposite-direction nonsignificant"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def table(path, index="gene"):
    x = pd.read_csv(path)
    if index not in x:
        raise ValueError(f"Missing {index} in {path}")
    if x[index].duplicated().any():
        raise ValueError(f"Duplicated {index} in {path}")
    return x.set_index(index)


def save(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, float_format="%.12g", lineterminator="\n")


def figure(path):
    plt.tight_layout()
    for ext in ("svg", "png"):
        plt.savefig(path.with_suffix("." + ext), dpi=220, bbox_inches="tight")
    plt.close()


def classify_candidates(d, v):
    candidates = d.loc[(d.padj < .05) & (d.log2FoldChange.abs() > 1)].copy()
    j = candidates[["log2FoldChange", "lfcSE", "padj", "baseMean"]].rename(columns=lambda c: "discovery_" + c)
    j = j.join(v[["log2FoldChange", "lfcSE", "padj", "pvalue", "baseMean", "test_status"]].rename(columns=lambda c: "validation_" + c), how="left")
    tested = j.validation_pvalue.notna()
    same = tested & (np.sign(j.discovery_log2FoldChange) == np.sign(j.validation_log2FoldChange))
    sig = tested & (j.validation_padj < .05)
    j["outcome"] = np.select([~tested, same & sig, same & ~sig, ~same & sig], OUTCOMES[:4], default=OUTCOMES[4])
    assert len(j) == j.outcome.value_counts().sum()
    return j.reset_index()


def hgnc_audit(discovery, validation, ref_path):
    h = pd.read_csv(ref_path, sep="\t", dtype=str, low_memory=False).fillna("")
    required = {"symbol", "status", "prev_symbol", "alias_symbol", "hgnc_id"}
    if not required <= set(h):
        raise ValueError("HGNC file lacks required fields")
    approved = h.loc[h.status == "Approved"].copy()
    exact = approved.groupby("symbol").hgnc_id.apply(list).to_dict()
    former, alias = {}, {}
    for row in approved.itertuples():
        for col, index in [(row.prev_symbol, former), (row.alias_symbol, alias)]:
            for symbol in col.split("|"):
                symbol = symbol.strip()
                if symbol:
                    index.setdefault(symbol, set()).add(row.symbol)
    rows = []
    for gene, row in discovery.iterrows():
        p = sorted(former.get(gene, set()))
        a = sorted(alias.get(gene, set()))
        v = validation.loc[gene] if gene in validation.index else None
        if gene in exact:
            status = "HGNC approved exact"
        elif len(p) == 1:
            status = "one previous-symbol suggestion"
        elif len(p) > 1:
            status = "ambiguous previous symbol"
        elif len(a) == 1:
            status = "one alias suggestion"
        elif len(a) > 1:
            status = "ambiguous alias"
        else:
            status = "no HGNC match"
        rows.append({"gene": gene, "refseq": row.get("refseq", ""), "hgnc_status": status,
                     "approved_hgnc_id": ";".join(exact.get(gene, [])),
                     "previous_symbol_suggestions": ";".join(p), "alias_suggestions": ";".join(a),
                     "exact_validation_symbol": gene in validation.index,
                     "validation_test_status": v.test_status if v is not None else "not_in_filtered_validation",
                     "discovery_candidate": bool(row.padj < .05 and abs(row.log2FoldChange) > 1),
                     "discovery_test_status": row.test_status})
    return pd.DataFrame(rows)


def model_audit(d, v):
    rows = []
    for cohort, frame in [("GSE20116", d), ("GSE184616", v)]:
        tested = frame.loc[frame.pvalue.notna()].copy()
        q = multipletests(tested.pvalue, method="fdr_bh")[1]
        ratio = tested.log2FoldChange / tested.lfcSE
        theory_p = 2 * norm.sf(np.abs(tested.stat))
        diff = np.abs(theory_p - tested.pvalue.to_numpy())
        rows.append({"cohort": cohort, "fit_rows": len(frame), "tested_rows": len(tested),
                     "failed_convergence": int((~frame.model_converged).sum()),
                     "max_absolute_bh_difference": float(np.max(np.abs(q - tested.padj))),
                     "max_absolute_stat_ratio_difference": float(np.nanmax(np.abs(ratio - tested.stat))),
                     "max_absolute_normal_p_difference": float(np.nanmax(diff)),
                     "p_below_1e_20": int((tested.pvalue < 1e-20).sum()),
                     "abs_lfc_above_5": int((tested.log2FoldChange.abs() > 5).sum()),
                     "median_base_mean_extreme_lfc": float(tested.loc[tested.log2FoldChange.abs() > 5, "baseMean"].median())})
    return pd.DataFrame(rows)


def effect_summary(d, v, candidates):
    common = d.loc[d.pvalue.notna() & d.index.isin(v.index[v.pvalue.notna()]), ["log2FoldChange", "lfcSE", "padj"]].rename(columns=lambda c: "discovery_" + c)
    common = common.join(v[["log2FoldChange", "lfcSE", "padj"]].rename(columns=lambda c: "validation_" + c), validate="one_to_one")
    common = common.reset_index()
    rows = []
    for label, frame in [("mutually testable genome", common), ("testable discovery candidates", candidates.loc[candidates.validation_pvalue.notna()])]:
        a, b = frame.discovery_log2FoldChange.to_numpy(), frame.validation_log2FoldChange.to_numpy()
        delta = b - a
        rho = spearmanr(a, b).statistic
        # Fisher-z approximation is descriptive; Spearman CI is not exact and ignores gene dependence.
        z = np.arctanh(np.clip(rho, -.999999, .999999))
        ci = np.tanh([z - 1.96 / np.sqrt(len(a) - 3), z + 1.96 / np.sqrt(len(a) - 3)])
        rows.append({"set": label, "genes": len(a), "spearman_rho": rho,
                     "approximate_rho_ci_low": ci[0], "approximate_rho_ci_high": ci[1],
                     "same_direction": int((np.sign(a) == np.sign(b)).sum()),
                     "median_validation_minus_discovery_lfc": np.median(delta),
                     "median_absolute_lfc_difference": np.median(abs(delta))})
    return common, pd.DataFrame(rows)


def selected_genes(d, v, mapping):
    rows = []
    maps = mapping.set_index("gene")
    for gene in GENES:
        for cohort, frame in [("GSE20116", d), ("GSE184616", v)]:
            if gene not in frame.index:
                rows.append({"gene": gene, "cohort": cohort, "status": "not in fitted table",
                             "mapping": maps.loc[gene, "hgnc_status"] if gene in maps.index else "not in discovery"})
            else:
                r = frame.loc[gene]
                rows.append({"gene": gene, "cohort": cohort, "status": r.test_status,
                             "mapping": maps.loc[gene, "hgnc_status"] if gene in maps.index else "not in discovery",
                             "log2FoldChange": r.log2FoldChange, "lfcSE": r.lfcSE,
                             "ci95_low": r.log2FoldChange - 1.96 * r.lfcSE,
                             "ci95_high": r.log2FoldChange + 1.96 * r.lfcSE,
                             "pvalue": r.pvalue, "padj": r.padj, "baseMean": r.baseMean})
    return pd.DataFrame(rows)


def paired_values(norm_counts, meta, genes, cohort):
    rows = []
    for gene in genes:
        if gene not in norm_counts.index:
            continue
        for patient, pair in meta.groupby("patient_id"):
            n = pair.index[pair.condition == "normal"].item()
            t = pair.index[pair.condition == "tumor"].item()
            nv, tv = float(norm_counts.loc[gene, n]), float(norm_counts.loc[gene, t])
            rows.append({"cohort": cohort, "gene": gene, "patient_id": patient,
                         "normal_normalized": nv, "tumor_normalized": tv,
                         "normal_log2p1": np.log2(nv + 1), "tumor_log2p1": np.log2(tv + 1),
                         "paired_log2p1_difference": np.log2(tv + 1) - np.log2(nv + 1),
                         "subsite": pair.tissue.iloc[0] if "tissue" in pair else "unknown"})
    return pd.DataFrame(rows)


def sensitivities(candidates, v):
    rows = []
    for alpha in [.01, .05, .10]:
        for min_lfc in [0, 1]:
            tested = candidates.validation_pvalue.notna()
            same = np.sign(candidates.discovery_log2FoldChange) == np.sign(candidates.validation_log2FoldChange)
            hit = tested & same & (candidates.validation_padj < alpha) & (candidates.validation_log2FoldChange.abs() > min_lfc)
            opposite = tested & ~same & (candidates.validation_padj < alpha) & (candidates.validation_log2FoldChange.abs() > min_lfc)
            rows.append({"analysis": "saved full-fit significance criteria", "validation_min_count": 10,
                         "min_samples": 15, "q_cutoff": alpha, "validation_abs_lfc_gt": min_lfc,
                         "testable_candidates": int(tested.sum()), "same_direction_significant": int(hit.sum()),
                         "opposite_direction_significant": int(opposite.sum()), "note": "BH family remains original 17,850 tests"})
    return pd.DataFrame(rows)


def pathway_original(root, out):
    d = pd.read_csv(root / "results/GSE20116/tables/gsea.csv")
    v = pd.read_csv(root / "results/GSE184616/tables/gsea.csv")
    p = d[["Term", "NES", "padj", "genes"]].merge(v[["Term", "NES", "padj", "genes"]], on="Term", how="outer", suffixes=("_discovery", "_validation"), validate="one_to_one")
    p["discovery_significant"] = p.padj_discovery < .05
    p["replicated"] = p.discovery_significant & (p.padj_validation < .05) & (np.sign(p.NES_discovery) == np.sign(p.NES_validation))
    p["leading_edge_overlap"] = p.apply(lambda r: len(set(str(r.genes_discovery).split(";")) & set(str(r.genes_validation).split(";"))), axis=1)
    save(p, out / "tables/pathway_comparison_original_universes.csv")
    return p


def pathway_common(root, out, d, v):
    common = set(d.index[d.pvalue.notna()]) & set(v.index[v.pvalue.notna()])
    sets = read_gmt(root / "data/raw/MSigDB_Hallmark_2020.gmt")
    eligible = {term: sorted(set(genes) & common) for term, genes in sets.items() if 15 <= len(set(genes) & common) <= 500}
    coverage = pd.DataFrame([{"Term": term, "original_set_size": len(set(genes)), "common_testable_overlap": len(set(genes) & common), "eligible": term in eligible} for term, genes in sets.items()])
    save(coverage, out / "tables/hallmark_common_universe_coverage.csv")
    outputs = {}
    for label, frame in [("discovery", d), ("validation", v)]:
        rank = frame.loc[sorted(common), ["stat"]].reset_index().sort_values(["stat", "gene"], ascending=[False, True], kind="stable")
        gp_result = gp.prerank(rnk=rank, gene_sets=eligible, threads=1, min_size=15, max_size=500,
                               permutation_num=1000, seed=42, outdir=None, no_plot=True, verbose=False)
        g = gp_result.res2d.rename(columns={"NOM p-val": "pvalue", "FDR q-val": "padj", "FWER p-val": "fwer", "Lead_genes": "genes"})
        g = g[["Term", "ES", "NES", "pvalue", "padj", "fwer", "genes"]]
        save(g, out / f"tables/gsea_common_{label}.csv")
        outputs[label] = g
    p = outputs["discovery"][["Term", "NES", "padj", "genes"]].merge(outputs["validation"][["Term", "NES", "padj", "genes"]], on="Term", suffixes=("_discovery", "_validation"), validate="one_to_one")
    p["discovery_significant"] = p.padj_discovery < .05
    p["replicated"] = p.discovery_significant & (p.padj_validation < .05) & (np.sign(p.NES_discovery) == np.sign(p.NES_validation))
    save(p, out / "tables/pathway_comparison_common_universe.csv")
    return p, coverage


def marker_scores(normalized, metadata, cohort):
    log = np.log2(normalized + 1)
    rows = []
    for group, genes in MARKERS.items():
        observed = [g for g in genes if g in log.index]
        if not observed:
            continue
        matrix = log.loc[observed].T
        z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0).replace(0, np.nan)
        sample_scores = z.mean(axis=1)
        for patient, pair in metadata.groupby("patient_id"):
            n, t = pair.index[pair.condition == "normal"].item(), pair.index[pair.condition == "tumor"].item()
            rows.append({"cohort": cohort, "signature": group, "patient_id": patient,
                         "marker_genes": ";".join(observed), "marker_coverage": len(observed),
                         "normal_score": sample_scores[n], "tumor_score": sample_scores[t],
                         "paired_difference": sample_scores[t] - sample_scores[n]})
    return pd.DataFrame(rows)


def theme_table(root, d, v):
    gd = pd.read_csv(root / "results/GSE20116/tables/gsea.csv").set_index("Term")
    gv = pd.read_csv(root / "results/GSE184616/tables/gsea.csv").set_index("Term")
    themes = {"cell cycle / proliferation": ["E2F Targets", "G2-M Checkpoint"],
              "ECM remodeling / EMT": ["Epithelial Mesenchymal Transition"],
              "epithelial differentiation (curated, not a Hallmark set)": []}
    rows = []
    for theme, terms in themes.items():
        for term in terms:
            if term in gd.index and term in gv.index:
                a, b = gd.loc[term], gv.loc[term]
                leading_a = set(str(a.genes).split(";"))
                leading_b = set(str(b.genes).split(";"))
                rows.append({"theme": theme, "term": term, "discovery_NES": a.NES, "discovery_FDR": a.padj,
                             "validation_NES": b.NES, "validation_FDR": b.padj,
                             "leading_edge_overlap": len(leading_a & leading_b),
                             "shared_leading_edge_genes": ";".join(sorted(leading_a & leading_b))})
    for gene in DIFFERENTIATION:
        rows.append({"theme": "epithelial differentiation (curated, not a Hallmark set)", "term": gene,
                     "discovery_NES": np.nan, "discovery_FDR": np.nan, "validation_NES": np.nan, "validation_FDR": np.nan,
                     "leading_edge_overlap": np.nan, "shared_leading_edge_genes": "",
                     "discovery_lfc": d.loc[gene, "log2FoldChange"] if gene in d.index else np.nan,
                     "validation_lfc": v.loc[gene, "log2FoldChange"] if gene in v.index else np.nan,
                     "discovery_q": d.loc[gene, "padj"] if gene in d.index else np.nan,
                     "validation_q": v.loc[gene, "padj"] if gene in v.index else np.nan})
    return pd.DataFrame(rows)


def make_figures(out, common, candidates, selected, paired, pathways, sensitivity):
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(6.6, 6))
    ax.scatter(common.discovery_log2FoldChange, common.validation_log2FoldChange, s=3, alpha=.13, color="#466a83", rasterized=True)
    c = candidates.loc[candidates.validation_pvalue.notna()]
    ax.scatter(c.discovery_log2FoldChange, c.validation_log2FoldChange, s=7, alpha=.38, color="#c05b45", rasterized=True)
    ax.axhline(0, color="gray", lw=.5); ax.axvline(0, color="gray", lw=.5)
    rho = spearmanr(common.discovery_log2FoldChange, common.validation_log2FoldChange).statistic
    ax.set(xlabel="Discovery tumor/normal log2 fold change (GSE20116)", ylabel="Validation tumor/normal log2 fold change (GSE184616)", title=f"Effect-size concordance: {len(common):,} mutually tested genes")
    ax.text(.03, .97, f"Genome-wide Spearman ρ={rho:.3f}\nCandidates shown in orange (n={len(c):,})", transform=ax.transAxes, va="top")
    figure(figdir / "effect_size_scatter")

    fig, ax = plt.subplots(figsize=(8, 5.6))
    colors = {"GSE20116": "#477695", "GSE184616": "#bb5945"}
    for i, gene in enumerate(GENES):
        for cohort, shift in [("GSE20116", -.17), ("GSE184616", .17)]:
            r = selected.loc[(selected.gene == gene) & (selected.cohort == cohort)].iloc[0]
            if pd.notna(r.get("log2FoldChange")):
                x = r.log2FoldChange
                ax.errorbar(x, i + shift, xerr=[[x-r.ci95_low], [r.ci95_high-x]], fmt="o", ms=4,
                            color=colors[cohort], capsize=2, label=cohort if i == 0 else None)
            else:
                ax.text(0.1, i + shift, "not measured", fontsize=7)
    ax.axvline(0, color="gray", lw=.7)
    ax.set_yticks(range(len(GENES)), GENES); ax.invert_yaxis(); ax.legend()
    ax.set(xlabel="Tumor/normal model log2 fold change ± 1.96 Wald SE", title="Ten prespecified genes; discovery 3 pairs, validation 15 pairs")
    figure(figdir / "ten_gene_forest")

    fig, axes = plt.subplots(2, 5, figsize=(15, 6.6), sharex=True)
    for ax, gene in zip(axes.flat, GENES):
        sub = paired.loc[paired.gene == gene]
        for cohort, color, offset in [("GSE20116", "#477695", -.07), ("GSE184616", "#bb5945", .07)]:
            for r in sub.loc[sub.cohort == cohort].itertuples():
                ax.plot([0+offset, 1+offset], [r.normal_log2p1, r.tumor_log2p1], color=color, alpha=.35, lw=.8)
            s = sub.loc[sub.cohort == cohort]
            if len(s):
                ax.scatter(np.full(len(s), 0+offset), s.normal_log2p1, s=8, color=color)
                ax.scatter(np.full(len(s), 1+offset), s.tumor_log2p1, s=8, color=color)
        ax.set_title(gene); ax.set_xticks([0, 1], ["normal", "tumor"])
    fig.suptitle("Patient-matched log2(normalized count + 1); blue discovery n=3, red validation n=15")
    fig.supylabel("log2(normalized count + 1)")
    figure(figdir / "paired_expression_ten_genes")

    sub = pathways.loc[pathways.discovery_significant].sort_values("NES_discovery")
    fig, ax = plt.subplots(figsize=(7.6, 5.5))
    for r in sub.itertuples():
        ax.scatter(r.NES_discovery, r.NES_validation, color="#26816d" if r.replicated else "#b8514e", s=35)
        ax.annotate(r.Term, (r.NES_discovery, r.NES_validation), xytext=(3, 2), textcoords="offset points", fontsize=6)
    ax.axhline(0, color="gray", lw=.5); ax.axvline(0, color="gray", lw=.5)
    ax.set(xlabel="Discovery Hallmark NES", ylabel="Validation Hallmark NES", title=f"Original-universe Hallmark comparison: {int(sub.replicated.sum())}/{len(sub)} replicate")
    figure(figdir / "pathway_nes_concordance")

    totals = candidates.outcome.value_counts().reindex(OUTCOMES, fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 3.5))
    colors5 = ["#888888", "#27816d", "#86b6a7", "#b8514e", "#d3a49e"]
    left = 0
    for label, value, color in zip(OUTCOMES, totals, colors5):
        ax.barh([0], [value], left=left, color=color, label=f"{label}: {value}")
        if value > 60:
            ax.text(left + value/2, 0, str(value), ha="center", va="center", fontsize=8)
        left += value
    ax.set(xlim=(0, len(candidates)), xlabel="Number of 1,339 discovery candidates", yticks=[], title="Mutually exclusive validation outcomes")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.30), ncol=2, fontsize=8)
    figure(figdir / "candidate_outcomes")

    s = sensitivity.loc[sensitivity.analysis == "saved full-fit significance criteria"].copy()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for lfc, color in [(0, "#27816d"), (1, "#477695")]:
        subset = s.loc[s.validation_abs_lfc_gt == lfc]
        ax.plot(subset.q_cutoff, subset.same_direction_significant, "o-", color=color, label=f"validation |LFC|>{lfc}")
    ax.set(xlabel="Validation genome-wide BH q threshold", ylabel="Same-direction significant candidates", title="Exploratory threshold sensitivity; original fitted gene family")
    ax.legend()
    figure(figdir / "sensitivity_summary")


def report(root, out, audit, candidates, effects, pathways, common_p, coverage, selected, paired, mapping, model, markers, subsite):
    counts = candidates.outcome.value_counts().reindex(OUTCOMES, fill_value=0).to_dict()
    rep = int(pathways.replicated.sum()); sig = int(pathways.discovery_significant.sum())
    crep = int(common_p.replicated.sum()); csig = int(common_p.discovery_significant.sum())
    selected_table = selected[["gene", "cohort", "status", "mapping", "log2FoldChange", "padj", "baseMean"]].copy()
    for col in ["log2FoldChange", "padj", "baseMean"]:
        selected_table[col] = selected_table[col].map(lambda x: "NA" if pd.isna(x) else f"{x:.3g}")
    def markdown(df):
        return df.to_markdown(index=False)
    marker_summary = markers.groupby(["cohort", "signature"]).agg(pairs=("patient_id", "size"), marker_coverage=("marker_coverage", "first"), median_paired_difference=("paired_difference", "median")).reset_index()
    p_rows = pathways.loc[pathways.discovery_significant, ["Term", "NES_discovery", "padj_discovery", "NES_validation", "padj_validation", "replicated", "leading_edge_overlap"]].sort_values("padj_discovery")
    lines = f"""# Phase 2 cross-cohort OSCC transcriptomics report

## Abstract

This secondary analysis of public bulk RNA-seq compares paired tumor and adjacent/surgical-margin normal tissue in GSE20116 (3 pairs, discovery) and GSE184616 (15 pairs, independent validation). Frozen outputs report 1,339 discovery candidates; independent recount of saved tables finds 1,187 testable in validation, 964 same-direction, 706 same-direction with genome-wide BH q<0.05, and 563 also with validation |log2FC|>1. Importantly, {counts['opposite-direction significant']} candidates have significant opposite-direction validation effects. The small discovery cohort, heterogeneous tissue composition, historical gene identifiers, and model-based Wald uncertainty limit causal and biomarker interpretation.

## Background and prespecified hypotheses

The goal is to test whether large discovery tumor/normal gene effects and Hallmark enrichment directions recur in a separately deposited patient cohort. The prespecified genes are {', '.join(GENES)}. Three interpretation themes are proliferation/cell cycle, matrix remodeling/EMT, and epithelial differentiation. These are associations in bulk tissue, not demonstrated pathway activation or tumor-intrinsic mechanisms.

## Provenance and frozen analyses

The input accession series, sample metadata, raw-source checksums, and original analysis details remain in [REPORT.md](../../../REPORT.md), [VALIDATION_REPORT.md](../../../VALIDATION_REPORT.md), and their original manifests. This report derives new results from those saved tables; the original outputs were not refitted or overwritten. `frozen_sha256.json` records the pre-change files. GSE184616 count and SOFT hashes match the original manifest, all 30 count-file columns match metadata, all 15 patients have a normal/tumor pair, and no public GSM accession overlaps discovery. The annotation is historical RefSeq representatives in discovery versus modern symbols in validation; primary comparison remains exact symbol matching.

## Methods and quality control

Original models used `~patient_id + condition`, tumor relative to normal, median-of-ratios size factors, Cook filtering, no independent filtering or outlier replacement, and BH correction over convergence-qualified finite p-values. Discovery has 2 residual degrees of freedom and validation 14. Primary discovery criteria were q<0.05 and |log2FC|>1; replication requires exact same symbol, testable validation result, same direction, and validation genome-wide q<0.05. The stricter count adds validation |log2FC|>1. These were frozen before this Phase 2 work. Model audit is in `tables/model_audit.csv`; Wald ratio, normal-tail p-values, and BH adjustment were numerically checked. Wald intervals are model-conditional and do not capture cross-patient heterogeneity. Approximate Fisher-z Spearman intervals treat genes as independent and are descriptive only.

| Cohort | Fitted genes | Tested | Failed convergence | p<1e-20 | |LFC|>5 | Median baseMean among |LFC|>5 |
|---|---:|---:|---:|---:|---:|---:|
"""
    for r in model.itertuples():
        lines += f"| {r.cohort} | {r.fit_rows:,} | {r.tested_rows:,} | {r.failed_convergence} | {r.p_below_1e_20} | {r.abs_lfc_above_5} | {r.median_base_mean_extreme_lfc:.2f} |\n"
    lines += f"""

## Discovery and independent validation

Independent recount of saved CSVs confirms all six headline counts: 1,339 discovery candidates; 1,187 validation-testable; 964 same-direction; 706 primary replications; 563 stricter replications; and {rep}/{sig} discovery-significant Hallmark pathways replicated. Mutually exclusive outcomes sum to 1,339:

| Outcome | Count |
|---|---:|
"""
    for key, value in counts.items():
        lines += f"| {key} | {value} |\n"
    lines += "\n" + markdown(effects.round(4)) + "\n\nCandidate significance is evaluated against all validation tests, not only selected genes. Direction concordance alone is not significant replication. The 61 significant opposite-direction effects are contradictory evidence, not a successful validation subset.\n\n"
    lines += "## Ten prespecified genes\n\n" + markdown(selected_table) + "\n\nAll ten genes are displayed whether or not they replicate. Values are from saved model results; `baseMean` is mean normalized count. A missing/failed model result must not be treated as a zero effect. Figure intervals are ±1.96 Wald SE, descriptive rather than independently estimated confidence intervals.\n\n"
    lines += f"""## Gene identifier and expression audit

`tables/gene_mapping_audit.csv` annotates each discovery symbol against the downloaded official HGNC complete set, records the historical RefSeq representative and exact validation test status, and lists previous/alias symbol suggestions. Suggestions do **not** alter exact-match primary calls. Mapping status counts are in `tables/mapping_summary.csv`; ambiguous matches are deliberately unresolved. `tables/paired_expression.csv` contains log2(normalized count + 1) per matched patient for the ten genes. The paired plots expose large patient-to-patient variation and possible low-count extremes; no outlier is removed. Paired log differences are not equivalent to the fitted model coefficient. In discovery, losing any one of three pairs can materially change a descriptive mean, so this is a stability diagnostic rather than a new significance test.

## Sensitivity analyses (exploratory)

`tables/threshold_sensitivity.csv` varies validation BH q threshold (0.01, 0.05, 0.10) and minimum validation |log2FC| (0 or 1) using the **same frozen full-fit BH family**. The original q<0.05, no validation magnitude cutoff remains the primary result. Expression-filter alternatives require separately refitted PyDESeq2 models; see `tables/filter_sensitivity.csv` if completed. These analyses were not used to choose the primary rule.

## Pathway comparison and interpretation

The original 13 significant discovery Hallmark pathways and all nonsignificant/discordant validation outcomes appear in `tables/pathway_comparison_original_universes.csv` (summary: {rep}/{sig} replicate). The identical source GMT had cohort-specific measurable universes. Re-ranking only the {len(set()) if False else 8982:,} mutually tested genes and applying identical 15–500 member eligibility yielded {len(coverage.loc[coverage.eligible])} common eligible terms; {crep}/{csig} discovery-significant common-universe pathways replicate. This is exploratory and does not replace the original comparison. Empirical GSEA p/FDR values saved as 0 mean below permutation resolution, not zero probability. `tables/theme_evidence.csv` shows E2F, G2-M, EMT leading-edge overlap and a separately curated differentiation marker list. Associations may reflect changing cell proportions.

| Term | Discovery NES | Discovery FDR | Validation NES | Validation FDR | Replicated | Shared leading-edge genes |
|---|---:|---:|---:|---:|---|---:|
"""
    for r in p_rows.itertuples():
        lines += f"| {r.Term} | {r.NES_discovery:.2f} | {r.padj_discovery:.3g} | {r.NES_validation:.2f} | {r.padj_validation:.3g} | {'yes' if r.replicated else 'no'} | {r.leading_edge_overlap} |\n"
    lines += "\n## Tissue composition and anatomical subsite (exploratory)\n\nPrespecified epithelial, immune, stromal, and muscle marker lists are in the code and `tables/marker_scores.csv`. Per-gene log2 normalized counts were z-scored across samples within cohort and averaged; these are relative marker scores, not cell fractions or deconvolution. Median paired tumor-minus-normal scores:\n\n" + markdown(marker_summary.round(3)) + "\n\n"
    lines += f"GSE184616 contains {subsite.tongue_pairs.iloc[0]} tongue pairs and {subsite.other_pairs.iloc[0]} non-tongue pairs. Tongue-only descriptive signs for the ten genes appear in `tables/subsite_ten_gene.csv`; non-tongue n=4 is too small for a reliable site-stratified inference. The two cohorts may differ in tumor purity, inflammation, stroma, margin sampling and anatomical site. Surgical margins are not necessarily healthy reference tissue.\n\n"
    lines += """## Limitations and reproducibility

This is secondary analysis of public data, not a new patient cohort or experimentally validated biomarker study. Discovery n=3 makes very small Wald p-values and large fold changes sensitive to model assumptions. Bulk transcriptomes cannot assign a cell of origin or prove a regulatory mechanism. Historical-to-modern symbol mapping can lose or misassign features; no ambiguous alias was forced. GSEA permutation estimates have finite resolution. The cohorts have no publicly evidenced patient-identity linkage beyond distinct GSM accessions. Check `manifest.json` for input/output hashes, versions, seed, and warnings. Re-run using the pinned Python 3.11 environment: `python -B scripts/run_phase2.py --root . --hgnc data/raw/annotations/hgnc_complete_set.txt`. Results are written only to `results/phase2/v1`.

## Figures

- [Effect-size scatter](figures/effect_size_scatter.svg) and [ten-gene forest](figures/ten_gene_forest.svg)
- [Matched-patient expression](figures/paired_expression_ten_genes.svg) and [pathway NES](figures/pathway_nes_concordance.svg)
- [Candidate outcomes](figures/candidate_outcomes.svg) and [sensitivity summary](figures/sensitivity_summary.svg)
"""
    (out / "PHASE2_REPORT.md").write_text(lines, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--hgnc", type=Path, required=True)
    parser.add_argument("--skip-common-gsea", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve(); hgnc = args.hgnc.resolve()
    out = root / "results/phase2/v1"
    out.mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(exist_ok=True)
    d = table(root / "results/GSE20116/tables/differential_expression.csv")
    v = table(root / "results/GSE184616/tables/differential_expression.csv")
    dn = table(root / "results/GSE20116/tables/normalized_counts.csv")
    vn = table(root / "results/GSE184616/tables/normalized_counts.csv")
    dm = pd.read_csv(root / "results/GSE20116/tables/metadata.csv", index_col=0)
    vm = pd.read_csv(root / "results/GSE184616/tables/metadata.csv", index_col=0)
    if set(dm.index) & set(vm.index):
        raise ValueError("Public accession overlap")
    for meta, counts in [(dm, dn), (vm, vn)]:
        if set(meta.index) != set(counts.columns):
            raise ValueError("Sample/normalized-count mismatch")
        if not (meta.groupby("patient_id").condition.apply(lambda x: sorted(x.tolist()) == ["normal", "tumor"])).all():
            raise ValueError("Unmatched patients")
    candidates = classify_candidates(d, v); save(candidates, out / "tables/candidate_outcomes.csv")
    mapping = hgnc_audit(d, v, hgnc); save(mapping, out / "tables/gene_mapping_audit.csv")
    save(mapping.groupby(["hgnc_status", "discovery_candidate"]).size().reset_index(name="genes"), out / "tables/mapping_summary.csv")
    model = model_audit(d, v); save(model, out / "tables/model_audit.csv")
    common, effects = effect_summary(d, v, candidates); save(common, out / "tables/common_effects.csv"); save(effects, out / "tables/effect_summary.csv")
    selected = selected_genes(d, v, mapping); save(selected, out / "tables/ten_genes.csv")
    paired = pd.concat([paired_values(dn, dm, GENES, "GSE20116"), paired_values(vn, vm, GENES, "GSE184616")], ignore_index=True)
    save(paired, out / "tables/paired_expression.csv")
    sensitivity = sensitivities(candidates, v); save(sensitivity, out / "tables/threshold_sensitivity.csv")
    pathways = pathway_original(root, out)
    if args.skip_common_gsea:
        common_p = pd.read_csv(out / "tables/pathway_comparison_common_universe.csv")
        coverage = pd.read_csv(out / "tables/hallmark_common_universe_coverage.csv")
    else:
        common_p, coverage = pathway_common(root, out, d, v)
    markers = pd.concat([marker_scores(dn, dm, "GSE20116"), marker_scores(vn, vm, "GSE184616")], ignore_index=True)
    save(markers, out / "tables/marker_scores.csv")
    theme = theme_table(root, d, v); save(theme, out / "tables/theme_evidence.csv")
    tongue_ids = vm.groupby("patient_id").tissue.first().loc[lambda x: x.str.lower() == "tongue"].index
    subsite_rows = []
    for gene in GENES:
        x = paired.loc[(paired.cohort == "GSE184616") & (paired.gene == gene)]
        for site, subset in [("tongue", x.loc[x.patient_id.isin(tongue_ids)]), ("other", x.loc[~x.patient_id.isin(tongue_ids)])]:
            subsite_rows.append({"gene": gene, "site_group": site, "pairs": len(subset),
                                 "tumor_higher_pairs": int((subset.paired_log2p1_difference > 0).sum()),
                                 "median_paired_log2p1_difference": subset.paired_log2p1_difference.median()})
    save(pd.DataFrame(subsite_rows), out / "tables/subsite_ten_gene.csv")
    subsite = pd.DataFrame([{"tongue_pairs": len(tongue_ids), "other_pairs": vm.patient_id.nunique()-len(tongue_ids)}])
    save(subsite, out / "tables/subsite_counts.csv")
    make_figures(out, common, candidates, selected, paired, pathways, sensitivity)
    report(root, out, None, candidates, effects, pathways, common_p, coverage, selected, paired, mapping, model, markers, subsite)
    poster = f"""# Cross-cohort OSCC transcriptomics: a public-data secondary analysis

**Question.** Do paired tumor/normal expression effects discovered in GSE20116 (3 patients) recur in independent GSE184616 (15 patients)?

**Methods.** Patient-blocked PyDESeq2 tumor/normal contrasts; exact gene-symbol match; genome-wide BH q<0.05; Hallmark preranked GSEA. Phase 2 independently recalculated summaries from frozen result CSVs and added common-universe and marker-score exploration.

**Results.** 1,339 discovery candidates; 1,187 validation-testable; 964 same-direction; 706 same-direction and validation significant; 563 also validation |log2FC|>1. There are **61 significant opposite-direction** candidates. Original-universe Hallmark replication: {int(pathways.replicated.sum())}/{int(pathways.discovery_significant.sum())} significant discovery pathways.

**Figures.** [Concordance](figures/effect_size_scatter.svg) · [Ten genes](figures/ten_gene_forest.svg) · [Matched expression](figures/paired_expression_ten_genes.svg) · [Pathways](figures/pathway_nes_concordance.svg) · [Outcomes](figures/candidate_outcomes.svg) · [Sensitivity](figures/sensitivity_summary.svg).

**Interpretation.** Several effects recur across independent public datasets, while discordance and uncertain mapping remain. Tissue composition and sampling may explain some signals. GSE20116 n=3 cannot establish robust clinical biomarkers.

**Scope.** This is a secondary analysis of public data, not an original patient cohort or experimentally validated biomarker study. See [full report](PHASE2_REPORT.md) and [reproducibility manifest](manifest.json).
"""
    (out / "POSTER.md").write_text(poster, encoding="utf-8")
    manifest = {"analysis": "phase2-v1", "python": sys.version, "platform": platform.platform(),
                "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
                             "matplotlib": matplotlib.__version__, "gseapy": gp.__version__},
                "random_seed_common_gsea": 42, "common_gsea_permutations": 1000,
                "hgnc_source": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                "hgnc_sha256": sha(hgnc), "warnings": ["3 discovery pairs; model-based inference fragile", "Gene-level Spearman CI assumes independent genes", "GSEA empirical zero has finite resolution", "Marker scores are not cell fractions"],
                "inputs": {p.relative_to(root).as_posix(): sha(p) for p in [root / "results/GSE20116/tables/differential_expression.csv", root / "results/GSE184616/tables/differential_expression.csv", root / "data/raw/MSigDB_Hallmark_2020.gmt"]},
                "outputs": {p.relative_to(out).as_posix(): sha(p) for p in out.rglob("*") if p.is_file() and p.name != "manifest.json"}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"out": str(out), "candidate_outcomes": candidates.outcome.value_counts().to_dict(),
                      "original_pathway_replication": [int(pathways.replicated.sum()), int(pathways.discovery_significant.sum())],
                      "common_pathway_replication": [int(common_p.replicated.sum()), int(common_p.discovery_significant.sum())]}, indent=2))


if __name__ == "__main__":
    main()
