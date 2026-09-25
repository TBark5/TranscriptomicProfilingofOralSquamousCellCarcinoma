"""Publication figures generated solely from computed analysis outputs."""
from pathlib import Path
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from .differential import classify

COLORS = {"normal": "#327D9D", "tumor": "#B74B44", "up": "#B74B44", "down": "#327D9D", "not_significant": "#B8BEC6"}


def save(fig, directory, name):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(directory / f"{name}.{extension}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def safe_neglog(p):
    return -np.log10(pd.Series(p).clip(lower=1e-300, upper=1))


def top_genes(results, n=30, alpha=.05, lfc=1.):
    ranked = results.loc[results.padj.notna()].sort_values(["padj", "pvalue"], kind="stable")
    hits = ranked.loc[classify(ranked, alpha, lfc) != "not_significant"]
    if len(hits) >= 2:
        return hits.head(n).index.tolist(), "Top threshold-passing genes"
    return ranked.head(n).index.tolist(), "Top ranked genes (exploratory; fewer than two pass thresholds)"


def make_figures(counts, log_expr, metadata, results, gsea, curve, output, alpha=.05, lfc=1., dataset="GSE20116"):
    sns.set_theme(style="ticks", context="paper", font="DejaVu Sans", font_scale=1.1)
    plt.rcParams.update({"svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
    output = Path(output)
    sample_labels = metadata.patient_id + " " + metadata.condition
    large_cohort = len(metadata) > 12
    patient_colors = dict(zip(metadata.patient_id.unique(), sns.color_palette("tab20" if large_cohort else "Set2", metadata.patient_id.nunique())))
    captions = []

    fig, axes = plt.subplots(1, 2, figsize=(20, 6) if large_cohort else (11, 4))
    axes[0].bar(range(len(metadata)), counts.sum(axis=0) / 1e6, color=metadata.condition.map(COLORS))
    axes[0].set(xticks=range(len(metadata)), xticklabels=sample_labels, ylabel="Count sum (millions)", title="Retained representative-count totals")
    sns.boxplot(data=log_expr.rename(columns=sample_labels.to_dict()), ax=axes[1], color="#DDE6EC", showfliers=False)
    axes[1].set(ylabel="log2(normalized count + 1)", title="Expression distributions")
    for ax in axes:
        ax.tick_params(axis="x", rotation=90 if large_cohort else 45, labelsize=7 if large_cohort else 9)
    fig.tight_layout()
    save(fig, output, "quality_control")
    captions.append("quality_control: Count sums across retained gene representatives are not full sequencing-library sizes. Boxplots show log2(size-factor-normalized counts + 1); outliers hidden for readability.")

    fig, ax = plt.subplots(figsize=(12, 10) if large_cohort else (6, 5))
    correlation = log_expr.corr()
    sns.heatmap(correlation, annot=not large_cohort, fmt=".2f", cmap="Blues", vmin=0, vmax=1, ax=ax,
                xticklabels=sample_labels, yticklabels=sample_labels, cbar_kws={"label": "Pearson r"})
    ax.set_title("Sample correlation · log normalized expression")
    save(fig, output, "sample_correlation")
    captions.append("sample_correlation: Pearson correlations across retained genes using log normalized expression, without patient or condition regression.")

    variable = log_expr.var(axis=1).nlargest(min(2000, len(log_expr))).index
    pca = PCA(n_components=2, svd_solver="full")
    scores = pca.fit_transform(log_expr.loc[variable].T)
    pca_table = pd.DataFrame(scores, index=metadata.index, columns=["PC1", "PC2"])
    fig, ax = plt.subplots(figsize=(11, 8) if large_cohort else (7, 5))
    for patient, group in metadata.groupby("patient_id", sort=False):
        positions = metadata.index.get_indexer(group.index)
        ax.plot(scores[positions, 0], scores[positions, 1], color=patient_colors[patient], alpha=.7, zorder=1)
    for i, (_, sample) in enumerate(metadata.iterrows()):
        ax.scatter(*scores[i], color=COLORS[sample.condition], s=80, marker="o" if sample.condition == "normal" else "^", zorder=2)
        label = f"{sample.patient_id.replace('OSCC_', '')} {'N' if sample.condition == 'normal' else 'T'}" if large_cohort else f"{sample.patient_id} {sample.condition}"
        ax.annotate(label, scores[i], xytext=(5, 5), textcoords="offset points", fontsize=8 if large_cohort else 9)
    ax.margins(.25)
    ax.set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", ylabel=f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", title="PCA · patient-matched samples")
    save(fig, output, "pca")
    captions.append("pca: Centered, unscaled PCA on up to 2,000 most variable log normalized genes; lines connect matched patients. Selection does not use DE p-values. PCA is descriptive, not a significance test.")

    fig, ax = plt.subplots(figsize=(8, 6))
    labels = classify(results, alpha, lfc)
    for group in ("not_significant", "down", "up"):
        subset = results.loc[labels == group]
        ax.scatter(subset.log2FoldChange, safe_neglog(subset.padj), s=9, alpha=.65, color=COLORS[group], label=f"{group.replace('_', ' ')} ({len(subset):,})", rasterized=False)
    ax.axhline(-np.log10(alpha), color="#555555", ls="--", lw=.8)
    for threshold in (-lfc, lfc):
        ax.axvline(threshold, color="#555555", ls="--", lw=.8)
    for direction, label_x in [("down", .025), ("up", .84)]:
        informative = results.loc[labels == direction].sort_values("padj").head(4)
        for i, (gene, row) in enumerate(informative.iterrows()):
            ax.annotate(gene, (row.log2FoldChange, safe_neglog([row.padj]).iloc[0]),
                        xytext=(label_x, .95 - .08 * i), textcoords="axes fraction", fontsize=9,
                        arrowprops={"arrowstyle": "-", "color": "#777777", "lw": .6})
    ax.set(xlabel="log2 fold change · tumor / normal", ylabel="−log10(BH-adjusted p-value)", title=f"Paired differential expression · {dataset}")
    ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(.56, 1.))
    save(fig, output, "volcano")
    captions.append(f"volcano: Paired PyDESeq2 Wald tests of tumor versus normal; thresholds BH-adjusted p < {alpha} and |log2 fold change| > {lfc}. Missing adjusted p-values are omitted; zeros are clipped to 1e-300 for plotting only. Labels select up to four actual threshold-passing genes per direction. Fold changes are unshrunk.")

    genes, heat_title = top_genes(results, alpha=alpha, lfc=lfc)
    if len(genes) >= 2:
        values = log_expr.loc[genes]
        z = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
        col_colors = pd.DataFrame({"Condition": metadata.condition.map(COLORS), "Patient": metadata.patient_id.map(patient_colors)}, index=metadata.index)
        grid = sns.clustermap(z, method="average", metric="euclidean", col_colors=col_colors,
                              cmap="vlag", center=0, vmin=-2, vmax=2, figsize=(14, 11) if large_cohort else (9, 10),
                              xticklabels=sample_labels, yticklabels=True, cbar_kws={"label": "Row z-score"})
        grid.fig.suptitle(heat_title, y=1.02)
        handles = [Line2D([0], [0], marker="s", ls="", color=color, label=str(name)) for name, color in {**{k: COLORS[k] for k in ("normal", "tumor")}, **patient_colors}.items()]
        grid.ax_col_dendrogram.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, 1.), ncol=5, frameon=False, fontsize=8)
        save(grid.fig, output, "heatmap")
    else:
        fig, ax = plt.subplots(figsize=(7, 4)); ax.axis("off")
        ax.text(.5, .5, "Fewer than two genes have valid adjusted p-values", ha="center")
        save(fig, output, "heatmap")
    captions.append("heatmap: Up to 30 top genes, ranked by adjusted p-value; row z-scores of log normalized expression, average-linkage Euclidean clustering on rows and columns. Color strips identify condition and patient. This DE-selected view is not independent evidence of separation; any non-significant fallback is identified in the title.")

    candidates = [gene for direction in ["up", "down"] for gene in results.loc[labels == direction].sort_values("padj").head(2).index]
    if not candidates:
        candidates = genes[:4]
    fig, axes = plt.subplots(1, max(1, len(candidates)), figsize=(max(6, 3 * len(candidates)), 7 if large_cohort else 4), squeeze=False)
    for ax, gene in zip(axes[0], candidates):
        for patient, group in metadata.groupby("patient_id", sort=False):
            ids = [group.index[group.condition == c][0] for c in ("normal", "tumor")]
            ax.plot([0, 1], log_expr.loc[gene, ids], marker="o", color=patient_colors[patient], label=patient)
        row = results.loc[gene]
        ax.set(xticks=[0, 1], xticklabels=["Normal", "Tumor"], title=f"{gene}\nLFC={row.log2FoldChange:.2f}; q={row.padj:.2g}")
    if candidates:
        axes[0, 0].set_ylabel("log2(normalized count + 1)")
        axes[0, -1].legend(title="Patient", frameon=False)
    else:
        axes[0, 0].text(.5, .5, "No eligible candidate genes", ha="center")
    fig.tight_layout()
    save(fig, output, "paired_expression")
    captions.append("paired_expression: Individual patient trajectories for up to four data-selected candidates. Points are log normalized expression, while displayed LFC and q are from the paired count model.")

    fig, ax = plt.subplots(figsize=(9, 6))
    if not gsea.empty:
        table = gsea.copy()
        for col in ["NES", "padj"]:
            table[col] = pd.to_numeric(table[col])
        table = table.sort_values("padj").head(15).sort_values("NES")
        y = np.arange(len(table))
        scatter = ax.scatter(table.NES, y, c=safe_neglog(table.padj).clip(upper=10), cmap="viridis", s=65)
        ax.set_yticks(y, [textwrap.fill(t, 37) for t in table.Term])
        ax.axvline(0, color="grey", lw=.7)
        fig.colorbar(scatter, ax=ax, label="−log10(GSEA FDR), capped at 10")
        ax.set(xlabel="Normalized enrichment score (NES)", title="Top preranked pathways · exploratory")
    else:
        ax.text(.5, .5, "No pathway results available; see enrichment_status.json", ha="center", transform=ax.transAxes)
    save(fig, output, "pathways")
    captions.append("pathways: Up to 15 lowest-FDR preranked GSEA terms, including non-significant terms if necessary. Positive NES indicates enrichment toward positive tumor/normal Wald statistics. FDR is GSEA's permutation-derived estimate, not BH. Gene-set enrichment does not establish activation or causality.")
    if curve:
        fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
        axes[0].plot(curve["running_es"], color="#327D9D")
        axes[0].axhline(0, color="grey", lw=.7)
        axes[0].set(ylabel="Running enrichment score", title=f"{curve['term']}\nNES={curve['nes']:.2f}; GSEA FDR={curve['fdr']:.3g}")
        axes[1].vlines(curve["hits"], 0, 1, color="#555555", lw=.5)
        axes[1].set(xlabel="Gene rank · positive to negative tumor/normal Wald statistic", yticks=[])
        save(fig, output, "enrichment_curve")
        captions.append("enrichment_curve: Running weighted enrichment score and gene-set hit positions for the first library's lowest-FDR term; nominal zero permutation p-values/FDR are resolution-limited, not certainty.")
    else:
        for extension in ("png", "svg"):
            (output / f"enrichment_curve.{extension}").unlink(missing_ok=True)
    (output / "CAPTIONS.md").write_text("# Figure captions\n\n" + "\n\n".join(captions) + "\n", encoding="utf-8")
    return pca_table, pca.explained_variance_ratio_.tolist(), correlation
