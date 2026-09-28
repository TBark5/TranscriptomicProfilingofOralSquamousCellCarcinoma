"""Write a research report from measured results, never template findings."""
from pathlib import Path
import pandas as pd
from .differential import classify


def markdown_table(frame):
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in row) + " |")
    return "\n".join(lines)


def write_report(root, manifest, result, gsea, qc, ora=None, report_path=None):
    if manifest["dataset"] != "GSE20116":
        raise ValueError("Discovery report requires GSE20116; use the independent validation report writer")
    config = manifest["config"]
    calls = classify(result, config["alpha"], config["lfc"])
    leading = pd.concat([result.loc[calls == d].sort_values("padj").head(5) for d in ["up", "down"]])
    gene_table = markdown_table(leading.reset_index()[["gene", "refseq", "log2FoldChange", "pvalue", "padj"]]) if len(leading) else "No genes pass both configured cutoffs."
    pathway_table = "No GSEA results were produced; consult results/tables/enrichment_status.json."
    if not gsea.empty:
        display = gsea.copy()
        for column in ["NES", "pvalue", "padj"]:
            display[column] = pd.to_numeric(display[column])
        display = display.sort_values("padj").head(8)
        pathway_table = markdown_table(display[["library", "Term", "NES", "pvalue", "padj"]])
    ora_table = "No eligible ORA result tables were produced."
    if ora is not None and not ora.empty:
        ora_table = markdown_table(ora.sort_values("padj").head(6)[["direction", "Term", "overlap", "fold_enrichment", "padj"]])
    comparison = []
    for symbol in ["MMP1", "INHBA", "HMGA2", "CASQ1"]:
        if symbol in result.index:
            row = result.loc[symbol]
            comparison.append(f"{symbol}: LFC {row.log2FoldChange:.2f}, BH q {row.padj:.3g}, current call `{calls.loc[symbol]}`")
    comparison_text = "; ".join(comparison)
    text = f"""# Transcriptomic Profiling of Oral Squamous Cell Carcinoma

## Background and hypothesis
This exploratory reanalysis asks which gene representatives and pathways differ between OSCC and matched normal oral tissue. The hypothesis is that patient-adjusted expression differences are present. Differential expression measures association; it does not identify causal drivers or establish clinical utility.

## Dataset and provenance
GSE20116, GPL9442, SRP002009, PRJNA124251; three matched pairs (patients 8, 33, 51), six bulk RNA-seq libraries. The actual GEO SOFT records supply sample identities and describe hg18 alignments in MAX format. Instead of downloading the approximately 3.9 GB alignment archive, this project reads raw count columns D–I from Tuch et al. Table S1 (supplement s009). These raw sums are distinct from the adjacent normalized columns and from the publication's statistical results. No published p-values are reused. Source URLs and SHA256 hashes are in `results/manifest.json`; all sample identities and available GEO characteristics are in `results/tables/metadata.csv`.

Counts were originally quantified from uniquely aligned reads over RefSeq exons using the AB SOLiD whole-transcriptome pipeline and hg18. The paper describes a prefilter of at least 50 reads in either tissue in each patient. Actual workbook inspection finds 47 transcript rows with a pair maximum of 49; all rows meet 49 in every pair. The original measurements are preserved; the boundary discrepancy may relate to the source pseudocount but its cause is not established. This is an ascertainment-limited universe, not all human genes. This project selects the transcript with most exons per source gene symbol, breaking ties lexically by RefSeq ID, and never sums overlapping isoform counts. These are **gene representatives**, not modern gene-union quantifications. Source symbols remain historical; exact gene-set matching, rejected mappings and redundant isoforms are recorded. No current-symbol alias guesses are made. Age, sex, stage and HPV status are unavailable in the GEO records and remain missing. The paper describes normal specimens collected at negative surgical margins [1].

## Methods
Validation checks finite nonnegative integer counts, unique feature/sample identities, exact metadata alignment, complete pairs and full-rank design. Of {manifest['source_transcripts']:,} input transcripts, {manifest['representatives']:,} gene representatives were selected and {len(result):,} passed an additional filter of at least {config['min_count']} reads in {config['min_samples']} samples. That abundance-only filter does not use tumor/normal fold changes.

PyDESeq2 0.5.2 fits negative-binomial GLMs with design `~ patient_id + condition`, contrast tumor/normal, and median-of-ratios size factors [2,3]. Raw integer counts enter the model. Normalized counts and log2(normalized counts + 1) are stored separately; the latter is a visualization transform, not a variance-stabilizing transformation. The design rank is {manifest['diagnostics']['design_rank']} with {manifest['diagnostics']['residual_df']} residual degrees of freedom. Cook filtering is enabled; count replacement is disabled for this small paired design. Independent filtering is disabled to keep the BH family explicit. P-values for any failed gene-wise, MAP or LFC optimization are excluded before BH correction; the complete audit retains their estimates and reasons. Wald tests assess a zero log fold change, with an additional descriptive effect-size cutoff (not a formal test against |LFC| = 1). LFCs are unshrunk and can be unstable.

Default calls require adjusted p < {config['alpha']} and |log2 fold change| > {config['lfc']}. Preranked GSEA uses all finite, convergence-qualified, Cook-filtered Wald statistics, with {config['permutations']} gene-set permutations and seed {config['seed']}. Gene sets are restricted to 15–500 tested members. Overrepresentation uses separate up/down lists of at least five genes and the tested background; zero-overlap eligible terms remain in the hypergeometric BH family. ORA BH is within each library/direction; GSEA FDR is the permutation-derived estimate within each library. These are distinct corrections, with no global adjustment across exploratory collections. GMT sources and overlap audits are recorded [4,5].

## Quality control
All six libraries passed structural validation. Retained representative-count totals ranged from {qc['library_min']:,} to {qc['library_max']:,}; these totals do not equal total sequenced reads. Sample correlations ranged from {qc['correlation_min']:.3f} to {qc['correlation_max']:.3f} on log normalized expression. PCA uses the 2,000 most variable genes (or all if fewer), centered without gene scaling: PC1 explained {qc['pca_variance'][0]:.1%} and PC2 {qc['pca_variance'][1]:.1%}. No sample was excluded based on visual appearance. Distribution, correlation, PCA and matched-expression figures accompany this report. Raw-read quality and alignment quality were not recomputed from the historical reads.

## Differential-expression results
{int((calls == 'up').sum()):,} representatives were upregulated and {int((calls == 'down').sum()):,} downregulated at the configured cutoffs. {manifest['diagnostics']['failed_convergence']} features failed at least one convergence check; {manifest['diagnostics']['missing_pvalue']} features have no eligible p-value in total. Leading candidates below are selected only from actual threshold-passing results, ordered by adjusted p-value:

{gene_table}

The candidate table includes up to five genes per direction. The volcano plot uses computed p-values only. The heatmap is selected using the DE results and hence is not independent confirmation of group separation. Patient trajectories expose heterogeneity hidden by the common modeled condition effect. Complete results include estimates for non-significant and excluded features.

## Pathway results
Lowest-FDR terms are shown below regardless of significance; inspect the displayed FDR before interpreting a term as enriched. Positive NES means enrichment toward tumor-upregulated statistics; negative NES means enrichment toward normal-upregulated statistics. Neither implies pathway activation or inhibition.

{pathway_table}

Leading overrepresentation terms (separate up/down queries; BH within each query/library):

{ora_table}

Nominal zero permutation p-values or FDR estimates reflect finite resolution, not a probability proven to be zero. Preranking uses gene-set rather than patient-label permutations and does not account fully for gene correlation. All pathways and leading-edge genes are saved in `gsea.csv`; ORA overlaps are in `ora.csv`.

## Discussion and limitations
Tuch et al. reported expression differences involving matrix remodeling and differentiation [1]. For the four expression examples highlighted in that paper, the current model yields: {comparison_text}. Agreement with that publication is a same-cohort comparison, **not external replication**. Historical candidate lists are not evidence of significance in the current fit. Cell-cycle and muscle-related enrichment, if supported by the tables above, describes patterns in the ranked measurements; loss of muscle-rich normal tissue is a plausible compositional explanation that this design cannot distinguish from tumor-intrinsic regulation.

Three patients provide limited power and only two residual degrees of freedom. Strong observed effects do not remove uncertainty about broader OSCC populations. Bulk tissue composition, stromal/immune admixture, muscle content and tumor heterogeneity can produce expression differences. Normal surgical margins may contain field effects and differ in cell composition. Patient adjustment preserves pairing but cannot separate condition from an unrecorded condition-confounded technical batch. There is no basis here for survival prediction, diagnosis or treatment recommendations.

Legacy SOLiD chemistry, color-space alignment, hg18 annotation, source prefiltering and historical symbols limit comparability with modern RNA-seq. A single representative can miss isoform changes; overlapping genomic genes may still share reads. Size-factor normalization assumes most genes are unchanged or changes balance sufficiently. The p-values rely on negative-binomial and asymptotic Wald assumptions in a very small cohort. Functional enrichment also depends on gene-set definitions, mapping coverage, unequal gene detectability and the restricted background.

## Follow-up
Validate candidate direction and effect sizes in an independent, anatomically restricted OSCC cohort with known pairs and sufficient patients. Pre-specify candidate genes, quantify tumor purity and relevant covariates, and test with qPCR using validated reference genes and biological replicates. Mechanistic claims require appropriately controlled perturbation experiments and replication. No external cohort or laboratory validation was performed here.

## Reproducibility and references
Run `oscc run --offline` after caching the documented sources. Configuration, installed versions, source and output checksums, model warnings and enrichment status are recorded. See README for installation and the notebook for a walkthrough.

1. Tuch BB et al. (2010). [Tumor Transcriptome Sequencing Reveals Allelic Expression Imbalances Associated with Copy Number Alterations](https://doi.org/10.1371/journal.pone.0009317). PLOS ONE. [GEO GSE20116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20116).
2. Love MI, Huber W, Anders S (2014). [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8).
3. Muzellec B et al. (2023). [PyDESeq2: a python package for bulk RNA-seq differential expression analysis](https://doi.org/10.1093/bioinformatics/btad547).
4. Subramanian A et al. (2005). [Gene set enrichment analysis](https://doi.org/10.1073/pnas.0506580102); Liberzon A et al. (2015). [The Molecular Signatures Database hallmark gene set collection](https://doi.org/10.1016/j.cels.2015.12.004).
5. Fang Z et al. (2023). [GSEApy: a comprehensive package for performing gene set enrichment analysis in Python](https://doi.org/10.1093/bioinformatics/btac757). [Enrichr gene-set libraries](https://maayanlab.cloud/Enrichr/#libraries).
"""
    destination = Path(report_path) if report_path else Path(root, "results/GSE20116/REPORT.md")
    if destination.resolve() == Path(root, "REPORT.md").resolve():
        raise ValueError("The original REPORT.md is preserved; choose a dataset output report")
    destination.parent.mkdir(parents=True, exist_ok=True)
    text = text.replace("results/tables/", "tables/").replace("results/manifest.json", "manifest.json")
    destination.write_text(text, encoding="utf-8")
