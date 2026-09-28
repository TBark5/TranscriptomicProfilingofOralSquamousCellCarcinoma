"""Measured external-validation report, separate from the historical REPORT.md."""
from pathlib import Path
from .report import markdown_table


def write_validation_report(root, manifest, comparison, pathways, metadata):
    summary = manifest["summary"]
    diagnostics = manifest["diagnostics"]
    candidate = comparison.loc[comparison.discovery_candidate]
    featured = comparison.reindex(["PTHLH", "LAMC2", "COL4A6", "TMPRSS11B", "PTGFR", "PYGM"])
    columns = ["gene", "discovery_log2FoldChange", "validation_log2FoldChange", "validation_padj", "validation_outcome"]
    featured_table = markdown_table(featured.reset_index()[columns])
    outcomes = candidate.validation_outcome.value_counts().rename_axis("Outcome").reset_index(name="Count")
    mapping = candidate.mapping_status.value_counts().rename_axis("Mapping / expression status").reset_index(name="Count")
    pathway_columns = ["Term", "discovery_NES", "validation_NES", "discovery_padj", "validation_padj", "replicated"]
    pathway_table = markdown_table(pathways.loc[pathways.discovery_significant].sort_values("discovery_padj")[pathway_columns])
    patients = metadata.drop_duplicates("patient_id")
    sites = markdown_table(patients.tissue.value_counts().rename_axis("Anatomic site").reset_index(name="Pairs"))
    source_table = markdown_table(__import__("pandas").DataFrame([
        {"File": name, "SHA256": entry["sha256"], "Bytes": entry["bytes"]}
        for name, entry in manifest["sources"].items() if isinstance(entry, dict)]))
    concordance = summary["concordant_candidates"] / summary["testable_candidates"] if summary["testable_candidates"] else float("nan")
    report = f"""# Independent OSCC validation: GSE184616

## Execution and scope
The independent analysis completed using real deposited counts. GSE20116 remains the discovery cohort and GSE184616 is the external validation cohort. No sample pooling or joint normalization was performed. The original REPORT.md, discovery tables, figures, processed data and manifest were checked by SHA256 before and after validation and remained unchanged.

This is replication of differential-expression associations, not validation of a diagnostic classifier, prognosis model, causal mechanism or clinical biomarker. Sample accessions do not overlap. Separate studies, institutions and collection periods support cohort independence; individual identities cannot be audited from public anonymized records.

## Verified source and cohort
Source: [GEO GSE184616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184616), not a versioned accession ending in `.1`. We inspected `GSE184616_unnormalisedGeneCounts.txt.gz` and `GSE184616_family.soft.gz`. The deposited FPKM file was not used for count modeling. The count matrix has {manifest['source_genes']:,} unique gene-symbol rows and {manifest['samples']} samples, with finite, nonnegative integer values and exact GEO-title-to-count-column alignment. No rounding, imputation, or duplicate summation was needed.

The {manifest['pairs']} matched patients have IDs OSCC_1 through OSCC_14 and OSCC_16; no OSCC_15 is deposited. Each has one adjacent-normal oral mucosa and one primary tumor. The metadata-derived age range is {patients.age.min()}–{patients.age.max()} years; {int((patients.sex == 'Male').sum())} male and {int((patients.sex == 'Female').sum())} female patients. All samples were retained. Sites are heterogeneous:

{sites}

The series describes an HPV-negative cohort. Sample GSM5593760 (OSCC_7-P) gives only “Oral Squamous Cell Carcinoma” in its diagnosis field, while its normal counterpart explicitly says HPV-negative. The sample-level omission is retained as “not stated in sample”, not silently imputed. Age, sex, smoking and site agree within every pair. Patient blocking absorbs patient-level covariates; they are not added as collinear fixed effects.

GEO describes NovaSeq 6000, hg38, STAR 2.7.2 and GENCODE 31. The raw count file records strand-specific read-pair counts; RSEM-derived FPKMs are a separate resource. Raw sequencing/alignment quality was not recomputed. [Sample processing record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5593748).

{source_table}

## Analysis fixed before validation model fitting
Discovery candidates were frozen from the original saved discovery results using their saved thresholds (BH q < {manifest['candidate_selection']['alpha']}, absolute log2 fold change > {manifest['candidate_selection']['lfc']}). The {summary['discovery_candidates']:,}-gene candidate list and its checksum were saved before fitting the validation model. This is a retrospective external validation plan, not a prospectively registered study.

For validation, an abundance-only rule requires at least 10 reads in at least 15 of 30 samples, fixed at the size of one condition group without using DE outcomes. {summary['genes']:,} genes passed. PyDESeq2 {manifest['versions']['pydeseq2']} fits `~ patient_id + condition` with tumor/normal contrast, median-of-ratios normalization, no count replacement, Cook filtering, and no independent filtering. Design rank: {diagnostics['design_rank']}; residual degrees of freedom: {diagnostics['residual_df']}. {diagnostics['failed_convergence']} genes failed at least one optimizer convergence check. {summary['tested']:,} genes have finite eligible p-values; BH correction is across these validation tests, not just selected candidates. Model warnings are preserved in the manifest and diagnostics file.

Primary gene replication requires a frozen discovery candidate, the same effect direction, and validation genome-wide BH q < 0.05. A stricter descriptive count additionally requires validation absolute log2 fold change > 1. The Wald null is zero effect; the effect-size cutoff is not a formal test against a twofold-change boundary. No thresholds were tuned to maximize replication.

## Identifier coverage and candidate concordance
Both deposited tables already use gene symbols. Mapping is exact and case-sensitive, with no speculative alias conversion. Historical RefSeq representatives in GSE20116 are compared with modern gene-level counts in GSE184616; these are different quantification units. Unmatched, filtered, or failed-test genes remain visible in the candidate table and are not counted as evidence of biological non-replication.

{markdown_table(mapping)}

Of {summary['discovery_candidates']:,} frozen candidates, {summary['testable_candidates']:,} were testable and {summary['concordant_candidates']:,} had concordant signs ({concordance:.1%} of testable candidates). **{summary['replicated_candidates']:,} replicated** under the primary criterion; **{summary['replicated_large_effect']:,}** also exceeded absolute validation log2 fold change 1.

{markdown_table(outcomes)}

The following six candidates were highlighted in the original discovery README before external validation; they are displayed regardless of validation outcome:

{featured_table}

Across all {summary['common_tested_genes']:,} mutually testable genes, the descriptive Spearman correlation between discovery and validation log2 fold changes is **{summary['lfc_spearman']:.4f}**. This is not a classification accuracy or a formal test of equality of effects. The effect comparison table contains both standard errors and the difference between cohort estimates.

Validation-wide DE calls (q < 0.05 and absolute LFC > 1): {summary['up']:,} upregulated and {summary['down']:,} downregulated genes. These are distinct from candidate replication counts.

![Cohort effect comparison](results/GSE184616/figures/effect_size_comparison.png)

![Previously highlighted candidates](results/GSE184616/figures/candidate_effects.png)

Intervals are approximate unshrunk 95% Wald intervals. Discovery selection creates winner's-curse bias, and cross-platform differences limit literal equality of effect sizes.

## Pathway replication
Both analyses use the same checksum-locked human MSigDB_Hallmark_2020 definitions. Validation GSEA uses all finite, convergence-qualified signed Wald statistics, {manifest['config']['permutations']:,} gene-set permutations, seed {manifest['config']['seed']}, and 15–500 tested members per set. A pathway replicates when it has discovery and validation GSEA FDR < 0.05 and the same NES sign. This criterion refers to permutation-derived within-library FDR, not gene-level BH.

Of {summary['discovery_significant_pathways']} discovery-significant pathways, {summary['testable_discovery_pathways']} were testable in both cohorts and **{summary['replicated_pathways']} replicated**. Every discovery-significant term is shown, including discordant and nonsignificant terms:

{pathway_table}

![Pathway replication](results/GSE184616/figures/pathway_replication.png)

Background gene coverage differs between cohorts; NES values are not measurements on an identical scale. This comparison does not restrict both analyses to an identical common gene universe. Gene-set permutations do not fully model gene correlation, overlapping pathways are dependent, and zero FDR estimates reflect finite permutation resolution. Separate up/down ORA tables are exploratory supporting outputs, not the primary pathway replication criterion. Enrichment direction is not evidence of pathway activation or inhibition.

## Quality-control figures
![Validation PCA](results/GSE184616/figures/pca.png)

![Validation differential expression](results/GSE184616/figures/volcano.png)

PCA uses the most variable genes without DE selection. The saved DE heatmap and top-gene patient trajectories are descriptive and selected from validation DE results; they are not independent confirmation of those same results. The fixed six-gene comparison above uses discovery-selected candidates. All PNG figures also have SVG counterparts.

## Limitations
GSE20116 has only three patients and two residual degrees of freedom; its selected effects may be unstable. Legacy SOLiD/hg18 RefSeq-representative counts and prefiltering differ from NovaSeq/hg38 gene-level quantification. Exact symbol matching misses renamed genes. The 15-patient external series covers several oral subsites and an age-restricted, series-described HPV-negative population, limiting generalization. Bulk tissue composition, muscle content, stromal/immune admixture, tumor purity and normal-margin field effects can drive concordant signals without tumor-intrinsic regulation. No purity adjustment, anatomic-subsite sensitivity analysis, survival analysis, qPCR or perturbation experiments were performed. Patient-level clinical characteristics cannot remove condition-confounded technical effects. Replication supports associations in this external cohort, not clinical deployment.

## Reproduction and outputs
Run `.venv/Scripts/python.exe -m oscc.cli external-validation --offline` after source caching, or omit `--offline` to acquire checksum-locked missing inputs. Validation settings are fixed except seed and permutation count. Use `python -m oscc.cli validate --dataset GSE184616` for output-integrity and alignment checks. Do not run concurrent writers for the same dataset.

The original discovery snapshot remains under `results/tables`, `results/figures`, and `REPORT.md`. New discovery reruns write to `results/GSE20116`; external validation writes to `results/GSE184616` and this report. The comparison always reads the frozen original discovery snapshot, even if a newer discovery rerun exists.

- [Validation manifest](results/GSE184616/manifest.json): source/output hashes, frozen candidate hash, software versions, warnings and original-snapshot checks.
- [Candidate concordance](results/GSE184616/tables/candidate_concordance.csv): all frozen candidates, coverage and outcomes.
- [All discovery gene comparisons](results/GSE184616/tables/gene_concordance.csv).
- [Pathway replication](results/GSE184616/tables/pathway_replication.csv).
- [Validation differential expression](results/GSE184616/tables/differential_expression.csv).
- [Verified sample metadata](results/GSE184616/tables/metadata.csv).

Study reference: [Satgunaseelan et al., 2021, Oral Squamous Cell Carcinoma in Young Patients Show Higher Rates of EGFR Amplification](https://doi.org/10.3389/fonc.2021.750852).
"""
    Path(root, "VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
