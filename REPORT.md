# Transcriptomic Profiling of Oral Squamous Cell Carcinoma

## Background and hypothesis
This exploratory reanalysis asks which gene representatives and pathways differ between OSCC and matched normal oral tissue. The hypothesis is that patient-adjusted expression differences are present. Differential expression measures association; it does not identify causal drivers or establish clinical utility.

## Dataset and provenance
GSE20116, GPL9442, SRP002009, PRJNA124251; three matched pairs (patients 8, 33, 51), six bulk RNA-seq libraries. The actual GEO SOFT records supply sample identities and describe hg18 alignments in MAX format. Instead of downloading the approximately 3.9 GB alignment archive, this project reads raw count columns D–I from Tuch et al. Table S1 (supplement s009). These raw sums are distinct from the adjacent normalized columns and from the publication's statistical results. No published p-values are reused. Source URLs and SHA256 hashes are in `results/manifest.json`; all sample identities and available GEO characteristics are in `results/tables/metadata.csv`.

Counts were originally quantified from uniquely aligned reads over RefSeq exons using the AB SOLiD whole-transcriptome pipeline and hg18. The paper describes a prefilter of at least 50 reads in either tissue in each patient. Actual workbook inspection finds 47 transcript rows with a pair maximum of 49; all rows meet 49 in every pair. The original measurements are preserved; the boundary discrepancy may relate to the source pseudocount but its cause is not established. This is an ascertainment-limited universe, not all human genes. This project selects the transcript with most exons per source gene symbol, breaking ties lexically by RefSeq ID, and never sums overlapping isoform counts. These are **gene representatives**, not modern gene-union quantifications. Source symbols remain historical; exact gene-set matching, rejected mappings and redundant isoforms are recorded. No current-symbol alias guesses are made. Age, sex, stage and HPV status are unavailable in the GEO records and remain missing. The paper describes normal specimens collected at negative surgical margins [1].

## Methods
Validation checks finite nonnegative integer counts, unique feature/sample identities, exact metadata alignment, complete pairs and full-rank design. Of 15,668 input transcripts, 10,541 gene representatives were selected and 10,541 passed an additional filter of at least 10 reads in 3 samples. That abundance-only filter does not use tumor/normal fold changes.

PyDESeq2 0.5.2 fits negative-binomial GLMs with design `~ patient_id + condition`, contrast tumor/normal, and median-of-ratios size factors [2,3]. Raw integer counts enter the model. Normalized counts and log2(normalized counts + 1) are stored separately; the latter is a visualization transform, not a variance-stabilizing transformation. The design rank is 4 with 2 residual degrees of freedom. Cook filtering is enabled; count replacement is disabled for this small paired design. Independent filtering is disabled to keep the BH family explicit. P-values for any failed gene-wise, MAP or LFC optimization are excluded before BH correction; the complete audit retains their estimates and reasons. Wald tests assess a zero log fold change, with an additional descriptive effect-size cutoff (not a formal test against |LFC| = 1). LFCs are unshrunk and can be unstable.

Default calls require adjusted p < 0.05 and |log2 fold change| > 1.0. Preranked GSEA uses all finite, convergence-qualified, Cook-filtered Wald statistics, with 1000 gene-set permutations and seed 42. Gene sets are restricted to 15–500 tested members. Overrepresentation uses separate up/down lists of at least five genes and the tested background; zero-overlap eligible terms remain in the hypergeometric BH family. ORA BH is within each library/direction; GSEA FDR is the permutation-derived estimate within each library. These are distinct corrections, with no global adjustment across exploratory collections. GMT sources and overlap audits are recorded [4,5].

## Quality control
All six libraries passed structural validation. Retained representative-count totals ranged from 7,365,969 to 21,502,237; these totals do not equal total sequenced reads. Sample correlations ranged from 0.613 to 0.897 on log normalized expression. PCA uses the 2,000 most variable genes (or all if fewer), centered without gene scaling: PC1 explained 59.9% and PC2 16.7%. No sample was excluded based on visual appearance. Distribution, correlation, PCA and matched-expression figures accompany this report. Raw-read quality and alignment quality were not recomputed from the historical reads.

## Differential-expression results
363 representatives were upregulated and 976 downregulated at the configured cutoffs. 35 features failed at least one convergence check; 35 features have no eligible p-value in total. Leading candidates below are selected only from actual threshold-passing results, ordered by adjusted p-value:

| gene | refseq | log2FoldChange | pvalue | padj |
| --- | --- | --- | --- | --- |
| PTHLH | NM_198965 | 3.914 | 1.065e-21 | 7.992e-19 |
| LAMC2 | NM_005562 | 3.787 | 2.084e-21 | 1.288e-18 |
| COL4A6 | NM_001847 | 3.723 | 7.296e-20 | 3.194e-17 |
| MMP11 | NM_005940 | 4.502 | 5.826e-18 | 1.53e-15 |
| SPP1 | NM_001040058 | 5.225 | 2.018e-17 | 4.511e-15 |
| TMPRSS11B | NM_182502 | -7.451 | 4.374e-31 | 4.596e-27 |
| PTGFR | NM_001039585 | -5.199 | 2.695e-26 | 1.415e-22 |
| PYGM | NM_005609 | -5.49 | 4.092e-26 | 1.433e-22 |
| CRNN | NM_016190 | -7.304 | 5.997e-26 | 1.575e-22 |
| MAL | NM_002371 | -6.965 | 6.647e-25 | 1.397e-21 |

The candidate table includes up to five genes per direction. The volcano plot uses computed p-values only. The heatmap is selected using the DE results and hence is not independent confirmation of group separation. Patient trajectories expose heterogeneity hidden by the common modeled condition effect. Complete results include estimates for non-significant and excluded features.

## Pathway results
Lowest-FDR terms are shown below regardless of significance; inspect the displayed FDR before interpreting a term as enriched. Positive NES means enrichment toward tumor-upregulated statistics; negative NES means enrichment toward normal-upregulated statistics. Neither implies pathway activation or inhibition.

| library | Term | NES | pvalue | padj |
| --- | --- | --- | --- | --- |
| MSigDB_Hallmark_2020 | E2F Targets | 2.605 | 0 | 0 |
| MSigDB_Hallmark_2020 | Myogenesis | -2.565 | 0 | 0 |
| MSigDB_Hallmark_2020 | Myc Targets V1 | 2.358 | 0 | 0 |
| MSigDB_Hallmark_2020 | mTORC1 Signaling | 2.169 | 0 | 0 |
| MSigDB_Hallmark_2020 | G2-M Checkpoint | 2.111 | 0 | 0 |
| MSigDB_Hallmark_2020 | Myc Targets V2 | 1.88 | 0 | 0.00132 |
| MSigDB_Hallmark_2020 | Glycolysis | 1.775 | 0 | 0.003299 |
| MSigDB_Hallmark_2020 | Adipogenesis | -1.773 | 0 | 0.003941 |

Leading overrepresentation terms (separate up/down queries; BH within each query/library):

| direction | Term | overlap | fold_enrichment | padj |
| --- | --- | --- | --- | --- |
| down | Myogenesis | 67 | 4.776 | 1.594e-28 |
| up | Epithelial Mesenchymal Transition | 23 | 4.16 | 2.85e-07 |
| up | Glycolysis | 20 | 3.911 | 4.075e-06 |
| up | Estrogen Response Late | 18 | 3.748 | 2.147e-05 |
| down | KRAS Signaling Up | 28 | 2.45 | 0.0001457 |
| up | KRAS Signaling Up | 15 | 3.53 | 0.0002575 |

Nominal zero permutation p-values or FDR estimates reflect finite resolution, not a probability proven to be zero. Preranking uses gene-set rather than patient-label permutations and does not account fully for gene correlation. All pathways and leading-edge genes are saved in `gsea.csv`; ORA overlaps are in `ora.csv`.

## Discussion and limitations
Tuch et al. reported expression differences involving matrix remodeling and differentiation [1]. For the four expression examples highlighted in that paper, the current model yields: MMP1: LFC 3.61, BH q 5.17e-09, current call `up`; INHBA: LFC 3.63, BH q 3.46e-09, current call `up`; HMGA2: LFC 3.83, BH q 2.32e-05, current call `up`; CASQ1: LFC -7.86, BH q 0.000302, current call `down`. Agreement with that publication is a same-cohort comparison, **not external replication**. Historical candidate lists are not evidence of significance in the current fit. Cell-cycle and muscle-related enrichment, if supported by the tables above, describes patterns in the ranked measurements; loss of muscle-rich normal tissue is a plausible compositional explanation that this design cannot distinguish from tumor-intrinsic regulation.

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
