# Output data dictionary

All research measurements here are derived from GSE20116 Table S1; none are synthetic. Source raw transcript rows are reduced to gene representatives as documented in `data/PROVENANCE.md`.

| File | Meaning |
|---|---|
| differential_expression.csv | One row per retained source-symbol representative. `baseMean`: mean size-factor-normalized count; `log2FoldChange`: unshrunk tumor/normal LFC; `lfcSE`: SE on log2 scale; `stat`: signed Wald statistic; `pvalue_model`: PyDESeq2 p after its Cook filter but before this project's convergence exclusion; `pvalue`: convergence-qualified p; `padj`: BH over finite p-values; `model_converged`, `test_status`: inference eligibility; `refseq`, `exons`: selected annotation; `direction`: saved cutoff call. |
| upregulated.csv, downregulated.csv | Strict saved cutoffs, q < alpha and absolute LFC > threshold; missing q never passes. |
| counts.csv | Genes × samples, raw nonnegative integer sums for retained representatives. |
| normalized_counts.csv | Same dimensions, raw counts divided by fitted sample size factors. Not TPM/FPKM. |
| log_expression.csv | log2(normalized count + 1), visualization only. |
| metadata.csv | GEO accession, patient, condition, workbook column, tissue, platform, genome and available clinical fields. Empty clinical fields mean unavailable, not normal/negative. |
| transcript_mapping.csv | Every source transcript, count/annotation fields and representative/redundant/ambiguous decision. |
| filter_audit.csv | Additional abundance filter applied to representatives. |
| convergence.csv, model_diagnostics.json | Per-feature optimizer flags/dispersions and design/normalization/warning audit. |
| size_factors.csv, library_sizes.csv | Fitted size factors and sums over retained representatives (not full sequencing-library totals). |
| pca.csv, sample_correlations.csv | Descriptive PCA coordinates and Pearson sample correlations on log expression. PCA explained variance is in manifest. |
| ranking.csv | Unique tested symbols and signed Wald statistic, descending; stable lexical order for ties without random jitter. |
| gsea.csv | `ES`: enrichment score; `NES`: normalized ES; `pvalue`: nominal permutation p; `padj`: GSEA permutation FDR; `fwer`: permutation family-wise p; `genes`: leading-edge genes. Positive NES points toward tumor-upregulated statistics. |
| ora.csv | `direction`: up/down query; `overlap`: query/set intersection; `set_size`: tested set members; `query_size`, `background_size`: actual universe counts; `fold_enrichment`: observed / expected query fraction; hypergeometric p and within-library/direction BH `padj`; overlapping `genes`. |
| mapping_*.csv | Tested symbols and exact membership in each gene-set library; absence can mean lack of set membership or historical-symbol mismatch. |
| enrichment_status.json, enrichment_curve.json | Library source/status/coverage and the plotted running-score/hit arrays for the first library's lowest-FDR term. |

An empty pathway CSV has headers and an explicit status record. Missing gene p-values are not replaced with ones or zeros. GSEA zero estimates are finite-permutation outcomes, not exact probabilities. Output CSV line endings are LF to keep manifest hashes portable through Git.
