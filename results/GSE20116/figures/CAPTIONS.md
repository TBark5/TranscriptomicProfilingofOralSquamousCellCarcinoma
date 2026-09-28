# Figure captions

quality_control: Count sums across retained gene representatives are not full sequencing-library sizes. Boxplots show log2(size-factor-normalized counts + 1); outliers hidden for readability.

sample_correlation: Pearson correlations across retained genes using log normalized expression, without patient or condition regression.

pca: Centered, unscaled PCA on up to 2,000 most variable log normalized genes; lines connect matched patients. Selection does not use DE p-values. PCA is descriptive, not a significance test.

volcano: Paired PyDESeq2 Wald tests of tumor versus normal; thresholds BH-adjusted p < 0.05 and |log2 fold change| > 1.0. Missing adjusted p-values are omitted; zeros are clipped to 1e-300 for plotting only. Labels select up to four actual threshold-passing genes per direction. Fold changes are unshrunk.

heatmap: Up to 30 top genes, ranked by adjusted p-value; row z-scores of log normalized expression, average-linkage Euclidean clustering on rows and columns. Color strips identify condition and patient. This DE-selected view is not independent evidence of separation; any non-significant fallback is identified in the title.

paired_expression: Individual patient trajectories for up to four data-selected candidates. Points are log normalized expression, while displayed LFC and q are from the paired count model.

pathways: Up to 15 lowest-FDR preranked GSEA terms, including non-significant terms if necessary. Positive NES indicates enrichment toward positive tumor/normal Wald statistics. FDR is GSEA's permutation-derived estimate, not BH. Gene-set enrichment does not establish activation or causality.

enrichment_curve: Running weighted enrichment score and gene-set hit positions for the first library's lowest-FDR term; nominal zero permutation p-values/FDR are resolution-limited, not certainty.
