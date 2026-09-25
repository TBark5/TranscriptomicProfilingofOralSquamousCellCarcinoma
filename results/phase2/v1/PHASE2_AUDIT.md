# Phase 2 pre-change audit

The original discovery and validation outputs were inventoried and SHA-256 hashed before Phase 2 edits. The baseline is `results/phase2/v1/frozen_sha256.json`. The repository had no commits before the `oscc-phase2-validation` branch was created, so this manifest is the file-level baseline.

## Confirmed strengths

- Discovery GSE20116 has three matched tumor/normal pairs; independent GSE184616 has 15 matched pairs. The public GSM accession sets do not overlap. Validation sample names and pairing were checked against the deposited GEO SOFT and count-file column names; the saved source hashes match the acquired files.
- Both fitted analyses use a patient-blocked PyDESeq2 model, explicit tumor/normal contrast, convergence flags, Cook filtering, median-of-ratios normalization, and BH correction on finite convergence-qualified tests. Same locked Hallmark GMT file and GSEA seed were used.
- Recalculation from saved tables reproduces 1,339 discovery candidates, 1,187 testable in validation, 964 same-direction, 706 primary replications, 563 stricter replications, and 10/13 originally significant Hallmark pathways replicated. Recomputed BH values agree numerically with saved values.

## Methodological concerns and missing evidence

- Three discovery pairs give only two residual degrees of freedom. Wald p-values and extreme fold changes are fragile to dispersion assumptions and individual patients. They are not proof of stable population-level effects.
- The validation table's five mutually exclusive candidate outcomes are: 152 untestable, 706 same-direction significant, 258 same-direction nonsignificant, **61 opposite-direction significant**, and 162 opposite-direction nonsignificant. The 61 contradictory significant effects need prominent reporting.
- Historical discovery RefSeq representatives were assigned gene symbols, while validation counts use modern annotation symbols. Exact symbol intersection avoids speculative remapping but loses renamed genes and may mask many-to-one or ambiguous assignments. The original mapping lacks a modern authoritative alias audit.
- Hallmark files match, but discovery and validation GSEA backgrounds differ. Original pathway replication is thus a comparison over cohort-specific measured universes. GSEA's saved zero empirical p/FDR values have finite permutation resolution and must not be described as exact zeros.
- Bulk tumor versus surgical-margin normal differences may reflect tissue composition, tumor purity, and site. The public sample metadata do not establish individual identity across cohorts, cell of origin, or experimentally validated biomarkers. One validation sample lacks a sample-level HPV-negative label although its series is described as HPV-negative.
- The original validation report highlights only six genes and does not show all ten now prespecified. No original paired expression plots, common-universe enrichment sensitivity, tissue-marker exploration, or subsite sensitivity cover the full requested scope.

## Planned corrections

Keep original results frozen. Derive Phase 2 tables and figures solely in a versioned directory; verify arithmetic and model statistics against saved outputs; audit HGNC exact/previous/alias mappings without changing primary calls; show every prespecified gene and all outcome categories; compare effects and pathways under a common measurable universe; run exploratory threshold, tissue-marker, and tongue-subset checks; and report uncertainty and infeasible analyses explicitly.
