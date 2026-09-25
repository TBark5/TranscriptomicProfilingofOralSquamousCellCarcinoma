# Data provenance and acquisition decision

## Source inspection

Inspected on 2026-09-23: the actual downloaded GEO family SOFT and original publication workbook. GEO reports zero sample expression-table rows, the AB whole-transcriptome pipeline, hg18 alignments in MAX format, six per-sample `.sorted.max.txt.gz` supplements, and the series MAX specification/archive. The roughly 3.9 GB archive was not downloaded.

The original publication provides [Table S1, supplement s009](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0009317.s009&type=supplementary). This workbook was verified to contain **15,668 data rows and 37 columns**, with three header rows, in sheet `counts.annotated.byTranscript`. It provides raw transcript counts, annotations, normalized expression and the original analysis side by side. The pipeline uses **A:C annotations and D:I raw count sums only**; it does not infer counts by reversing normalization or copy the original p-values.

The [edgeR user guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) independently documents using these raw counts in its paired oral-carcinoma example. The old `https://bioinf.wehi.edu.au/edgeR/TableS1.txt` URL returned HTTP 404, so the project uses the original publication supplement directly. The first observed count row is NM_182502/TMPRSS11B, 8N=2592, 8T=3, 33N=7805, 33T=321, 51N=3372, 51T=9, agreeing with the worked example. The related microarray accession GSE19089 is **not used**.

## Accession and sample mapping

Series [GSE20116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20116); platform GPL9442; SRA SRP002009; BioProject PRJNA124251. Samples are parsed from GEO titles and checked against patient characteristics, not inferred from column position alone.

| Sample accession | Patient ID | Condition | Workbook raw-count column | GEO title |
|---|---|---|---|---|
| GSM515513 | 8 | normal | D / 8N | normal_tissue_patient_8 |
| GSM515514 | 8 | tumor | E / 8T | tumor_tissue_patient_8 |
| GSM515515 | 33 | normal | F / 33N | normal_tissue_patient_33 |
| GSM515516 | 33 | tumor | G / 33T | tumor_tissue_patient_33 |
| GSM515517 | 51 | normal | H / 51N | normal_tissue_patient_51 |
| GSM515518 | 51 | tumor | I / 51T | tumor_tissue_patient_51 |

GEO supplies tissue, patient ID, instrument, extraction/processing descriptions and library strategy. Age, sex, stage and HPV status are not supplied in the SOFT records and are recorded as missing. The original paper describes negative surgical-margin normal tissue. No per-patient anatomic subsites, purity estimates or clinical stages are invented.

## Original processing and this project's changes

The source experiment used rRNA-depleted total RNA and AB SOLiD 3.0 sequencing. Its pipeline quantified uniquely aligned reads over RefSeq exons on hg18. The paper describes filtering to at least 50 reads in at least one tissue from every patient. **Actual workbook audit:** 47 transcript rows have a within-patient maximum of 49 in at least one pair; every row meets 49 in every pair. This boundary discrepancy is consistent with, but does not prove, use of the paper's pseudocount during filtering. We retain the supplied raw measurements unchanged, rather than imposing a correction to reproduce the prose description. The published normalized values include a pseudocount; the raw `sum` columns used here include genuine zeros.

This project retains historical source symbols. It selects maximum-exon-count transcripts, resolving ties by lexical RefSeq ID independently of expression/statistics. There are 10,541 selected symbols, 5,126 redundant transcript rows, and one rejected ambiguous symbol (`NOP5/NOP58`). No isoform counts are added. Thus, the result unit is one transcript representative per gene symbol. This selection is inspired by the documented edgeR case but does not reproduce its old Bioconductor reannotation: modern aliases are not guessed, and counts are not presented as modern gene-union quantifications.

The additional default rule requires >=10 raw reads in >=3 samples; all 10,541 representatives in this source snapshot pass. `transcript_mapping.csv` preserves the selected/rejected/redundant rows and original RefSeq IDs, and `filter_audit.csv` records the abundance decision.

## Download and integrity

```bash
oscc download
oscc run
oscc run --offline
```

The first command caches workbook/SOFT; `run` also caches requested GMTs. Each source download is limited to 30 MB, written atomically and verified. The checked-in source hashes below are enforced in code, including when using cached files. Upstream changes fail visibly rather than silently changing the analysis.

| File | URL | SHA256 |
|---|---|---|
| TableS1.xls | https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0009317.s009&type=supplementary | `5aa416fc86ebf716252f7df55e34c8d33dbe035b92889cfdfd879332eb289ffd` |
| GSE20116_family.soft.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20116/soft/GSE20116_family.soft.gz | `7bfb1cdc0299e2529eb2a79ab22d519201aa64d17fb5d1844ff039fabf057af0` |
| MSigDB_Hallmark_2020.gmt | https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020 | `4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81` |

The completed run's manifest includes exact source byte sizes, hashes, verification time and software versions. A verification timestamp is not necessarily the first download time. Raw sources are deliberately ignored by Git. If downloading manually, save the URLs above under `data/raw/` with exactly the listed filenames, then run offline. HTML error pages or changed bytes fail checksum checks.

## Gene sets

Default: Enrichr's human `MSigDB_Hallmark_2020` GMT snapshot, 50 source sets. In this run, 49 sets have 15–500 tested members. Of 10,506 eligible tested source symbols, 2,932 occur in at least one Hallmark set; this is **set-membership coverage**, not evidence that every other identifier failed annotation. Exact historical-symbol matching can miss renamed genes; `mapping_MSigDB_Hallmark_2020.csv` records inclusion without speculative remapping. The RefSeq-to-symbol conversion itself comes from the original workbook.

Optional human libraries `GO_Biological_Process_2023` and `Reactome_2022` use the same endpoint. They are not part of the completed analysis snapshot. Their retrieved hashes are recorded at runtime, but are not pre-locked in code; preserve their cached GMTs if using them for later exact reproduction. Gene-set unavailability is explicit in `enrichment_status.json`. Any results from optional collections must be interpreted with their own testing multiplicity and overlap.

## Attribution

Tuch BB et al. (2010), *Tumor Transcriptome Sequencing Reveals Allelic Expression Imbalances Associated with Copy Number Alterations*, PLOS ONE 5:e9317, https://doi.org/10.1371/journal.pone.0009317. The article and supplement are published under the source's Creative Commons Attribution terms. Project MIT licensing does not replace source-data or gene-set terms. Gene-set definition files are not redistributed in Git.
