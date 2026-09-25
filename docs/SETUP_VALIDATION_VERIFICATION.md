# Setup and validation verification

Date: 2026-09-24. Project root: `C:\Users\dingi\OneDrive\Documents\PycharmProjects\TranscriptomicProfilingofOralSquamousCellCarcinoma`. All commands used Python 3.11.13 through `.venv/Scripts/python.exe`; dependencies came from `requirements-lock.txt`, followed by an editable `--no-deps` project install. The prior broken environment remains as `.venv-before-repair`. PyCharm itself was not controlled automatically; set the interpreter to `C:\Users\dingi\OneDrive\Documents\PycharmProjects\TranscriptomicProfilingofOralSquamousCellCarcinoma\.venv\Scripts\python.exe` (see [setup](SETUP_AND_VALIDATION.md)).

## Checks

| Check | Result |
| --- | --- |
| Locked dependencies and editable project install | Passed |
| `pip check` | Passed: no broken requirements |
| Original test suite before changes | Passed: 25 tests |
| Extended test suite | Passed: 35 tests |
| Final test suite | Passed: 35 tests; third-party Matplotlib/seaborn deprecation warnings |
| Original GSE20116 integrity and source hashes | Passed |
| GSE184616 count and metadata audit | Passed: 59,050 unique-symbol count rows, 30 aligned libraries, 15 pairs |
| Paired GSE184616 PyDESeq2 and Hallmark analysis | Completed successfully |
| Validation and discovery output integrity | Passed |
| Dashboard smoke test | Passed for both cohorts, including 30 / 15 sample/pair metric |
| New separate-directory discovery rerun | Passed; all 18 result CSV hashes match the original |
| Original files | All 47 checksummed REPORT, results, and processed-data files remain byte-for-byte unchanged |

The first pytest attempt used Windows' shared `%LOCALAPPDATA%\Temp\pytest-of-dingi`; permissions prevented seven fixtures from being created (18 passed). A fresh `.tools/setup-repair` temporary directory resolved this and all 25 original tests passed. The first extended dashboard smoke test indexed the gene-set library as the gene dropdown; setting `MMP1` failed. A stable `gene-selector` key fixed the selector and the discovery and validation dashboard checks passed. Both were corrected test/setup issues; neither blocked statistical validation.

## Executed GSE184616 analysis

The GEO supplementary counts and SOFT metadata were inspected and SHA256 locked. All sample accessions, integer nonnegative finite counts, title/condition correspondence, matched patients, and exact count-column to GEO-sample alignment were checked before modeling. No rounding, imputation, duplicate summation, or pooled counts were used. The series describes HPV-negative patients; sample GSM5593760 does not repeat HPV-negative in its tumor diagnosis characteristic, which remains explicitly “not stated in sample”.

- 59,050 source genes; 17,921 passed the 10-read / 15-sample filter; 17,850 finite convergence-qualified, Cook-filtered tests.
- Model design rank 16; residual degrees of freedom 14. 71 convergence failures were excluded from BH significance calls.
- Frozen discovery candidates 1,339; testable in validation 1,187; same direction 964; replicated under same direction and validation genome-wide BH q < 0.05 706; also with |validation LFC| > 1 563.
- Discovery-significant Hallmark pathways 13; testable in both cohorts 13; replicated by matching direction and FDR < 0.05 in both 10.
- Descriptive LFC Spearman correlation across 8,982 mutually tested genes: 0.5243.

PyDESeq2 recorded 71 feature-level optimizer failures and a pandas dtype deprecation warning; those failed fits were not treated as significant. GSEApy reported ties for 0.13% of ranked genes; alphabetical tie-breaking was used without jitter. Warnings and status are retained in the validation manifest.

See `VALIDATION_REPORT.md` and `results/GSE184616/` for measured tables, figures, diagnostics, and source hashes. This reports independent expression and pathway replication in the external cohort, not causal or clinical validation.
