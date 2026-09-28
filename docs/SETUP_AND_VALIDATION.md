# Repaired setup and independent validation

Current project root:
`C:\Users\dingi\OneDrive\Documents\PycharmProjects\TranscriptomicProfilingofOralSquamousCellCarcinoma`

The environment was recreated with Python 3.11.13, matching the project's `>=3.11,<3.12` requirement. The old environment is retained as `.venv-before-repair`; no original analysis files were removed. Dependencies are installed from the existing `requirements-lock.txt` and the project is installed editable with `--no-deps`.

## PyCharm interpreter
Select **Settings > Project > Python Interpreter > Add Interpreter > Add Local Interpreter > Existing**, then use this exact executable:

`C:\Users\dingi\OneDrive\Documents\PycharmProjects\TranscriptomicProfilingofOralSquamousCellCarcinoma\.venv\Scripts\python.exe`

The IDE's global interpreter registry has not been changed automatically. Set the working directory to the project root above. Avoid using the interpreter from `.venv-before-repair`.

## Separate datasets and commands
Run these commands in PowerShell from the project root:

```powershell
.venv\Scripts\python.exe -m pip check
.venv\Scripts\python.exe -m oscc.cli validate --dataset GSE20116
.venv\Scripts\python.exe -m oscc.cli download --dataset GSE184616
.venv\Scripts\python.exe -m oscc.cli external-validation --offline
.venv\Scripts\python.exe -m oscc.cli validate --dataset GSE184616
.venv\Scripts\python.exe -m pytest -q
.venv\Scripts\python.exe -m streamlit run app/streamlit_app.py
```

`run --dataset GSE184616` is equivalent to `external-validation`. Validation uses fixed thresholds: 10 reads in 15 samples, q < 0.05, descriptive |LFC| > 1, and the locked Hallmark library. Only seed and permutations are customizable; incompatible options fail rather than being ignored. Discovery candidates use the thresholds recorded in the original discovery manifest. Validation primary replication uses whole-validation-family BH q < 0.05 and the discovery effect direction.

| Dataset / artifact | Location |
| --- | --- |
| Original discovery snapshot (preserved) | `results/tables`, `results/figures`, `results/manifest.json`, `REPORT.md` |
| New discovery reruns | `results/GSE20116/tables`, `results/GSE20116/figures`, `results/GSE20116/REPORT.md` |
| Validation cohort | `results/GSE184616/tables`, `results/GSE184616/figures`, `results/GSE184616/manifest.json` |
| Validation report | `VALIDATION_REPORT.md` |
| Validation inputs | `data/raw/GSE184616` |
| New discovery processed matrices | `data/processed/GSE20116` |

The dashboard has a cohort selector and reads matching manifests and labels. For discovery it reads the preserved snapshot until a new complete discovery run exists. An incomplete newer run is rejected, not silently replaced. Independent validation always reads the original discovery snapshot. Cohorts are normalized and fitted separately; count matrices are not concatenated.

Re-running `oscc run --offline` now writes only the new discovery directory and its report; it does not overwrite the original REPORT.md. Do not run concurrent writers for one cohort. Check manifest status before using outputs. A failed validation attempt writes a failure report rather than leaving a stale success report.

Execution evidence and any failures are recorded in `docs/SETUP_VALIDATION_VERIFICATION.md` after checks complete. The historical `docs/VERIFICATION.md` describes the original pre-extension analysis.
