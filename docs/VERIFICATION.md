# Executed verification

Environment: Windows, Python 3.11.13, local isolated environment. Verification date: 2026-09-23. The default system `python` alias was not runnable and the installed Python 3.14 was unsuitable for the pinned stack; a local Python 3.11 runtime was installed. `.tools/` and `.venv/` are excluded from Git.

| Check actually executed | Outcome |
|---|---|
| Install pinned requirements into the local environment | Successful (initial installation with uv; `python -m pip install -r requirements.txt` subsequently succeeded with requirements satisfied). |
| Editable project installation | Successful; `oscc` entry point runs. |
| `python -m pip check` | No broken requirements found. |
| `oscc download --root .tools/download-check` | Successfully downloaded a fresh publication workbook and GEO SOFT through the project's own downloader, with both expected SHA256 hashes. |
| `oscc run --offline` | Completed the entire real GSE20116 pipeline, including Hallmark GSEA/ORA, tables, report and eight PNG/SVG figure types. |
| `oscc run` | Completed using the verified source cache. |
| `oscc validate` | Completed output hashes, metadata alignment and expression consistency validated. |
| `python -m pytest -q --junitxml=results/test-results.xml` | **25 passed**, 15 third-party plotting deprecation warnings; 16.28 seconds in the recorded final test run. |
| `python scripts/execute_notebook.py` | **10 code cells executed without errors**; outputs saved in the notebook. |
| `python scripts/check_dashboard.py` | Initial render, both gene-threshold sliders, gene search, all-gene table, candidate selection, ORA mode and pathway search passed without app exceptions. |
| `streamlit run app/streamlit_app.py --server.headless true --server.address 127.0.0.1 --server.port 8501` | Started successfully; localhost `/_stcore/health` returned HTTP 200 and `ok`. |
| `python -m compileall -q src app scripts` | Completed successfully. |
| Figure inspection | PNG volcano, clustered heatmap, PCA and pathway plots visually inspected; overlapping volcano labels were revised. All eight PNG and SVG pairs were generated. |
| Deterministic rerun | All **18 generated CSV SHA256 hashes** matched exactly between the final two complete runs. Manifest timestamps intentionally differ. CSVs contain LF line endings. |
| Original prefilter audit | All 15,668 source rows meet a 49-read within-patient maximum in every pair; 47 fall below the paper's stated 50-read boundary. Documented without modifying counts. |

The tests use explicitly synthetic fixtures or generated NB counts in temporary test directories. None are loaded by the production CLI. One initial test failed because pandas attached an index name during metadata alignment; its assertion was corrected to compare data/order without incidental index names. A dashboard smoke script initially selected sliders by the wrong widget-list order; named widget keys now make the intended controls explicit. These issues were corrected and the checks rerun successfully.

## Scientific execution limits

The real fit retained 10,541 gene representatives, with 10,506 finite convergence-qualified tests and 35 excluded optimizer failures. The PyDESeq2 warning about dispersion estimation with fewer than three residual degrees of freedom is retained in `results/tables/model_diagnostics.json`. A pandas compatibility warning emitted by PyDESeq2 is also retained. The analysis does not suppress these scientific/model warnings or count failed fits as significant.

Hallmark is the only completed real-data gene-set collection. Optional GO and Reactome branches are implemented but were not run against real data. Gene-set-service failure and honest empty output handling are covered by a test. The 3.9 GB MAX archive was not downloaded, and raw-read/alignment QC was not redone. No external-cohort or laboratory validation was performed.

Network downloads required the environment's explicit network permission; they succeeded when allowed. No credentials are stored. The GitHub Actions workflow is supplied but has not run remotely. Linux/macOS installation, browser-engine rendering of Plotly, and public hosting have not been verified. Streamlit checks cover its Python execution/widget interactions plus HTTP server startup, not browser screenshot testing.

This project is initialized as a local Git repository with files ready for review. No GitHub remote, commit, public deployment or push was created.
