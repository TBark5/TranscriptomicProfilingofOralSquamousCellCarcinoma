"""One-command real-data pipeline. A completion manifest gates dashboard loading."""
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from .acquisition import acquire, read_metadata, read_publication_counts, SAMPLE_COLUMNS, sha256
from .preprocessing import select_representatives, filter_counts, log_expression
from .validation import validate_counts
from .differential import fit_paired, classify
from .enrichment import run_enrichment
from .visualization import make_figures
from .report import write_report
from .datasets import get_dataset


def run(root, alpha=.05, lfc=1., min_count=10, min_samples=None, seed=42, permutations=1000,
        libraries=None, offline=False, dataset="GSE20116"):
    cohort = get_dataset(dataset)
    if cohort.role == "validation":
        if alpha != .05 or lfc != 1. or min_count != 10 or min_samples not in (None, 15) or libraries not in (None, ["MSigDB_Hallmark_2020"]):
            raise ValueError("External validation uses fixed thresholds, 10 reads in 15 samples, and Hallmark; only seed/permutations may vary")
        from .external_validation import run_validation
        return run_validation(root, offline=offline, seed=seed, permutations=permutations)
    min_samples = cohort.min_samples if min_samples is None else min_samples
    if not 0 < alpha < 1 or lfc < 0 or permutations < 100:
        raise ValueError("Require 0 < alpha < 1, lfc >= 0, and at least 100 GSEA permutations")
    root = Path(root).resolve()
    output = cohort.output(root)
    tables, figures = output / "tables", output / "figures"
    processed = root / "data/processed" / dataset
    for directory in (tables, figures, processed):
        directory.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "manifest.json"
    manifest = {"status": "running", "dataset": dataset, "role": cohort.role, "data_kind": "real",
                "started_utc": datetime.now(timezone.utc).isoformat(),
                "config": {"alpha": alpha, "lfc": lfc, "min_count": min_count, "min_samples": min_samples,
                           "seed": seed, "permutations": permutations, "libraries": libraries or ["MSigDB_Hallmark_2020"]}}
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    try:
        np.random.seed(seed)
        print("Acquiring and validating GSE20116 publication counts...", flush=True)
        manifest["sources"] = acquire(root, offline, dataset)
        metadata = read_metadata(root / "data/raw/GSE20116_family.soft.gz")
        source = read_publication_counts(root / "data/raw/TableS1.xls")
        selected, mapping = select_representatives(source)
        manifest.update(source_transcripts=len(source), representatives=len(selected), samples=len(metadata), pairs=metadata.patient_id.nunique())
        counts = selected[SAMPLE_COLUMNS].rename(columns=dict(zip(metadata.source_column, metadata.index)))
        counts, metadata = validate_counts(counts, metadata)
        counts, filter_audit = filter_counts(counts, min_count, min_samples)
        counts.index.name = "gene"
        mapping.to_csv(tables / "transcript_mapping.csv", index=False, lineterminator='\n')
        filter_audit.to_csv(tables / "filter_audit.csv", index_label="gene", lineterminator='\n')
        counts.to_csv(processed / "counts.csv", lineterminator='\n')
        metadata.to_csv(processed / "metadata.csv", lineterminator='\n')
        print(f"Fitting {len(counts):,} representatives across {len(metadata)} matched samples...", flush=True)
        result, normalized, sizes, convergence, diagnostics = fit_paired(counts, metadata)
        result = result.join(selected[["refseq", "exons"]])
        result.index.name = "gene"
        result["direction"] = classify(result, alpha, lfc)
        result.to_csv(tables / "differential_expression.csv", lineterminator='\n')
        for direction in ("up", "down"):
            result.loc[result.direction == direction].sort_values("padj").to_csv(tables / f"{direction}regulated.csv", lineterminator='\n')
        log_expr = log_expression(normalized)
        for name, frame in [("counts", counts), ("normalized_counts", normalized), ("log_expression", log_expr),
                             ("metadata", metadata), ("convergence", convergence)]:
            frame.to_csv(tables / f"{name}.csv", lineterminator='\n')
        sizes.to_csv(tables / "size_factors.csv", lineterminator='\n')
        pd.DataFrame({"retained_count_sum": counts.sum(), "size_factor": sizes}).to_csv(tables / "library_sizes.csv", index_label="sample_accession", lineterminator='\n')
        manifest["diagnostics"] = diagnostics
        (tables / "model_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
        print("Running local gene-set enrichment...", flush=True)
        ora, gsea, enrichment_status, curve = run_enrichment(result, root, tables, manifest["config"]["libraries"], alpha, lfc, seed, permutations, offline)
        manifest["enrichment"] = enrichment_status
        print("Rendering PNG/SVG figures and research report...", flush=True)
        pca, explained, correlations = make_figures(counts, log_expr, metadata, result, gsea, curve, figures, alpha, lfc, dataset=dataset)
        pca.to_csv(tables / "pca.csv", index_label="sample_accession", lineterminator='\n')
        correlations.to_csv(tables / "sample_correlations.csv", lineterminator='\n')
        off_diagonal = correlations.to_numpy()[~np.eye(len(metadata), dtype=bool)]
        qc = {"library_min": int(counts.sum().min()), "library_max": int(counts.sum().max()),
              "correlation_min": float(off_diagonal.min()), "correlation_max": float(off_diagonal.max()), "pca_variance": explained}
        manifest["qc"] = qc
        manifest["summary"] = {"genes": len(result), "tested": int(result.pvalue.notna().sum()),
                               "up": int((result.direction == "up").sum()), "down": int((result.direction == "down").sum())}
        manifest["versions"] = {x: importlib.metadata.version(x) for x in ["numpy", "pandas", "scipy", "statsmodels", "pydeseq2", "gseapy", "matplotlib", "streamlit"]}
        manifest["python"] = platform.python_version()
        write_report(root, manifest, result, gsea, qc, ora, report_path=output / "REPORT.md")
        manifest["output_sha256"] = {str(p.relative_to(root)).replace("\\", "/"): sha256(p) for p in tables.glob("*.csv")}
        manifest["status"] = "complete"
        manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
        print(json.dumps(manifest["summary"], indent=2), flush=True)
    except Exception as exc:
        manifest["status"], manifest["error"] = "failed", str(exc)
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
