"""Validated dashboard data access, independently testable without Streamlit."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .acquisition import sha256
from .validation import validate_counts
from .differential import classify
from .datasets import bundle_directory


def load_results(root, dataset="GSE20116", legacy=False):
    root = Path(root).resolve()
    output = bundle_directory(root, dataset, legacy)
    manifest_path = output / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError("No analysis manifest. Run `oscc run` from the repository first.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("data_kind") != "real":
        raise ValueError("A completed real-data analysis is required. Inspect the manifest and rerun `oscc run`.")
    if manifest.get("dataset") != dataset:
        raise ValueError("Manifest dataset does not match requested cohort")
    filenames = ["differential_expression", "metadata", "counts", "normalized_counts", "log_expression", "gsea", "ora"]
    bundle = {"manifest": manifest, "output_directory": output}
    for name in filenames:
        relative = (output / "tables" / f"{name}.csv").relative_to(root).as_posix()
        path = root / relative
        if not path.exists():
            raise FileNotFoundError(f"Missing {relative}; rerun the pipeline")
        expected = manifest.get("output_sha256", {}).get(relative)
        if not expected or sha256(path) != expected:
            raise ValueError(f"Output checksum mismatch for {relative}; rerun the complete pipeline")
        bundle[name] = pd.read_csv(path, index_col=None if name in ["gsea", "ora"] else 0,
                                   dtype={"patient_id": str} if name == "metadata" else None)
    for relative, expected in manifest.get("output_sha256", {}).items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Output checksum mismatch for {relative}")
    counts, metadata = validate_counts(bundle["counts"], bundle["metadata"])
    result = bundle["differential_expression"]
    if result.index.has_duplicates or not result.index.equals(counts.index):
        raise ValueError("DE result identifiers do not match the expression matrix")
    for name in ["normalized_counts", "log_expression"]:
        expression = bundle[name]
        if not expression.index.equals(counts.index) or not expression.columns.equals(counts.columns):
            raise ValueError("Expression output ordering differs from counts")
        if not np.isfinite(expression.to_numpy()).all():
            raise ValueError("Nonfinite expression output")
    if not np.allclose(np.log2(bundle["normalized_counts"] + 1), bundle["log_expression"]):
        raise ValueError("Log-expression transformation is inconsistent")
    bundle["metadata"] = metadata
    return bundle


def filtered_results(result, alpha, lfc, search=""):
    out = result.copy()
    out["direction"] = classify(out, alpha, lfc)
    out = out.loc[out.direction != "not_significant"]
    if search.strip():
        out = out.loc[out.index.str.contains(search.strip(), case=False, regex=False)]
    return out.sort_values("padj")


def gene_expression(bundle, gene):
    if gene not in bundle["log_expression"].index:
        raise KeyError(f"Gene {gene} is not in the analyzed source-symbol universe")
    out = bundle["metadata"][["patient_id", "condition"]].copy()
    out["log_expression"] = bundle["log_expression"].loc[gene]
    out["normalized_count"] = bundle["normalized_counts"].loc[gene]
    return out.reset_index()
