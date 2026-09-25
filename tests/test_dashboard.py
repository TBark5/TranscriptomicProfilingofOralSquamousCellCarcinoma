import json
import numpy as np
import pandas as pd
import pytest
from oscc.acquisition import sha256
from oscc.dashboard import load_results, filtered_results, gene_expression


def write_test_bundle(tmp_path, fixture_data):
    """Synthetic schema fixture for loader tests ONLY; never used by the app."""
    counts, metadata = fixture_data
    tables = tmp_path / "results/tables"
    tables.mkdir(parents=True)
    result = pd.DataFrame({"log2FoldChange": [2.] * len(counts), "padj": [.01] * len(counts)}, index=counts.index)
    frames = {"counts": counts, "normalized_counts": counts.astype(float), "log_expression": np.log2(counts + 1),
              "metadata": metadata, "differential_expression": result, "gsea": pd.DataFrame(columns=["Term", "padj"]), "ora": pd.DataFrame(columns=["Term", "padj"])}
    for name, frame in frames.items():
        frame.to_csv(tables / f"{name}.csv", index=name not in ["gsea", "ora"])
    manifest = {"status": "complete", "data_kind": "real", "dataset": "GSE20116", "test_fixture_only": True,
                "output_sha256": {f"results/tables/{p.name}": sha256(p) for p in tables.glob("*.csv")}}
    (tmp_path / "results/manifest.json").write_text(json.dumps(manifest))
    return tables


def test_missing_outputs(tmp_path):
    with pytest.raises(FileNotFoundError, match="manifest"):
        load_results(tmp_path)


def test_loader_and_tampering(tmp_path, fixture_data):
    tables = write_test_bundle(tmp_path, fixture_data)
    bundle = load_results(tmp_path)
    assert len(gene_expression(bundle, "SYN_UP")) == 6
    assert len(filtered_results(bundle["differential_expression"], .05, 1, "UP")) == 1
    with (tables / "counts.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="checksum"):
        load_results(tmp_path)


def test_incomplete_run_is_rejected(tmp_path, fixture_data):
    write_test_bundle(tmp_path, fixture_data)
    path = tmp_path / "results/manifest.json"
    manifest = json.loads(path.read_text())
    manifest["status"] = "failed"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="completed real-data"):
        load_results(tmp_path)
