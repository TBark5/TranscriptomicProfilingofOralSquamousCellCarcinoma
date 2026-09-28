"""Synthetic fixtures only; independent cohort routing and replication rules."""
import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from oscc.datasets import get_dataset, bundle_directory
from oscc.geo_validation import read_geo_metadata, read_geo_counts
from oscc.external_validation import compare_genes, compare_pathways, preservation_hashes
from oscc.dashboard import load_results
from oscc.acquisition import sha256
from test_dashboard import write_test_bundle


def soft_fixture(path, conflict=False):
    text = "^SERIES = GSE184616\n"
    for i, patient in enumerate([*range(1, 15), 16]):
        for j, (suffix, condition) in enumerate([("N", "Adjacent Normal"), ("P", "Primary Tumour")]):
            text += f"^SAMPLE = GSM{5593747 + 2*i + j}\n!Sample_title = OSCC_{patient}-{suffix}\n"
            text += "!Sample_source_name_ch1 = Tongue\n"
            fields = {"patient id": f"OSCC_{patient}", "condition": condition,
                      "age": "30", "gender": "Female", "smoking": "No",
                      "patient diagnosis": "HPV-negative Oral Squamous Cell Carcinoma"}
            if conflict and i == 0 and j == 1:
                fields["patient id"] = "OSCC_999"
            text += "".join(f"!Sample_characteristics_ch1 = {k}: {v}\n" for k, v in fields.items())
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_sample_pairing_from_metadata_and_conflicts(tmp_path):
    path = tmp_path / "synthetic.soft.gz"
    soft_fixture(path)
    meta = read_geo_metadata(path)
    assert len(meta) == 30 and meta.patient_id.nunique() == 15
    assert "OSCC_15" not in set(meta.patient_id)
    soft_fixture(path, conflict=True)
    with pytest.raises(ValueError, match="Conflicting"):
        read_geo_metadata(path)


def test_count_alignment_rejects_fractional_duplicate_and_missing_samples(tmp_path):
    soft = tmp_path / "synthetic.soft.gz"
    soft_fixture(soft)
    meta = read_geo_metadata(soft)
    frame = pd.DataFrame(np.arange(60).reshape(2, 30) + 1, index=["SYN_A", "SYN_B"], columns=meta.source_column)
    frame.index.name = "Gene"
    path = tmp_path / "synthetic_counts.txt.gz"
    frame.iloc[:, ::-1].to_csv(path, sep="\t")
    counts = read_geo_counts(path, meta)
    assert list(counts.columns) == list(meta.index)
    assert counts.iloc[0, 0] == 1
    frame.astype(float).add(.5).to_csv(path, sep="\t")
    with pytest.raises(ValueError, match="integer"):
        read_geo_counts(path, meta)
    pd.concat([frame, frame]).to_csv(path, sep="\t")
    with pytest.raises(ValueError, match="duplicate"):
        read_geo_counts(path, meta)
    frame.iloc[:, :-1].to_csv(path, sep="\t")
    with pytest.raises(ValueError, match="headers"):
        read_geo_counts(path, meta)


def synthetic_de(names, lfc, q):
    return pd.DataFrame({"log2FoldChange": lfc, "lfcSE": .2, "pvalue": q, "padj": q,
                         "test_status": "tested"}, index=names)


def test_replication_requires_frozen_candidate_testability_direction_and_q():
    names = list("ABCDEFG")
    discovery = synthetic_de(names, [2, -2, 2, 2, 2, 2, .5], [.01]*7)
    validation = synthetic_de(list("ABCFG"), [1.5, 2, 2, 2, 3], [.01, .01, .1, np.nan, .01])
    table = compare_genes(discovery, validation, list("ABCEFG"), list("ABCFG"))
    assert table.loc["A", "replicated"]
    assert table.loc["B", "validation_outcome"] == "significant_opposite_direction"
    assert table.loc["C", "validation_outcome"] == "not_significant"
    assert table.loc["D", "mapping_status"] == "not_in_validation_source"
    assert table.loc["E", "mapping_status"] == "filtered_low_expression"
    assert not table.loc["F", "testable"]
    assert not table.loc["G", "discovery_candidate"]
    assert table.replicated.sum() == 1


def test_pathway_replication_and_untestable_terms():
    discovery = pd.DataFrame({"library": ["H"]*3, "Term": list("ABC"), "NES": [2, -2, 2], "pvalue": [.01]*3, "padj": [.01]*3})
    validation = pd.DataFrame({"library": ["H"]*3, "Term": list("ABD"), "NES": [1.5, 2, 2], "pvalue": [.01]*3, "padj": [.01]*3})
    compared = compare_pathways(discovery, validation).set_index("Term")
    assert compared.loc["A", "replicated"]
    assert not compared.loc["B", "replicated"]
    assert not compared.loc["C", "testable"]
    assert not compared.loc["D", "discovery_significant"]


def test_dataset_write_paths_do_not_target_original_snapshot(tmp_path):
    discovery = get_dataset("GSE20116")
    validation = get_dataset("GSE184616")
    assert discovery.output(tmp_path) != tmp_path / "results"
    assert discovery.output(tmp_path) != validation.output(tmp_path)
    assert discovery.raw(tmp_path) != validation.raw(tmp_path)
    assert bundle_directory(tmp_path, "GSE20116", legacy=True) == tmp_path / "results"
    assert bundle_directory(tmp_path, "GSE184616") == validation.output(tmp_path)
    with pytest.raises(ValueError):
        get_dataset("GSE184616.1")


def test_validation_loader_never_falls_back_to_discovery(tmp_path, fixture_data):
    write_test_bundle(tmp_path, fixture_data)
    with pytest.raises(FileNotFoundError):
        load_results(tmp_path, "GSE184616")


def test_legacy_loader_checks_all_manifest_files(tmp_path, fixture_data):
    write_test_bundle(tmp_path, fixture_data)
    path = tmp_path / "results/manifest.json"
    manifest = json.loads(path.read_text())
    extra = tmp_path / "results/tables/extra.csv"
    extra.write_text("value\n1\n")
    manifest["output_sha256"]["results/tables/extra.csv"] = sha256(extra)
    path.write_text(json.dumps(manifest))
    load_results(tmp_path)
    extra.write_text("value\n2\n")
    with pytest.raises(ValueError, match="checksum"):
        load_results(tmp_path)


def test_preservation_ignores_new_outputs_but_detects_original_changes(tmp_path):
    original = tmp_path / "results/tables"
    original.mkdir(parents=True)
    (original / "original.csv").write_text("original")
    (tmp_path / "REPORT.md").write_text("original report")
    before = preservation_hashes(tmp_path)
    new = get_dataset("GSE184616").output(tmp_path)
    new.mkdir()
    (new / "new.csv").write_text("new")
    assert preservation_hashes(tmp_path) == before
    (tmp_path / "REPORT.md").write_text("changed")
    assert preservation_hashes(tmp_path) != before


def test_discovery_pipeline_failure_cannot_overwrite_originals(tmp_path, monkeypatch):
    from oscc import pipeline
    original = tmp_path / "results"
    original.mkdir()
    (original / "manifest.json").write_text("original manifest")
    (tmp_path / "REPORT.md").write_text("original report")
    before = preservation_hashes(tmp_path)
    def unavailable(*args, **kwargs):
        raise OSError("synthetic offline source missing")
    monkeypatch.setattr(pipeline, "acquire", unavailable)
    with pytest.raises(OSError, match="synthetic"):
        pipeline.run(tmp_path, offline=True)
    assert preservation_hashes(tmp_path) == before
    assert json.loads((original / "GSE20116/manifest.json").read_text())["status"] == "failed"


def test_failed_external_run_replaces_stale_success_report(tmp_path):
    from oscc.external_validation import run_validation
    (tmp_path / "VALIDATION_REPORT.md").write_text("stale success")
    with pytest.raises(FileNotFoundError):
        run_validation(tmp_path, offline=True)
    assert "**failed**" in (tmp_path / "VALIDATION_REPORT.md").read_text()
    manifest = json.loads((tmp_path / "results/GSE184616/manifest.json").read_text())
    assert manifest["status"] == "failed"
