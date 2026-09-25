import pandas as pd
import requests
from oscc.enrichment import run_enrichment


def test_unavailable_library_writes_honest_empty_outputs(tmp_path, monkeypatch):
    def unavailable(*args, **kwargs):
        raise requests.ConnectionError("synthetic unavailable-service fixture")
    monkeypatch.setattr("oscc.enrichment.download", unavailable)
    result = pd.DataFrame({"stat": [2., -2.], "pvalue": [.1, .1], "padj": [.2, .2],
                           "log2FoldChange": [1.5, -1.5], "model_converged": [True, True]}, index=["SYN_A", "SYN_B"])
    (tmp_path / "enrichment_curve.json").write_text("stale test fixture")
    ora, gsea, status, curve = run_enrichment(result, tmp_path, tmp_path, ["MSigDB_Hallmark_2020"])
    assert ora.empty and gsea.empty and curve is None
    assert status[0]["status"] == "unavailable"
    assert "padj" in pd.read_csv(tmp_path / "gsea.csv").columns
    assert not (tmp_path / "enrichment_curve.json").exists()
