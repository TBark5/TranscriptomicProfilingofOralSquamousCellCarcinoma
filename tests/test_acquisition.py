import numpy as np
import pandas as pd
import pytest
from oscc.acquisition import download, sha256, read_publication_counts, SAMPLE_COLUMNS


def test_cached_download_integrity(tmp_path):
    path = tmp_path / "source.bin"
    path.write_bytes(b"explicitly synthetic download fixture")
    info = download("https://example.invalid/fixture", path, sha256(path), offline=True)
    assert info["bytes"] == len(path.read_bytes())
    with pytest.raises(ValueError, match="Checksum"):
        download("https://example.invalid/fixture", path, "wrong", offline=True)
    with pytest.raises(FileNotFoundError, match="Offline"):
        download("https://example.invalid/fixture", tmp_path / "absent", offline=True)


def test_publication_parser_selects_raw_columns(monkeypatch):
    # Synthetic workbook-shaped data, deliberately different normalized columns.
    table = pd.DataFrame(np.nan, index=range(4), columns=range(15), dtype=object)
    table.iloc[0, 3:9] = SAMPLE_COLUMNS
    table.iloc[2, :3] = ["idRefSeq", "nameOfGene", "numberOfExons"]
    table.iloc[2, 3:9] = ["sum"] * 6
    table.iloc[3, :9] = ["NM_000001", "SYN_GENE", 4, 0, 2, 3, 4, 5, 6]
    table.iloc[3, 9:15] = [.001] * 6
    monkeypatch.setattr(pd, "read_excel", lambda *args, **kwargs: table)
    counts = read_publication_counts("synthetic-only.xls")
    assert counts[SAMPLE_COLUMNS].iloc[0].tolist() == [0, 2, 3, 4, 5, 6]
    table.iloc[2, 3] = "normalized_sum"
    with pytest.raises(ValueError, match="raw-count headers"):
        read_publication_counts("synthetic-only.xls")
