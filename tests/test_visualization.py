import numpy as np
import pandas as pd
from oscc.visualization import safe_neglog, top_genes, make_figures


def test_zero_and_missing_q_values():
    y = safe_neglog([0., .01, np.nan])
    assert y.iloc[0] == 300
    assert y.iloc[1] == 2
    assert np.isnan(y.iloc[2])


def test_empty_significance_figures(tmp_path, fixture_data):
    counts, metadata = fixture_data
    result = pd.DataFrame({"log2FoldChange": np.zeros(len(counts)), "padj": np.ones(len(counts)), "pvalue": np.ones(len(counts))}, index=counts.index)
    genes, title = top_genes(result)
    assert "exploratory" in title
    assert len(genes) == len(counts)
    make_figures(counts, np.log2(counts + 1), metadata, result, pd.DataFrame(), None, tmp_path)
    assert (tmp_path / "volcano.svg").exists()
    assert (tmp_path / "pathways.png").exists()
