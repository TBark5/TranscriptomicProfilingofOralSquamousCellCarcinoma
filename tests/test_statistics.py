import numpy as np
import pandas as pd
import pytest
from oscc.differential import bh_adjust, classify, fit_paired
from oscc.enrichment import prepare_ranking, overrepresentation, read_gmt


def test_bh_known_values_and_missing():
    actual = bh_adjust(pd.Series([.01, .04, .03, np.nan]))
    np.testing.assert_allclose(actual.iloc[:3], [.03, .04, .04])
    assert np.isnan(actual.iloc[3])
    assert bh_adjust([0., .5]).iloc[0] == 0
    with pytest.raises(ValueError):
        bh_adjust([1.1])


def test_strict_thresholds_direction():
    result = pd.DataFrame({"padj": [.01, .01, .05, .01, np.nan], "log2FoldChange": [2, -2, 2, 1, 10]})
    assert classify(result).tolist() == ["up", "down", "not_significant", "not_significant", "not_significant"]


def test_ranking_excludes_invalid_and_is_signed():
    result = pd.DataFrame({"pvalue": [.1, .2, np.nan, .4], "stat": [-3, 4, 8, 10], "model_converged": [True, True, True, False]}, index=["A", "B", "C", "D"])
    rank = prepare_ranking(result)
    assert rank.gene.tolist() == ["B", "A"]
    assert rank.stat.tolist() == [4, -3]
    with pytest.raises(ValueError, match="Duplicate"):
        prepare_ranking(pd.concat([result, result.iloc[:1]]))


def test_ora_uses_tested_background_and_zero_overlap():
    result = overrepresentation(["A", "B"], list("ABCDE"), {"hit": list("ABX"), "miss": list("CDX")}, min_size=1)
    assert result.background_size.eq(5).all()
    assert result.set_size.eq(2).all()
    assert result.loc[result.Term == "hit", "pvalue"].iloc[0] == pytest.approx(.1)
    assert result.loc[result.Term == "miss", "pvalue"].iloc[0] == 1.
    assert result.loc[result.Term == "hit", "padj"].iloc[0] == pytest.approx(.2)
    with pytest.raises(ValueError, match="subset"):
        overrepresentation(["X"], ["A"], {"x": ["X"]})


def test_gmt_deduplication(tmp_path):
    path = tmp_path / "synthetic.gmt"
    path.write_text("SYN_SET\tNA\tA\tA\tB\n")
    assert read_gmt(path) == {"SYN_SET": ["A", "B"]}


@pytest.mark.integration
def test_paired_count_model_fold_change_direction(fixture_data):
    _, metadata = fixture_data
    rng = np.random.default_rng(42)
    # Artificial NB measurements for an integration test only.
    mean = rng.uniform(100, 1000, size=(300, 1)) * np.array([.8, .8, 1.2, 1.2, 1., 1.])
    counts = rng.negative_binomial(50, 50 / (50 + mean))
    counts[0] = [100, 850, 130, 1050, 110, 900]
    counts[1] = [900, 100, 1080, 120, 850, 95]
    counts = pd.DataFrame(counts, index=[f"SYN_{i}" for i in range(300)], columns=metadata.index)
    result, normalized, sizes, audit, diagnostics = fit_paired(counts, metadata)
    assert result.loc["SYN_0", "log2FoldChange"] > 2
    assert result.loc["SYN_1", "log2FoldChange"] < -2
    assert diagnostics["design_rank"] == 4
    np.testing.assert_allclose(normalized.to_numpy(), counts.div(sizes, axis=1).to_numpy())
    assert result.loc[~audit.model_converged, "padj"].isna().all()
