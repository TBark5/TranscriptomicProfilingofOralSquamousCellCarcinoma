import numpy as np
import pytest
from oscc.validation import validate_counts, validate_metadata, design_matrix
from oscc.preprocessing import filter_counts, select_representatives
import pandas as pd


def test_pairing_and_rank(fixture_data):
    counts, metadata = fixture_data
    assert design_matrix(metadata).shape == (6, 4)
    assert validate_counts(counts, metadata)[0].shape == counts.shape


def test_missing_pair(fixture_data):
    _, metadata = fixture_data
    with pytest.raises(ValueError, match="exactly one"):
        validate_metadata(metadata.iloc[:-1])


def test_duplicate_pair_and_missing_patient(fixture_data):
    _, metadata = fixture_data
    metadata.iloc[1, metadata.columns.get_loc("condition")] = "normal"
    with pytest.raises(ValueError):
        validate_metadata(metadata)
    metadata.iloc[0, metadata.columns.get_loc("patient_id")] = None
    with pytest.raises(ValueError, match="Missing"):
        validate_metadata(metadata)


def test_metadata_alignment(fixture_data):
    counts, metadata = fixture_data
    aligned, _ = validate_counts(counts.iloc[:, ::-1], metadata)
    pd.testing.assert_frame_equal(aligned, counts, check_names=False)
    with pytest.raises(ValueError, match="identifiers differ"):
        validate_counts(counts.rename(columns={counts.columns[0]: "WRONG"}), metadata)


@pytest.mark.parametrize("value", [-1, .5, np.nan, np.inf])
def test_invalid_counts(fixture_data, value):
    counts, metadata = fixture_data
    counts = counts.astype(float)
    counts.iloc[0, 0] = value
    with pytest.raises(ValueError):
        validate_counts(counts, metadata)


def test_duplicate_features_and_nonnumeric(fixture_data):
    counts, metadata = fixture_data
    with pytest.raises(ValueError, match="duplicate"):
        validate_counts(pd.concat([counts, counts.iloc[:1]]), metadata)
    with pytest.raises(ValueError, match="numeric"):
        validate_counts(counts.astype(str), metadata)


def test_filter_boundary(fixture_data):
    counts, _ = fixture_data
    retained, audit = filter_counts(counts, 10, 3)
    assert "SYN_BOUNDARY" in retained.index
    assert "SYN_LOW" not in retained.index
    assert not audit.loc["SYN_ZERO", "retained"]
    with pytest.raises(ValueError, match="No features"):
        filter_counts(counts, 10000, 3)


def test_representatives_do_not_sum_isoforms():
    data = pd.DataFrame({"gene": ["A", "A", "B", None], "refseq": ["NM_002", "NM_001", "NM_003", "NM_004"], "exons": [5, 5, 4, 4], "sample": [20, 10, 30, 5]})
    selected, audit = select_representatives(data)
    assert selected.loc["A", "sample"] == 10
    assert selected.loc["A", "refseq"] == "NM_001"
    assert len(selected) == 2
    assert (audit.mapping_status == "missing_or_ambiguous_identifier").sum() == 1
