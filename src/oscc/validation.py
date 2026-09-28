"""Fail early on invalid counts, identifiers, or matched designs."""
import numpy as np
import pandas as pd


def validate_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    required = {"patient_id", "condition"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Metadata requires {required}")
    if metadata.empty or metadata.index.has_duplicates or metadata.index.isna().any():
        raise ValueError("Sample identifiers must be present and unique")
    out = metadata.copy()
    if out[list(required)].isna().any().any():
        raise ValueError("Missing patient or condition")
    out["patient_id"] = out.patient_id.astype(str)
    if not set(out.condition).issubset({"normal", "tumor"}):
        raise ValueError("Conditions must be normal or tumor")
    if any(not str(x).strip() for x in out.index) or out.patient_id.str.strip().eq("").any():
        raise ValueError("Blank identifiers")
    pairs = pd.crosstab(out.patient_id, out.condition).reindex(columns=["normal", "tumor"], fill_value=0)
    if len(pairs) < 2 or not (pairs == 1).all().all():
        raise ValueError("Require at least two patients, each with exactly one normal and one tumor")
    return out


def validate_counts(counts: pd.DataFrame, metadata: pd.DataFrame) -> tuple:
    """Return genes x samples integer counts in metadata order; never round inputs."""
    metadata = validate_metadata(metadata)
    if counts.empty or counts.index.has_duplicates or counts.columns.has_duplicates:
        raise ValueError("Empty matrix or duplicate gene/sample identifiers")
    if counts.index.isna().any() or any(not str(x).strip() for x in counts.index):
        raise ValueError("Missing gene identifiers")
    if set(counts.columns) != set(metadata.index):
        raise ValueError("Count columns and metadata sample identifiers differ")
    if not all(pd.api.types.is_numeric_dtype(t) for t in counts.dtypes):
        raise ValueError("Counts must have numeric dtypes")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Counts must be finite, nonnegative and complete")
    if not (values == np.floor(values)).all():
        raise ValueError("Raw integer counts required; normalized values cannot be modeled as counts")
    if (values >= np.iinfo(np.int64).max).any() or (counts.sum(axis=0) == 0).any():
        raise ValueError("Invalid count magnitude or empty library")
    return counts.loc[:, metadata.index].astype("int64"), metadata


def design_matrix(metadata: pd.DataFrame) -> pd.DataFrame:
    metadata = validate_metadata(metadata)
    x = pd.get_dummies(metadata[["patient_id", "condition"]], drop_first=True, dtype=float)
    x.insert(0, "Intercept", 1.0)
    if np.linalg.matrix_rank(x) != x.shape[1] or x.shape[0] <= x.shape[1]:
        raise ValueError("Paired design is not identifiable or has no residual degrees of freedom")
    return x
