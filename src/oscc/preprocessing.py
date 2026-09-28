"""Independent abundance filtering and explicit gene representative selection."""
import numpy as np
import pandas as pd


def filter_counts(counts, min_count=10, min_samples=3):
    if min_count < 1 or not 1 <= min_samples <= counts.shape[1]:
        raise ValueError("Invalid expression-filter threshold")
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    audit = pd.DataFrame({"samples_ge_min_count": (counts >= min_count).sum(axis=1), "retained": keep})
    if not keep.any():
        raise ValueError("No features pass the expression filter")
    return counts.loc[keep].copy(), audit


def select_representatives(frame):
    """Choose maximum exon count, then lexical RefSeq ID, independent of DE."""
    frame = frame.copy()
    valid = frame.gene.notna() & frame.gene.fillna("").str.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
    valid &= frame.refseq.fillna("").str.fullmatch(r"N[MR]_\d+(?:\.\d+)?")
    frame["mapping_status"] = np.where(valid, "eligible", "missing_or_ambiguous_identifier")
    candidates = frame.loc[valid].sort_values(["exons", "refseq"], ascending=[False, True], kind="stable")
    if candidates.refseq.duplicated().any():
        raise ValueError("Duplicate RefSeq identifiers in source table")
    chosen = candidates.drop_duplicates("gene")
    frame.loc[valid, "mapping_status"] = "redundant_transcript"
    frame.loc[chosen.index, "mapping_status"] = "selected"
    return chosen.sort_values("gene").set_index("gene"), frame


def log_expression(normalized):
    """Visualization only: log2(size-factor-normalized count + 1), not a VST."""
    return np.log2(normalized + 1)
