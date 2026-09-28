"""Strict parsers for the inspected GSE184616 deposited count matrix and SOFT."""
import csv
import gzip
import re
import pandas as pd
from .validation import validate_counts, validate_metadata

SOURCES = {
    "GSE184616_unnormalisedGeneCounts.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184616/suppl/GSE184616_unnormalisedGeneCounts.txt.gz",
        "sha256": "5d9e406952eb5c64e87774b336d5ed58fffca2fe83ead513e485df4e6c071451",
    },
    "GSE184616_family.soft.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184616/soft/GSE184616_family.soft.gz",
        "sha256": "31111ea040a0c3324c73ba4be2438d6c663e05922f204be359f2d768a6eb108c",
    },
}


def read_geo_metadata(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    if "^SERIES = GSE184616" not in text:
        raise ValueError("Expected GSE184616 SOFT")
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        accession = block.splitlines()[0].strip()
        fields = {}
        for line in block.splitlines():
            if line.startswith("!Sample_characteristics_ch1 = "):
                key, value = line.split(" = ", 1)[1].split(": ", 1)
                if key in fields:
                    raise ValueError(f"Duplicate characteristic {key}")
                fields[key] = value
        title = re.search(r"^!Sample_title = (.+)$", block, re.M).group(1)
        condition = {"Adjacent Normal": "normal", "Primary Tumour": "tumor"}.get(fields.get("condition"))
        patient = fields.get("patient id")
        if condition is None or title != f"{patient}-{'N' if condition == 'normal' else 'P'}":
            raise ValueError("Conflicting patient/title/condition annotations")
        rows.append({"sample_accession": accession, "source_column": title,
                     "patient_id": patient, "condition": condition,
                     "tissue": re.search(r"^!Sample_source_name_ch1 = (.+)$", block, re.M).group(1),
                     "age": int(fields["age"]), "sex": fields["gender"], "smoking": fields["smoking"],
                     "diagnosis": fields["patient diagnosis"],
                     "hpv_status": "negative" if "HPV-negative" in fields["patient diagnosis"] else "not stated in sample",
                     "platform": "GPL24676", "genome": "hg38", "annotation": "GENCODE 31"})
    metadata = validate_metadata(pd.DataFrame(rows).set_index("sample_accession"))
    if set(metadata.index) != {f"GSM{i}" for i in range(5593747, 5593777)}:
        raise ValueError("Unexpected deposited sample accessions")
    if set(metadata.patient_id) != {f"OSCC_{i}" for i in [*range(1, 15), 16]}:
        raise ValueError("Unexpected deposited patient identities")
    for column in ["age", "sex", "smoking", "tissue"]:
        if (metadata.groupby("patient_id")[column].nunique() != 1).any():
            raise ValueError(f"Conflicting within-patient {column}")
    return metadata


def read_geo_counts(path, metadata):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if header[0] != "Gene" or len(header) != len(set(header)):
        raise ValueError("Expected unique Gene/sample columns")
    counts = pd.read_csv(path, sep="\t", index_col=0, keep_default_na=False)
    if set(counts.columns) != set(metadata.source_column):
        raise ValueError("Raw count headers do not match GEO sample titles")
    counts = counts.rename(columns=dict(zip(metadata.source_column, metadata.index)))
    counts.index.name = "gene"
    return validate_counts(counts, metadata)[0]
