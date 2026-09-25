"""Retrieve the original publication supplement and inspect GEO SOFT metadata."""
import gzip
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

SOURCES = {
    "TableS1.xls": {
        "url": "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0009317.s009&type=supplementary",
        "sha256": "5aa416fc86ebf716252f7df55e34c8d33dbe035b92889cfdfd879332eb289ffd",
    },
    "GSE20116_family.soft.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20116/soft/GSE20116_family.soft.gz",
        "sha256": "7bfb1cdc0299e2529eb2a79ab22d519201aa64d17fb5d1844ff039fabf057af0",
    },
}
SAMPLE_COLUMNS = ["8N", "8T", "33N", "33T", "51N", "51T"]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, target, expected=None, offline=False, max_bytes=30_000_000):
    """Bounded, atomic download; verify cached bytes as well as new downloads."""
    target = Path(target)
    if not target.exists():
        if offline:
            raise FileNotFoundError(f"Offline source missing: {target}. Run oscc download online first.")
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_suffix(target.suffix + ".part")
        try:
            with requests.get(url, stream=True, timeout=(15, 90)) as response:
                response.raise_for_status()
                size = 0
                with partial.open("wb") as handle:
                    for chunk in response.iter_content(65536):
                        size += len(chunk)
                        if size > max_bytes:
                            raise ValueError(f"Source exceeds {max_bytes} byte download limit")
                        handle.write(chunk)
            if expected and sha256(partial) != expected:
                raise ValueError(f"Checksum mismatch: {url}; inspect upstream changes before updating source lock")
            partial.replace(target)
        finally:
            partial.unlink(missing_ok=True)
    digest = sha256(target)
    if expected and digest != expected:
        raise ValueError(f"Checksum mismatch for cached file {target}")
    return {"url": url, "sha256": digest, "bytes": target.stat().st_size}


def acquire(root, offline=False, dataset="GSE20116"):
    from .datasets import get_dataset
    from .geo_validation import SOURCES as VALIDATION_SOURCES
    raw = get_dataset(dataset).raw(root)
    sources = SOURCES if dataset == "GSE20116" else VALIDATION_SOURCES
    manifest = {name: download(spec["url"], raw / name, spec["sha256"], offline) for name, spec in sources.items()}
    manifest["verified_utc"] = datetime.now(timezone.utc).isoformat()
    return manifest


def read_metadata(path):
    text = gzip.open(path, "rt", encoding="utf-8").read()
    if "MAX format" not in text or "GSE20116_RAW.tar" not in text:
        raise ValueError("Unexpected GEO processing/supplement metadata")
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        accession = block.splitlines()[0].strip()
        title = re.search(r"!Sample_title = (.+)", block).group(1)
        match = re.fullmatch(r"(normal|tumor)_tissue_patient_(8|33|51)", title)
        if not match:
            raise ValueError(f"Unexpected sample title {title}")
        condition, patient = match.groups()
        characteristic = re.search(r"!Sample_characteristics_ch1 = patient: (\d+)", block).group(1)
        if characteristic != patient:
            raise ValueError("Conflicting patient annotations")
        rows.append({"sample_accession": accession, "patient_id": patient, "condition": condition,
                     "source_column": patient + ("N" if condition == "normal" else "T"),
                     "tissue": re.search(r"!Sample_source_name_ch1 = (.+)", block).group(1),
                     "platform": "GPL9442", "instrument": "AB SOLiD System 3.0", "genome": "hg18",
                     "age": pd.NA, "sex": pd.NA, "stage": pd.NA, "hpv_status": pd.NA,
                     "clinical_note": "Age, sex, stage and HPV not supplied in GEO SOFT; no imputation",
                     "sample_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"})
    metadata = pd.DataFrame(rows).set_index("sample_accession")
    expected = [f"GSM{i}" for i in range(515513, 515519)]
    if list(metadata.index) != expected or list(metadata.source_column) != SAMPLE_COLUMNS:
        raise ValueError("Unexpected GSE20116 sample identities/order")
    return metadata


def read_publication_counts(path):
    """Explicitly select Excel D:I raw 'sum' counts, never J:O normalized_sum."""
    table = pd.read_excel(path, sheet_name="counts.annotated.byTranscript", header=None, engine="xlrd")
    if table.iloc[0, 3:9].tolist() != SAMPLE_COLUMNS or table.iloc[2, 3:9].tolist() != ["sum"] * 6:
        raise ValueError("Publication raw-count headers changed")
    if table.iloc[2, :3].tolist() != ["idRefSeq", "nameOfGene", "numberOfExons"]:
        raise ValueError("Publication annotation headers changed")
    frame = table.iloc[3:, :9].copy()
    frame.columns = ["refseq", "gene", "exons"] + SAMPLE_COLUMNS
    frame = frame.dropna(how="all")
    for col in ["exons"] + SAMPLE_COLUMNS:
        frame[col] = pd.to_numeric(frame[col], errors="raise")
    if frame.exons.isna().any() or (frame.exons <= 0).any():
        raise ValueError("Invalid exon annotations")
    return frame
