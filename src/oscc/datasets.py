"""Explicit cohort identity, storage, and independently fixed analysis settings."""
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Dataset:
    accession: str
    role: str
    min_samples: int
    description: str

    def output(self, root):
        return Path(root) / "results" / self.accession

    def raw(self, root):
        base = Path(root) / "data/raw"
        return base if self.role == "discovery" else base / self.accession


DATASETS = {
    "GSE20116": Dataset("GSE20116", "discovery", 3,
        "Three matched pairs; historical SOLiD/hg18 RefSeq representatives."),
    "GSE184616": Dataset("GSE184616", "validation", 15,
        "Fifteen matched pairs; NovaSeq/hg38 gene counts; HPV-negative series."),
}


def get_dataset(accession):
    try:
        return DATASETS[accession]
    except KeyError:
        raise ValueError(f"Unsupported dataset: {accession}") from None


def bundle_directory(root, accession="GSE20116", legacy=False):
    """Legacy discovery snapshot is read-only; new writes always use output()."""
    dataset = get_dataset(accession)
    original = Path(root) / "results"
    if accession == "GSE20116" and (legacy or not (dataset.output(root) / "manifest.json").exists()):
        return original
    if legacy:
        raise ValueError("Only GSE20116 has an original legacy snapshot")
    return dataset.output(root)
