from pathlib import Path
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "oscc-test-mpl"))
import pandas as pd
import pytest


@pytest.fixture
def fixture_data():
    folder = Path(__file__).parent / "fixtures"
    return pd.read_csv(folder / "counts.csv", index_col=0), pd.read_csv(folder / "metadata.csv", index_col=0)
