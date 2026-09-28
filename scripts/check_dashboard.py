"""Exercise real dashboard widgets; requires a completed real-data pipeline."""
from pathlib import Path
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
app = AppTest.from_file(str(ROOT / "app/streamlit_app.py"), default_timeout=60).run()
assert not app.exception, [e.message for e in app.exception]
assert len(app.metric) == 4
app.slider(key="alpha").set_value(.01).run()
assert not app.exception
app.slider(key="lfc").set_value(2.).run()
assert not app.exception
app.text_input[0].set_value("MMP").run()
assert not app.exception
app.checkbox[0].set_value(True).run()
assert not app.exception
app.selectbox(key="gene-selector").set_value("MMP1").run()
assert not app.exception
app.radio(key="method").set_value("Overrepresentation").run()
assert not app.exception
app.text_input[1].set_value("Myogenesis").run()
assert not app.exception
print("Dashboard initial render, threshold sliders, searches, all-genes table, gene selection and ORA explorer passed")

app.radio(key="dataset").set_value("GSE184616").run()
assert not app.exception, [e.message for e in app.exception]
assert "30 / 15" in [m.value for m in app.metric]
print("Validation cohort render and sample/pair metrics passed")
