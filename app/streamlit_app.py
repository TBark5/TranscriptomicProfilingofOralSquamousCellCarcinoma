"""Run from the repository: streamlit run app/streamlit_app.py."""
from html import escape
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from oscc.dashboard import load_results, filtered_results, gene_expression
from oscc.differential import classify
from oscc.datasets import get_dataset, bundle_directory

st.set_page_config(page_title="OSCC Transcriptomics", page_icon="🧬", layout="wide")

# Palette shared by CSS and Plotly: tumor-up is terracotta, tumor-down is deep teal.
INK, MUTED, LINE, PAPER = "#16232E", "#5A6772", "#E3DED4", "#F6F4EF"
UP, DOWN, NEUTRAL = "#C0532F", "#1D5C7A", "#C9CDD0"
DIVERGING = [[0, DOWN], [.5, "#F4F1EA"], [1, UP]]
SEQUENTIAL = [[0, "#123E53"], [.5, "#6F9FB3"], [1, "#E4ECEF"]]
COHORTS = {"GSE20116": "Discovery", "GSE184616": "Validation"}

st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,400;9..144,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');
html, body, [class*="css"], .stApp, button, input, textarea, select {font-family: 'IBM Plex Sans', system-ui, sans-serif;}
.stApp {background: #F6F4EF;}
.block-container {padding-top: 2.5rem; padding-bottom: 3rem; max-width: 1400px;}
h1, h2, h3, .stApp h1, .stApp h2, .stApp h3 {font-family: 'Fraunces', Georgia, serif !important; letter-spacing: -.015em; color: #16232E;}
h2, h3 {font-weight: 600 !important;}
code, pre {font-family: 'IBM Plex Mono', monospace;}
[data-testid="stHeader"] {background: transparent;}

/* Hero */
.oscc-chips {display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 14px;}
.oscc-chip {padding: 5px 12px; border-radius: 999px; background: #EFEAE0; color: #5A4E3A; font-size: 12px; font-weight: 600; letter-spacing: .06em; text-transform: uppercase;}
.oscc-chip.accent {background: #E6EEF1; color: #1D5C7A;}
.oscc-title, .stApp h1.oscc-title {font-family: 'Fraunces', Georgia, serif !important; font-weight: 400 !important; padding: 0 !important; font-size: clamp(2rem, 3.4vw, 2.9rem); line-height: 1.08; letter-spacing: -.02em; margin: 0 0 1rem; color: #16232E; max-width: 900px;}
.oscc-lede {font-size: 1.06rem; line-height: 1.55; color: #4A5763; max-width: 760px; margin: 0 0 1.25rem;}
.oscc-footer {margin-top: 2.5rem; padding-top: 1rem; border-top: 1px solid #E3DED4; font-size: 12px; color: #6B7780;}

/* Banners */
.oscc-banner {display: flex; gap: 14px; align-items: flex-start; padding: 16px 20px; border-radius: 14px; font-size: 14px; line-height: 1.55; margin-bottom: 12px;}
.oscc-banner svg {flex-shrink: 0; margin-top: 1px;}
.oscc-banner.caution {background: #FBF1DE; border: 1px solid #EBD6AE; color: #5C4515;}
.oscc-banner.note {background: #E9F0F3; border: 1px solid #C9D9E0; color: #123E53;}

/* Cards: every bordered container is a white card */
[data-testid="stVerticalBlockBorderWrapper"]:has(> [class*="st-key-card"]) {background: #FFFFFF; border: 1px solid #E3DED4 !important; border-radius: 18px !important;}
.oscc-card-head {display: flex; flex-direction: column; gap: 4px; margin-bottom: 4px;}
.oscc-card-head h2 {margin: 0; padding: 0; font-size: 22px;}
.oscc-card-head span {font-size: 13px; color: #5A6772;}

/* Metric cards */
[data-testid="stMetric"] {padding: 4px 6px 0;}
[data-testid="stMetricLabel"] p {font-size: 13px; color: #5A6772;}
[data-testid="stMetricValue"] {font-family: 'Fraunces', Georgia, serif; font-size: 2.4rem; color: #16232E;}
.oscc-kpi-2 [data-testid="stMetricValue"] {color: #9E3F1F;}
.oscc-kpi-3 [data-testid="stMetricValue"] {color: #1D5C7A;}
[data-testid="stColumn"]:nth-of-type(2) [data-testid="stMetricValue"] {color: #9E3F1F;}
[data-testid="stColumn"]:nth-of-type(3) [data-testid="stMetricValue"] {color: #1D5C7A;}
[data-testid="stColumn"]:nth-of-type(2) [data-testid="stMetricLabel"] p::before,
[data-testid="stColumn"]:nth-of-type(3) [data-testid="stMetricLabel"] p::before {content: ""; display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 8px; vertical-align: 0;}
[data-testid="stColumn"]:nth-of-type(2) [data-testid="stMetricLabel"] p::before {background: #C0532F;}
[data-testid="stColumn"]:nth-of-type(3) [data-testid="stMetricLabel"] p::before {background: #1D5C7A;}
.oscc-kpi-foot {font-size: 12px; color: #5A6772; padding: 0 6px 4px;}
.oscc-meter {height: 6px; border-radius: 3px; margin: 6px 6px 6px;}
.oscc-meter > div {height: 6px; border-radius: 3px;}

/* Tabs as a segmented control */
.stTabs [data-baseweb="tab-list"] {gap: 4px; padding: 4px; background: #ECE8E0; border-radius: 12px; width: fit-content; max-width: 100%; overflow-x: auto;}
.stTabs [data-baseweb="tab"] {height: 44px; padding: 0 18px; border-radius: 9px; color: #4A5763; background: transparent;}
.stTabs [aria-selected="true"] {background: #FFFFFF; color: #16232E; box-shadow: 0 1px 2px rgba(20,33,43,.08);}
.stTabs [aria-selected="true"] p {font-weight: 600;}
.stTabs [data-baseweb="tab-highlight"], .stTabs [data-baseweb="tab-border"] {display: none;}
.stTabs [data-baseweb="tab-panel"] {padding-top: 1.5rem;}

/* Bar lists (leading genes, pathways) */
.oscc-group {font-size: 12px; font-weight: 600; letter-spacing: .08em; text-transform: uppercase; margin: 6px 0 10px;}
.oscc-bars {display: flex; flex-direction: column; gap: 10px; margin-bottom: 14px;}
.oscc-bar-row {display: grid; grid-template-columns: 96px minmax(0, 1fr) 64px; gap: 12px; align-items: center; font-size: 14px;}
.oscc-bar-row b {font-weight: 600; overflow: hidden; text-overflow: ellipsis;}
.oscc-track {height: 8px; border-radius: 4px;}
.oscc-track > div {height: 8px; border-radius: 4px;}
.oscc-num {font-family: 'IBM Plex Mono', monospace; font-size: 13px; text-align: right; color: #4A5763;}
.oscc-rule {height: 1px; background: #EFECE6; margin: 4px 0 14px;}
.oscc-note {font-size: 12px; color: #5A6772;}
.oscc-path-row {display: grid; grid-template-columns: minmax(120px, 190px) minmax(0, 1fr) minmax(0, 1fr) 56px; align-items: center; font-size: 13px; margin-bottom: 9px;}
.oscc-path-row > span:first-child {padding-right: 12px;}
.oscc-half {display: flex; height: 20px; align-items: center;}
.oscc-half.neg {justify-content: flex-end; border-right: 1px solid #B9BEC2;}
.oscc-half > div {height: 14px;}
.oscc-half.neg > div {border-radius: 4px 0 0 4px; background: #1D5C7A;}
.oscc-half.pos > div {border-radius: 0 4px 4px 0; background: #C0532F;}

/* Provenance */
.oscc-dl {display: grid; grid-template-columns: 120px minmax(0, 1fr); row-gap: 12px; column-gap: 12px; font-size: 14px; margin: 8px 0 18px;}
.oscc-dl dt {color: #5A6772;}
.oscc-dl dd {margin: 0;}
.oscc-link {display: flex; align-items: center; justify-content: center; gap: 8px; min-height: 44px; border-radius: 10px; border: 1px solid #C9D6DC; font-size: 14px; font-weight: 600; color: #1D5C7A !important; text-decoration: none !important;}
.oscc-link:hover {background: #E9F0F3;}

[data-testid="stAlert"] {border-radius: 14px;}
[data-testid="stDataFrame"] {border: 1px solid #E3DED4; border-radius: 12px; overflow: hidden;}
[data-testid="stImage"] img {border: 1px solid #E3DED4; border-radius: 16px; background: #FFFFFF;}
[data-testid="stExpander"] details {background: #FFFFFF; border: 1px solid #E3DED4; border-radius: 14px;}
.stDownloadButton button, .stButton button {border-radius: 10px; min-height: 44px; border: 1px solid #C9D6DC; font-weight: 600;}

/* Sidebar */
[data-testid="stSidebar"] {background: #14212B;}
[data-testid="stSidebar"] * {color: #E9E5DC;}
[data-testid="stSidebar"] h2 {font-size: 12px !important; font-family: 'IBM Plex Sans', sans-serif !important; font-weight: 600 !important; letter-spacing: .1em; text-transform: uppercase; color: #A9B4BC !important; margin-top: 1.25rem; padding-bottom: .25rem;}
[data-testid="stSidebar"] [data-testid="stCaptionContainer"] p {color: #A9B4BC; line-height: 1.55;}
[data-testid="stSidebar"] [data-testid="stSliderThumbValue"], [data-testid="stSidebar"] [data-testid="stSliderTickBarMin"], [data-testid="stSidebar"] [data-testid="stSliderTickBarMax"] {color: #E9A27F; font-family: 'IBM Plex Mono', monospace;}
[data-testid="stSidebar"] [role="slider"] {background: #E9A27F; border-color: #E9A27F;}
/* Cohort radio drawn as a two-button segmented control */
[data-testid="stSidebar"] [role="radiogroup"] {display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 4px; padding: 4px; background: #1F303D; border-radius: 12px;}
[data-testid="stSidebar"] [role="radiogroup"] > label {margin: 0; padding: 8px 4px; min-height: 56px; border-radius: 9px; justify-content: center; align-items: center; text-align: center; cursor: pointer;}
[data-testid="stSidebar"] [role="radiogroup"] > label > div:first-child {display: none;}
[data-testid="stSidebar"] [role="radiogroup"] > label p {font-weight: 600; font-size: 14px; color: #CBD3D9;}
[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) {background: #F6F4EF;}
[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) p {color: #14212B;}
[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) [data-testid="stCaptionContainer"] p {color: #5A6772;}
[data-testid="stSidebar"] [role="radiogroup"] [data-testid="stCaptionContainer"] p {font-size: 11px; font-weight: 500; color: #A9B4BC;}
.oscc-brand {display: flex; align-items: center; gap: 12px; margin: .25rem 0 1rem;}
.oscc-brand b {font-family: 'Fraunces', Georgia, serif; font-size: 20px; font-weight: 600; display: block;}
.oscc-brand span {font-size: 11px; letter-spacing: .08em; text-transform: uppercase; color: #A9B4BC;}
.oscc-model {margin-top: 1.5rem; padding: 16px; border-radius: 12px; background: #1F303D; display: flex; flex-direction: column; gap: 6px; font-size: 13px;}
.oscc-model span:first-child {font-size: 12px; letter-spacing: .08em; text-transform: uppercase; color: #A9B4BC;}
.oscc-model code {background: transparent; color: #E9E5DC; padding: 0; font-size: 13px;}
</style>""", unsafe_allow_html=True)

CAUTION_ICON = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="#8A6414" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M12 3 2 20h20L12 3z"/><path d="M12 10v4M12 17h.01"/></svg>'
INFO_ICON = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="#1D5C7A" stroke-width="2" stroke-linecap="round" aria-hidden="true"><circle cx="12" cy="12" r="9"/><path d="M12 11v6M12 7.5h.01"/></svg>'
CHECK_ICON = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#2F6B45" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" style="vertical-align:-3px;margin-right:6px"><path d="m5 12 5 5L20 7"/></svg>'


def style_figure(fig, height=None):
    """Apply the dashboard's typography, grid and palette to a Plotly figure."""
    fig.update_layout(template="plotly_white", paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font={"family": "IBM Plex Sans, system-ui, sans-serif", "color": INK, "size": 13},
                      margin={"l": 8, "r": 8, "t": 40, "b": 8},
                      legend={"title": None, "orientation": "h", "y": 1.02, "yanchor": "bottom", "x": 1, "xanchor": "right"},
                      hoverlabel={"bgcolor": "#FFFFFF", "bordercolor": LINE, "font_family": "IBM Plex Sans, sans-serif"})
    fig.update_xaxes(gridcolor="#EFECE6", zerolinecolor="#B9BEC2", linecolor=LINE)
    fig.update_yaxes(gridcolor="#EFECE6", zerolinecolor="#B9BEC2", linecolor=LINE)
    if fig.layout.title.text:
        fig.update_layout(title_font={"family": "Fraunces, Georgia, serif", "size": 20}, margin={"t": 64})
    if height:
        fig.update_layout(height=height)
    return fig


def card_head(title, subtitle=""):
    st.markdown(f'<div class="oscc-card-head"><h2>{escape(title)}</h2>{f"<span>{escape(subtitle)}</span>" if subtitle else ""}</div>', unsafe_allow_html=True)


def signed(value):
    return f"{'+' if value > 0 else '−'}{abs(value):.2f}"


def bar_rows(frame, color, track):
    scale = max(8., frame.log2FoldChange.abs().max()) if len(frame) else 8.
    rows = "".join(f'<div class="oscc-bar-row"><b>{escape(str(gene))}</b><div class="oscc-track" style="background:{track}"><div style="width:{abs(row.log2FoldChange) / scale * 100:.0f}%;background:{color}"></div></div><span class="oscc-num">{signed(row.log2FoldChange)}</span></div>'
                   for gene, row in frame.iterrows())
    return f'<div class="oscc-bars">{rows or "<span class=oscc-note>None pass the current thresholds.</span>"}</div>'


# ---------- Sidebar: brand and cohort ----------
st.sidebar.markdown("""<div class="oscc-brand"><svg width="36" height="36" viewBox="0 0 36 36" fill="none" stroke="#E9A27F" stroke-width="2" stroke-linecap="round" aria-hidden="true"><path d="M12 4c0 9 12 9 12 18s-12 9-12 10"/><path d="M24 4c0 9-12 9-12 18s12 9 12 10"/><path d="M14 10h8M13.5 26h9M16 18h4"/></svg><div><b>OSCC Explorer</b><span>Paired bulk RNA-seq</span></div></div>""", unsafe_allow_html=True)
st.sidebar.header("Cohort")
dataset = st.sidebar.radio("Cohort", list(COHORTS), captions=list(COHORTS.values()), key="dataset", label_visibility="collapsed")
cohort = get_dataset(dataset)
output = bundle_directory(ROOT, dataset)


@st.cache_data(show_spinner="Validating analysis outputs…")
def cached_load(root, signature, accession):
    return load_results(root, accession)


try:
    # Changes to any required output invalidate cached validation.
    paths = [output / "manifest.json"] + sorted(p for p in output.rglob("*") if p.is_file()) + ([ROOT / "VALIDATION_REPORT.md"] if dataset == "GSE184616" else [])
    signature = tuple((str(p), p.stat().st_mtime_ns, p.stat().st_size) for p in paths)
    bundle = cached_load(str(ROOT), signature, dataset)
except (OSError, ValueError, KeyError) as exc:
    st.error(str(exc))
    st.code("python -m oscc.cli run\nstreamlit run app/streamlit_app.py")
    st.info("This app requires genuine, completed pipeline outputs. It never substitutes demo measurements.")
    st.stop()

manifest = bundle["manifest"]
results = bundle["differential_expression"]
metadata = bundle["metadata"]
diagnostics = manifest["diagnostics"]
patients = list(dict.fromkeys(metadata.patient_id.astype(str)))

# ---------- Sidebar: thresholds and model ----------
st.sidebar.header("Gene thresholds")
alpha = st.sidebar.slider("BH-adjusted p-value cutoff", .001, .20, float(max(.001, min(.20, manifest["config"]["alpha"]))), .001, format="%.3f", key="alpha")
lfc = st.sidebar.slider("Absolute log2 fold-change cutoff", 0., 5., float(min(manifest["config"]["lfc"], 5)), .1, key="lfc")
st.sidebar.caption("Strict inequalities. LFC is tumor / normal. These controls filter the fitted results; they do not refit the model or recalculate pathways.")
st.sidebar.markdown(f"""<div class="oscc-model"><span>Model</span><code>{escape(diagnostics.get("design", "~ patient_id + condition"))}</code><span style="color:#A9B4BC">PyDESeq2 · Wald · BH</span></div>""", unsafe_allow_html=True)
calls = classify(results, alpha, lfc)
n_up, n_down = int(sum(calls == "up")), int(sum(calls == "down"))

# ---------- Header ----------
st.markdown(f"""<div class="oscc-chips"><span class="oscc-chip accent">{dataset}</span><span class="oscc-chip">{cohort.role} cohort</span><span class="oscc-chip">{len(patients)} matched patients</span></div>
<h1 class="oscc-title">Transcriptomic profiling of oral squamous cell carcinoma</h1>
<p class="oscc-lede">Explore tumor–normal expression differences while preserving each patient's matched comparison.</p>""", unsafe_allow_html=True)
st.markdown(f'<div class="oscc-banner caution">{CAUTION_ICON}<span>{escape(cohort.description)} {diagnostics["residual_df"]} residual degrees of freedom. These are exploratory associations, not validated biomarkers or clinical claims.</span></div>', unsafe_allow_html=True)
if dataset == "GSE184616":
    summary = manifest["summary"]
    st.markdown(f'<div class="oscc-banner note">{INFO_ICON}<span>Independent validation: {summary["replicated_candidates"]:,} of {summary["testable_candidates"]:,} testable discovery candidates replicated; {summary["replicated_pathways"]} discovery pathways replicated. See VALIDATION_REPORT.md for criteria and limitations.</span></div>', unsafe_allow_html=True)

# ---------- KPI cards ----------
if "source_transcripts" in manifest:
    gene_note = f"from {manifest['source_transcripts']:,} RefSeq transcripts"
elif "source_genes" in manifest:
    gene_note = f"from {manifest['source_genes']:,} source genes"
else:
    gene_note = f"{manifest['summary']['tested']:,} tested"
patient_note = f"patients {', '.join(patients[:-1])} and {patients[-1]}" if 1 < len(patients) <= 4 else f"{len(patients)} tumor–normal pairs"
share_up = n_up / max(n_up + n_down, 1)
kpis = [("Gene representatives" if dataset == "GSE20116" else "Genes", f"{len(results):,}", f'<div class="oscc-kpi-foot">{gene_note}</div>'),
        ("Upregulated in tumor", f"{n_up:,}", f'<div class="oscc-meter" style="background:#F3E4DC"><div style="width:{share_up:.0%};background:{UP}"></div></div>'),
        ("Downregulated in tumor", f"{n_down:,}", f'<div class="oscc-meter" style="background:#DDE8ED"><div style="width:{1 - share_up:.0%};background:{DOWN}"></div></div>'),
        ("Samples / pairs", f"{len(metadata)} / {metadata.patient_id.nunique()}", f'<div class="oscc-kpi-foot">{escape(patient_note)}</div>')]
for column, (label, value, foot) in zip(st.columns(4), kpis):
    with column.container(border=True, key=f"card-kpi-{label.split()[0].lower()}"):
        st.metric(label, value)
        st.markdown(foot, unsafe_allow_html=True)

genes_tab, overview, expression_tab, pathways_tab = st.tabs(["Differential expression", "Overview & provenance", "Gene & heatmap explorer", "Pathways"])

with genes_tab:
    left, right = st.columns([1.65, 1], gap="medium")
    with left.container(border=True, key="card-2"):
        card_head("Volcano plot", "log2 fold change (tumor / normal) against −log10 BH-adjusted p")
        plot = results.reset_index().copy()
        plot["direction"] = calls.to_numpy()
        plot["neglog10_q"] = -np.log10(plot.padj.clip(lower=1e-300))
        fig = px.scatter(plot, x="log2FoldChange", y="neglog10_q", color="direction", hover_name="gene",
                         hover_data=[c for c in ["refseq", "pvalue", "padj", "test_status"] if c in plot.columns],
                         color_discrete_map={"up": UP, "down": DOWN, "not_significant": NEUTRAL},
                         category_orders={"direction": ["up", "down", "not_significant"]},
                         labels={"log2FoldChange": "log2 fold change · tumor / normal", "neglog10_q": "−log10(BH-adjusted p)"})
        fig.add_hline(y=-np.log10(alpha), line_dash="dash", line_color="#8C959C", line_width=1)
        for x in [-lfc, lfc]:
            fig.add_vline(x=x, line_dash="dash", line_color="#8C959C", line_width=1)
        fig.update_traces(marker={"size": 6, "opacity": .75, "line": {"width": 0}})
        legend_names = {"up": "Up", "down": "Down", "not_significant": "Not significant"}
        fig.for_each_trace(lambda trace: trace.update(name=legend_names.get(trace.name, trace.name)))
        style_figure(fig, 440)
        st.plotly_chart(fig, use_container_width=True)
    with right.container(border=True, key="card-3"):
        card_head("Leading genes")
        significant = results.assign(direction=calls).dropna(subset=["padj"]).sort_values("padj")
        top_up = significant.loc[significant.direction == "up"].head(5)
        top_down = significant.loc[significant.direction == "down"].head(5)
        st.markdown(f'<div class="oscc-group" style="color:#9E3F1F">Up in tumor</div>{bar_rows(top_up, UP, "#F6EDE8")}<div class="oscc-rule"></div>'
                    f'<div class="oscc-group" style="color:{DOWN}">Down in tumor</div>{bar_rows(top_down, DOWN, "#E8EFF2")}'
                    '<span class="oscc-note">Ranked within direction by adjusted p · bar = |log2 FC|</span>', unsafe_allow_html=True)

    left, right = st.columns([1.65, 1], gap="medium")
    with left.container(border=True, key="card-4"):
        gsea = bundle["gsea"]
        libraries = sorted(gsea.library.unique()) if not gsea.empty else []
        library = next((lib for lib in libraries if "Hallmark" in lib), libraries[0] if libraries else None)
        card_head(f"{'Hallmark' if library and 'Hallmark' in library else 'Top'} pathways",
                  "Preranked GSEA · normalized enrichment score · positive = tumor-upregulated end of the ranking")
        if library is None:
            st.markdown('<span class="oscc-note">No GSEA table is available for this cohort.</span>', unsafe_allow_html=True)
        else:
            ranked = gsea.loc[gsea.library == library].sort_values("NES", ascending=False)
            shown = pd.concat([ranked.loc[ranked.NES > 0].head(5), ranked.loc[ranked.NES < 0].tail(4)])
            scale = max(2.8, shown.NES.abs().max())
            rows = "".join(f'<div class="oscc-path-row"><span>{escape(str(row.Term))}</span>'
                           f'<div class="oscc-half neg"><div style="width:{max(-row.NES, 0) / scale * 100:.0f}%"></div></div>'
                           f'<div class="oscc-half pos"><div style="width:{max(row.NES, 0) / scale * 100:.0f}%"></div></div>'
                           f'<span class="oscc-num">{signed(row.NES)}</span></div>' for row in shown.itertuples())
            st.markdown(rows + '<span class="oscc-note">Full tables, ORA and search are in the Pathways tab.</span>', unsafe_allow_html=True)
    with right.container(border=True, key="card-5"):
        card_head("Provenance")
        platform_parts = [str(metadata[c].iloc[0]) for c in ["instrument", "platform", "genome"] if c in metadata.columns and pd.notna(metadata[c].iloc[0])]
        platform = " · ".join(platform_parts[:1] + platform_parts[-1:]) if platform_parts else "—"
        source = "Tuch et al. 2010, raw count columns" if dataset == "GSE20116" else "GEO deposited unnormalized counts"
        st.markdown(f"""<dl class="oscc-dl">
<dt>Source</dt><dd>{source}</dd>
<dt>Platform</dt><dd>{escape(platform)}</dd>
<dt>Normalization</dt><dd>{escape(diagnostics.get("normalization", "Median-of-ratios size factors"))}</dd>
<dt>Excluded</dt><dd>{diagnostics["failed_convergence"]} features failed convergence</dd>
<dt>Integrity</dt><dd>{CHECK_ICON}Output hashes verified</dd>
</dl>
<a class="oscc-link" href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset}" target="_blank" rel="noopener">Open GEO record
<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" aria-hidden="true"><path d="M7 17 17 7M8 7h9v9"/></svg></a>""", unsafe_allow_html=True)

    st.subheader("Gene results")
    search = st.text_input("Search source gene symbols", placeholder="e.g. MMP")
    show_all = st.checkbox("Include genes that do not pass the thresholds", value=False)
    view = results.copy() if show_all else filtered_results(results, alpha, lfc)
    view["direction"] = classify(view, alpha, lfc)
    if search:
        view = view.loc[view.index.str.contains(search, case=False, regex=False)]
    st.dataframe(view.sort_values("padj"), use_container_width=True)
    st.caption("Missing adjusted p-values are omitted from the volcano; zero values are clipped to 1e−300 for display only.")
    st.download_button("Download displayed gene results", view.to_csv().encode(), "oscc_genes.csv", "text/csv")

with overview:
    st.subheader(f"{dataset}: {cohort.role} cohort")
    st.write(cohort.description)
    if dataset == "GSE20116":
        st.write("Original publication raw counts; one RefSeq representative per historical symbol, selected independently of differential expression.")
    else:
        st.write("Deposited unnormalized gene counts; sample titles and patient characteristics verified against GEO. Exact-symbol comparisons use frozen original discovery candidates.")
    st.dataframe(metadata, use_container_width=True)
    st.write("Model: PyDESeq2 negative-binomial regression, `~ patient_id + condition`, tumor versus normal Wald contrast, median-of-ratios normalization, convergence-qualified BH correction.")
    st.caption(f"{diagnostics['failed_convergence']} features failed an optimizer convergence check and are excluded from significance calls. See the report for assumptions and the complete diagnostics.")
    if (output / "figures/pca.png").exists():
        st.image(str(output / "figures/pca.png"), caption="PCA on up to 2,000 most variable log normalized genes; lines preserve patient pairing.", width=750)
    with st.expander("Source hashes, configuration and model diagnostics"):
        st.json(manifest)

with expression_tab:
    with st.container(border=True, key="card-6"):
        card_head("Patient-matched expression", "log2(normalized count + 1), one line per patient")
        options = results.sort_values("padj").index.tolist()
        gene = st.selectbox("Search or select a gene", options, key="gene-selector")
        frame = gene_expression(bundle, gene)
        fig = px.line(frame, x="condition", y="log_expression", color="patient_id", markers=True,
                      category_orders={"condition": ["normal", "tumor"]}, hover_data=["sample_accession", "normalized_count"],
                      labels={"log_expression": "log2(normalized count + 1)", "patient_id": "Patient", "condition": ""}, title=gene,
                      color_discrete_sequence=["#1D5C7A", "#C0532F", "#6B8E4E", "#8A6BB0", "#B38A2E", "#4A5763"])
        fig.update_traces(line={"width": 2.5}, marker={"size": 10})
        style_figure(fig, 440)
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(results.loc[[gene]], use_container_width=True)
    with st.container(border=True, key="card-7"):
        card_head("Top-gene heatmap", "Row z-scores of the genes passing the sidebar thresholds")
        n = st.slider("Number of heatmap genes", 5, 50, 25, key="heatmap_n")
        top = filtered_results(results, alpha, lfc).head(n).index
        if len(top) < 2:
            st.info("Fewer than two genes pass the selected cutoffs. Showing top ranked genes as an explicitly exploratory fallback.")
            top = results.loc[results.padj.notna()].sort_values("padj").head(n).index
        if len(top):
            matrix = bundle["log_expression"].loc[top]
            matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
            matrix.columns = metadata.patient_id + " " + metadata.condition
            fig = px.imshow(matrix, color_continuous_scale=DIVERGING, zmin=-2, zmax=2, aspect="auto", labels={"color": "Row z-score"})
            style_figure(fig, 650)
            st.plotly_chart(fig, use_container_width=True)
            st.caption("Samples remain in matched-pair order here. The saved publication heatmap uses hierarchical clustering. DE-selected heatmaps are not independent validation.")

with pathways_tab:
    st.markdown(f'<div class="oscc-banner note">{INFO_ICON}<span>Pathways were computed using the saved pipeline thresholds (q &lt; {manifest["config"]["alpha"]}, |LFC| &gt; {manifest["config"]["lfc"]}); sidebar sliders affect gene views only. Positive NES denotes the tumor-upregulated end of the ranking, not pathway activation.</span></div>', unsafe_allow_html=True)
    method = st.radio("Analysis", ["Preranked GSEA", "Overrepresentation"], horizontal=True, key="method")
    table = bundle["gsea" if method == "Preranked GSEA" else "ora"].copy()
    if table.empty:
        st.warning("No pathway table is available for this analysis. There may be too few genes or the gene-set source may have been unavailable.")
        st.json(manifest["enrichment"])
    else:
        lib = st.selectbox("Gene-set library", sorted(table.library.unique()))
        table = table.loc[table.library == lib]
        term_search = st.text_input("Search pathways")
        if term_search:
            table = table.loc[table.Term.str.contains(term_search, case=False, regex=False)]
        table = table.sort_values("padj")
        if not table.empty:
            with st.container(border=True, key="card-8"):
                x = "NES" if method == "Preranked GSEA" else "fold_enrichment"
                card_head("Functional enrichment", "Top 15 terms by adjusted p; colour = adjusted p")
                fig = px.bar(table.head(15).iloc[::-1], x=x, y="Term", orientation="h", color="padj",
                             hover_data=["genes"], color_continuous_scale=SEQUENTIAL)
                fig.update_traces(marker_line_width=0)
                style_figure(fig, 560).update_yaxes(title=None)
                st.plotly_chart(fig, use_container_width=True)
        st.caption("GSEA padj is a permutation FDR; ORA padj is BH within library/direction. Table includes non-significant terms. Zero permutation estimates are resolution-limited.")
        st.dataframe(table, use_container_width=True)
        st.download_button("Download displayed pathways", table.to_csv(index=False).encode(), "oscc_pathways.csv", "text/csv")
st.markdown('<div class="oscc-footer">Research portfolio · reproducible sources and computed results · differential expression is not causality</div>', unsafe_allow_html=True)
