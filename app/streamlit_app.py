"""Run from the repository: streamlit run app/streamlit_app.py."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
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

st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,400;9..144,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');
html, body, [class*="css"], .stApp, button, input, textarea, select {font-family: 'IBM Plex Sans', system-ui, sans-serif;}
.stApp {background: #F6F4EF;}
.block-container {padding-top: 2.25rem; padding-bottom: 3rem; max-width: 1400px;}
h1, h2, h3, .stApp h1, .stApp h2, .stApp h3 {font-family: 'Fraunces', Georgia, serif !important; letter-spacing: -.015em; color: #16232E;}
h2, h3 {font-weight: 600 !important;}
code, pre, [data-testid="stMetricValue"] {font-family: 'IBM Plex Mono', monospace;}

/* Hero */
.oscc-chips {display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 14px;}
.oscc-chip {padding: 5px 12px; border-radius: 999px; background: #EFEAE0; color: #5A4E3A; font-size: 12px; font-weight: 600; letter-spacing: .06em; text-transform: uppercase;}
.oscc-chip.accent {background: #E6EEF1; color: #1D5C7A;}
.oscc-title, .stApp h1.oscc-title {font-family: 'Fraunces', Georgia, serif !important; font-weight: 400 !important; padding: 0 !important; font-size: clamp(2rem, 3.6vw, 2.9rem); line-height: 1.08; letter-spacing: -.02em; margin: 0 0 1rem; color: #16232E;}
.oscc-lede {font-size: 1.06rem; line-height: 1.55; color: #4A5763; max-width: 760px; margin: 0 0 1.5rem;}
.oscc-footer {margin-top: 2.5rem; padding-top: 1rem; border-top: 1px solid #E3DED4; font-size: 12px; color: #6B7780;}

/* Metric cards */
[data-testid="stMetric"] {background: #FFFFFF; border: 1px solid #E3DED4; border-radius: 16px; padding: 18px 22px; box-shadow: 0 1px 2px rgba(20,33,43,.04);}
[data-testid="stMetricLabel"] p {font-size: 13px; color: #5A6772;}
[data-testid="stMetricValue"] {font-family: 'Fraunces', Georgia, serif; font-size: 2.3rem; color: #16232E;}
[data-testid="stColumn"]:nth-of-type(2) [data-testid="stMetricValue"] {color: #9E3F1F;}
[data-testid="stColumn"]:nth-of-type(3) [data-testid="stMetricValue"] {color: #1D5C7A;}

/* Tabs as a segmented control */
.stTabs [data-baseweb="tab-list"] {gap: 4px; padding: 4px; background: #ECE8E0; border-radius: 12px; width: fit-content; max-width: 100%; overflow-x: auto;}
.stTabs [data-baseweb="tab"] {height: 44px; padding: 0 18px; border-radius: 9px; color: #4A5763; background: transparent;}
.stTabs [aria-selected="true"] {background: #FFFFFF; color: #16232E; box-shadow: 0 1px 2px rgba(20,33,43,.08);}
.stTabs [aria-selected="true"] p {font-weight: 600;}
.stTabs [data-baseweb="tab-highlight"], .stTabs [data-baseweb="tab-border"] {display: none;}
.stTabs [data-baseweb="tab-panel"] {padding-top: 1.5rem;}

/* Cards around charts, tables and images */
[data-testid="stPlotlyChart"] {background: #FFFFFF; border: 1px solid #E3DED4; border-radius: 16px; overflow: hidden;}
[data-testid="stDataFrame"] {border: 1px solid #E3DED4; border-radius: 12px; overflow: hidden;}
[data-testid="stImage"] img {border: 1px solid #E3DED4; border-radius: 16px; background: #FFFFFF;}
[data-testid="stAlert"] {border-radius: 14px;}
[data-testid="stExpander"] details {background: #FFFFFF; border: 1px solid #E3DED4; border-radius: 14px;}
.stDownloadButton button, .stButton button {border-radius: 10px; min-height: 44px; border: 1px solid #C9D6DC; font-weight: 600;}

/* Sidebar */
[data-testid="stSidebar"] {background: #14212B;}
[data-testid="stSidebar"] * {color: #E9E5DC;}
[data-testid="stSidebar"] h2 {font-size: 12px !important; font-family: 'IBM Plex Sans', sans-serif !important; letter-spacing: .1em; text-transform: uppercase; color: #A9B4BC !important; margin-top: 1rem;}
[data-testid="stSidebar"] [data-testid="stCaptionContainer"] p {color: #A9B4BC; line-height: 1.55;}
[data-testid="stSidebar"] [data-baseweb="select"] > div {background: #1F303D; border-color: #2E4250;}
[data-testid="stSidebar"] [data-testid="stSliderThumbValue"], [data-testid="stSidebar"] [data-testid="stSliderTickBarMin"], [data-testid="stSidebar"] [data-testid="stSliderTickBarMax"] {color: #E9A27F; font-family: 'IBM Plex Mono', monospace;}
[data-testid="stSidebar"] [role="slider"] {background: #E9A27F; border-color: #E9A27F;}
.oscc-brand {display: flex; align-items: center; gap: 12px; margin: .25rem 0 1.5rem;}
.oscc-brand b {font-family: 'Fraunces', Georgia, serif; font-size: 20px; font-weight: 600; display: block;}
.oscc-brand span {font-size: 11px; letter-spacing: .08em; text-transform: uppercase; color: #A9B4BC;}
</style>""", unsafe_allow_html=True)


def style_figure(fig, height=None):
    """Apply the dashboard's typography, grid and palette to a Plotly figure."""
    fig.update_layout(template="plotly_white", paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF",
                      font={"family": "IBM Plex Sans, system-ui, sans-serif", "color": INK, "size": 13},
                      margin={"l": 16, "r": 16, "t": 56, "b": 16},
                      legend={"title": None, "orientation": "h", "y": 1.02, "yanchor": "bottom", "x": 1, "xanchor": "right"},
                      hoverlabel={"bgcolor": "#FFFFFF", "bordercolor": LINE, "font_family": "IBM Plex Sans, sans-serif"})
    fig.update_xaxes(gridcolor="#EFECE6", zerolinecolor="#B9BEC2", linecolor=LINE)
    fig.update_yaxes(gridcolor="#EFECE6", zerolinecolor="#B9BEC2", linecolor=LINE)
    if fig.layout.title.text:
        fig.update_layout(title_font={"family": "Fraunces, Georgia, serif", "size": 20})
    if height:
        fig.update_layout(height=height)
    return fig


st.sidebar.markdown("""<div class="oscc-brand"><svg width="36" height="36" viewBox="0 0 36 36" fill="none" stroke="#E9A27F" stroke-width="2" stroke-linecap="round" aria-hidden="true"><path d="M12 4c0 9 12 9 12 18s-12 9-12 10"/><path d="M24 4c0 9-12 9-12 18s12 9 12 10"/><path d="M14 10h8M13.5 26h9M16 18h4"/></svg><div><b>OSCC Explorer</b><span>Paired bulk RNA-seq</span></div></div>""", unsafe_allow_html=True)
dataset = st.sidebar.selectbox("Cohort", ["GSE20116", "GSE184616"], key="dataset")
cohort = get_dataset(dataset)
output = bundle_directory(ROOT, dataset)
st.markdown(f"""<div class="oscc-chips"><span class="oscc-chip accent">{dataset}</span><span class="oscc-chip">{cohort.role} cohort</span><span class="oscc-chip">Paired bulk RNA-seq</span></div>
<h1 class="oscc-title">Transcriptomic profiling of oral squamous cell carcinoma</h1>
<p class="oscc-lede">Explore tumor–normal expression differences while preserving each patient's matched comparison.</p>""", unsafe_allow_html=True)


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
st.warning(f"{cohort.description} {manifest['diagnostics']['residual_df']} residual degrees of freedom. Exploratory associations; no clinical claims.")
if dataset == "GSE184616":
    summary = manifest["summary"]
    st.info(f"Independent validation: {summary['replicated_candidates']:,} of {summary['testable_candidates']:,} testable discovery candidates replicated; {summary['replicated_pathways']} discovery pathways replicated. See VALIDATION_REPORT.md for criteria and limitations.")
st.sidebar.header("Gene thresholds")
alpha = st.sidebar.slider("BH-adjusted p-value cutoff", .001, .20, float(max(.001, min(.20, manifest["config"]["alpha"]))), .001, format="%.3f", key="alpha")
lfc = st.sidebar.slider("Absolute log2 fold-change cutoff", 0., 5., float(min(manifest["config"]["lfc"], 5)), .1, key="lfc")
st.sidebar.caption("Strict inequalities. LFC is tumor / normal. These controls filter the fitted results; they do not refit the model or recalculate pathways.")
calls = classify(results, alpha, lfc)
metrics = st.columns(4)
for column, label, value in zip(metrics, ["Gene representatives" if dataset == "GSE20116" else "Genes", "Upregulated", "Downregulated", "Samples / pairs"],
                               [f"{len(results):,}", f"{sum(calls == 'up'):,}", f"{sum(calls == 'down'):,}", f"{len(metadata)} / {metadata.patient_id.nunique()}"]):
    column.metric(label, value)

overview, genes_tab, expression_tab, pathways_tab = st.tabs(["Overview & provenance", "Differential expression", "Gene & heatmap explorer", "Pathways"])
with overview:
    st.subheader(f"{dataset}: {cohort.role} cohort")
    st.write(cohort.description)
    st.markdown(f"[GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset})")
    if dataset == "GSE20116":
        st.write("Original publication raw counts; one RefSeq representative per historical symbol, selected independently of differential expression.")
    else:
        st.write("Deposited unnormalized gene counts; sample titles and patient characteristics verified against GEO. Exact-symbol comparisons use frozen original discovery candidates.")
    st.dataframe(metadata, use_container_width=True)
    st.write("Model: PyDESeq2 negative-binomial regression, `~ patient_id + condition`, tumor versus normal Wald contrast, median-of-ratios normalization, convergence-qualified BH correction.")
    st.caption(f"{manifest['diagnostics']['failed_convergence']} features failed an optimizer convergence check and are excluded from significance calls. See the report for assumptions and the complete diagnostics.")
    if (output / "figures/pca.png").exists():
        st.image(str(output / "figures/pca.png"), caption="PCA on up to 2,000 most variable log normalized genes; lines preserve patient pairing.", width=750)
    with st.expander("Source hashes, configuration and model diagnostics"):
        st.json(manifest)

with genes_tab:
    st.subheader("Tumor-relative differential expression")
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
    legend_names = {"up": "Up in tumor", "down": "Down in tumor", "not_significant": "Not significant"}
    fig.for_each_trace(lambda trace: trace.update(name=legend_names.get(trace.name, trace.name)))
    style_figure(fig, 540)
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Missing adjusted p-values are omitted; zero values are clipped to 1e−300 for display only.")
    search = st.text_input("Search source gene symbols", placeholder="e.g. MMP")
    show_all = st.checkbox("Include genes that do not pass the thresholds", value=False)
    view = results.copy() if show_all else filtered_results(results, alpha, lfc)
    view["direction"] = classify(view, alpha, lfc)
    if search:
        view = view.loc[view.index.str.contains(search, case=False, regex=False)]
    st.dataframe(view.sort_values("padj"), use_container_width=True)
    st.download_button("Download displayed gene results", view.to_csv().encode(), "oscc_genes.csv", "text/csv")

with expression_tab:
    st.subheader("Patient-matched expression")
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
    st.subheader("Top-gene explorer")
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
    st.subheader("Functional enrichment")
    st.info(f"Pathways were computed using the saved pipeline thresholds (q < {manifest['config']['alpha']}, |LFC| > {manifest['config']['lfc']}); sidebar sliders affect gene views only. Positive NES denotes the tumor-upregulated end of the ranking, not pathway activation.")
    method = st.radio("Analysis", ["Preranked GSEA", "Overrepresentation"], horizontal=True)
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
        st.caption("GSEA padj is a permutation FDR; ORA padj is BH within library/direction. Table includes non-significant terms. Zero permutation estimates are resolution-limited.")
        if not table.empty:
            x = "NES" if method == "Preranked GSEA" else "fold_enrichment"
            fig = px.bar(table.head(15).iloc[::-1], x=x, y="Term", orientation="h", color="padj",
                         hover_data=["genes"], color_continuous_scale=SEQUENTIAL)
            fig.update_traces(marker_line_width=0)
            style_figure(fig, 600).update_yaxes(title=None)
            st.plotly_chart(fig, use_container_width=True)
        st.dataframe(table, use_container_width=True)
        st.download_button("Download displayed pathways", table.to_csv(index=False).encode(), "oscc_pathways.csv", "text/csv")
st.markdown('<div class="oscc-footer">Research portfolio · reproducible sources and computed results · differential expression is not causality</div>', unsafe_allow_html=True)
