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
st.markdown("""<style>
.block-container {padding-top:2rem; max-width:1450px;}
h1 {letter-spacing:-.035em;} [data-testid="stMetric"] {background:#f1f5f8;border-radius:8px;padding:16px;}
</style>""", unsafe_allow_html=True)
dataset = st.sidebar.selectbox("Cohort", ["GSE20116", "GSE184616"], key="dataset")
cohort = get_dataset(dataset)
output = bundle_directory(ROOT, dataset)
st.caption(f"RESEARCH EXPLORER  /  {dataset}  /  {cohort.role.upper()}  /  PAIRED BULK RNA-SEQ")
st.title("Transcriptomic Profiling of Oral Squamous Cell Carcinoma")
st.write("Explore tumor–normal expression differences while preserving each patient's matched comparison.")


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
                     color_discrete_map={"up": "#B74B44", "down": "#327D9D", "not_significant": "#B8BEC6"},
                     labels={"log2FoldChange": "log2 fold change · tumor / normal", "neglog10_q": "−log10(BH-adjusted p)"})
    fig.add_hline(y=-np.log10(alpha), line_dash="dash", line_color="grey")
    for x in [-lfc, lfc]:
        fig.add_vline(x=x, line_dash="dash", line_color="grey")
    fig.update_traces(marker={"size": 5, "opacity": .7})
    fig.update_layout(height=530, template="plotly_white")
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
                  labels={"log_expression": "log2(normalized count + 1)", "patient_id": "Patient"}, title=gene)
    fig.update_layout(template="plotly_white")
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
        fig = px.imshow(matrix, color_continuous_scale="RdBu_r", zmin=-2, zmax=2, aspect="auto", labels={"color": "Row z-score"})
        fig.update_layout(height=650)
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
                         hover_data=["genes"], color_continuous_scale="Viridis_r")
            fig.update_layout(height=600, template="plotly_white")
            st.plotly_chart(fig, use_container_width=True)
        st.dataframe(table, use_container_width=True)
        st.download_button("Download displayed pathways", table.to_csv(index=False).encode(), "oscc_pathways.csv", "text/csv")
st.caption("Research portfolio · reproducible sources and computed results · differential expression is not causality")
