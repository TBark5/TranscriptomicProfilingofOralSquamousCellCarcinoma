"""Local GMT enrichment: tested background ORA and signed Wald preranking."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
import requests
from scipy.stats import hypergeom
from .acquisition import download
from .differential import bh_adjust, classify

LIBRARIES = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "Reactome_2022"]
HALLMARK_SHA = "4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81"
ORA_COLUMNS = ["library", "direction", "Term", "overlap", "set_size", "query_size", "background_size", "fold_enrichment", "pvalue", "padj", "genes"]
GSEA_COLUMNS = ["library", "Term", "ES", "NES", "pvalue", "padj", "fwer", "genes"]


def read_gmt(path):
    sets = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.rstrip().split("\t")
        if len(fields) < 3:
            raise ValueError(f"Invalid GMT: {path}")
        if fields[0] in sets:
            raise ValueError("Duplicate GMT term")
        sets[fields[0]] = sorted(set(x.strip() for x in fields[2:] if x.strip()))
    if not sets:
        raise ValueError("Empty gene-set library")
    return sets


def prepare_ranking(results):
    """Finite tested Wald statistics; stable alphabetical tie-breaking, no jitter."""
    r = results.loc[results.pvalue.notna() & results.model_converged, ["stat"]].copy()
    r = r.loc[np.isfinite(r.stat)]
    r.index = r.index.astype(str)
    if r.index.has_duplicates:
        raise ValueError("Duplicate symbols: resolve representative features before enrichment")
    r.index.name = "gene"
    return r.reset_index().sort_values(["stat", "gene"], ascending=[False, True], kind="stable")


def overrepresentation(query, background, gene_sets, min_size=15, max_size=500):
    """Hypergeometric ORA includes zero-overlap sets in the BH family."""
    background, query = set(background), set(query)
    if not query <= background:
        raise ValueError("ORA query must be a subset of the tested background")
    rows = []
    for term, members in gene_sets.items():
        members = set(members) & background
        if not min_size <= len(members) <= max_size:
            continue
        overlap = query & members
        k, n, N, M = len(overlap), len(members), len(query), len(background)
        rows.append({"Term": term, "overlap": k, "set_size": n, "query_size": N,
                     "background_size": M, "fold_enrichment": (k / N) / (n / M) if N else 0.,
                     "pvalue": float(hypergeom.sf(k - 1, M, n, N)), "genes": ";".join(sorted(overlap))})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["padj"] = bh_adjust(result.pvalue)
    return result


def run_enrichment(results, root, output, libraries, alpha=.05, lfc=1., seed=42, permutations=1000, offline=False):
    output, raw = Path(output), Path(root) / "data/raw"
    ranking = prepare_ranking(results)
    ranking.to_csv(output / "ranking.csv", index=False, lineterminator='\n')
    background = set(ranking.gene)
    calls = classify(results, alpha, lfc)
    ora_tables, gsea_tables, status = [], [], []
    curve = None
    for library in libraries:
        if library not in LIBRARIES:
            raise ValueError(f"Unsupported human library {library}")
        path = raw / f"{library}.gmt"
        url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library}"
        try:
            provenance = download(url, path, HALLMARK_SHA if library == LIBRARIES[0] else None, offline)
            gene_sets = read_gmt(path)
        except (OSError, ValueError, requests.RequestException) as exc:
            status.append({"library": library, "status": "unavailable", "reason": str(exc)})
            continue
        eligible = {t: sorted(set(g) & background) for t, g in gene_sets.items()
                    if 15 <= len(set(g) & background) <= 500}
        members = set().union(*(set(g) for g in gene_sets.values()))
        pd.DataFrame({"gene": sorted(background), "in_library": [g in members for g in sorted(background)]}).to_csv(output / f"mapping_{library}.csv", index=False, lineterminator='\n')
        entry = {"library": library, "source": provenance, "status": "complete", "tested_background": len(background),
                 "symbols_in_library": len(background & members), "eligible_sets": len(eligible),
                 "ranking_ties": int(ranking.stat.duplicated().sum()), "ora": {}}
        for direction in ("up", "down"):
            query = set(results.index[calls == direction]) & background
            if len(query) < 5:
                entry["ora"][direction] = f"skipped: {len(query)} genes; require at least 5"
                continue
            ora = overrepresentation(query, background, gene_sets)
            if not ora.empty:
                ora["library"], ora["direction"] = library, direction
                ora_tables.append(ora)
            entry["ora"][direction] = f"tested {len(ora)} terms with {len(query)} query genes"
        if eligible and len(ranking) > 15:
            gsea = gp.prerank(rnk=ranking, gene_sets=eligible, threads=1, min_size=15, max_size=500,
                              permutation_num=permutations, seed=seed, outdir=None, no_plot=True, verbose=False)
            table = gsea.res2d.rename(columns={"NOM p-val": "pvalue", "FDR q-val": "padj", "FWER p-val": "fwer", "Lead_genes": "genes"})
            table["library"] = library
            gsea_tables.append(table[GSEA_COLUMNS])
            if curve is None and not table.empty:
                term = table.sort_values(["padj", "pvalue", "Term"]).iloc[0].Term
                data = gsea.results[term]
                curve = {"term": term, "library": library, "running_es": list(data["RES"]),
                         "hits": list(data["hits"]), "nes": float(data["nes"]), "fdr": float(data["fdr"])}
        else:
            entry["gsea"] = "skipped: no eligible gene sets"
        status.append(entry)
    ora = pd.concat(ora_tables, ignore_index=True) if ora_tables else pd.DataFrame(columns=ORA_COLUMNS)
    gsea = pd.concat(gsea_tables, ignore_index=True) if gsea_tables else pd.DataFrame(columns=GSEA_COLUMNS)
    ora.reindex(columns=ORA_COLUMNS).to_csv(output / "ora.csv", index=False, lineterminator='\n')
    gsea.to_csv(output / "gsea.csv", index=False, lineterminator='\n')
    (output / "enrichment_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    if curve:
        (output / "enrichment_curve.json").write_text(json.dumps(curve), encoding="utf-8")
    else:
        (output / "enrichment_curve.json").unlink(missing_ok=True)
    return ora, gsea, status, curve
