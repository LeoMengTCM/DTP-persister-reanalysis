#!/usr/bin/env python3
"""QC and pseudobulk for GSE282074 (sc-rDSeq)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse
sc.settings.verbosity = 0

ens = pd.read_csv(f"{PROC}/ensembl2symbol.csv"); ens2sym = dict(zip(ens.ensembl, ens.symbol))
LIBS = {
    "DMSO_1": ("GSM8635332", "DMSO"), "DMSO_2": ("GSM8635333", "DMSO"),
    "Gefitinib_1": ("GSM8635334", "Gefitinib"), "Gefitinib_2": ("GSM8635335", "Gefitinib"),
    "Osimertinib_1": ("GSM8635336", "Osimertinib"),
}
ads, rows = [], []
for lib, (gsm, trt) in LIBS.items():
    fn = f"{RAW282}/{gsm}_rsem_STRAND_lRNA_count_matrix_Ensembl_counts_PC9_{lib}.tsv.gz"
    df = pd.read_csv(fn, sep="\t", index_col=0)
    n_ens = df.shape[0]
    df.index = [ens2sym.get(i, i) for i in df.index]
    df = df[~df.index.str.startswith("ENSG")]
    df = df.groupby(level=0).sum()
    a = ad.AnnData(sparse.csr_matrix(df.T.values.astype(np.float32)), obs=pd.DataFrame(index=[f"{lib}_{c}" for c in df.columns]),
                   var=pd.DataFrame(index=df.index))
    a.obs["library"] = lib; a.obs["treatment"] = trt; a.obs["group"] = "DMSO" if trt == "DMSO" else "TKI"
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    n0 = a.n_obs
    q = a.obs
    print(lib, "pre-QC", n0, "| genes/cell quantiles", q.n_genes_by_counts.quantile([.01, .5, .99]).round(0).tolist(),
          "| pct_mt quantiles", q.pct_counts_mt.quantile([.5, .95, .99]).round(1).tolist(), flush=True)

    print("   cells that would pass a 20% mito cut-off:", int((q.pct_counts_mt < 20).sum()), flush=True)
    keep = (q.n_genes_by_counts >= 200) & (q.n_genes_by_counts <= 8000)
    a = a[keep].copy()
    rows.append(dict(library=lib, GSM=gsm, treatment=trt, ensembl_ids=n_ens, symbols=a.n_vars, barcodes_deposited=n0,
                     cells_passing_QC=a.n_obs, median_genes=float(a.obs.n_genes_by_counts.median()),
                     median_UMI=float(a.obs.total_counts.median()), median_pct_mt=round(float(a.obs.pct_counts_mt.median()), 2)))
    ads.append(a)
qc = pd.DataFrame(rows); qc.to_csv(f"{TBL}/S_QC_GSE282074.csv", index=False); print(qc.to_string())
A = ad.concat(ads, join="outer", fill_value=0)
A.layers["counts"] = A.X.copy()
pb = pd.DataFrame({l: np.asarray(A[A.obs.library == l].layers["counts"].sum(0)).ravel() for l in LIBS}, index=A.var_names)
pb.round().astype(int).to_csv(f"{CACHE}/GSE282074_pseudobulk_counts.csv")
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.write_h5ad(f"{CACHE}/GSE282074_qc.h5ad")
print(A)

sets = {k: v for k, v in METABOLIC_PATHWAYS.items()}
sets["Defense_score_9genes"] = DEFENSE_GENES
json.dump({"sets": sets, "aliases": ALIASES}, open(f"{CACHE}/gene_sets.json", "w"), indent=1)
