#!/usr/bin/env python3
"""QC of the 11 Drop-seq libraries (GSE149383) and per-library pseudobulk counts."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse

SAMPLES = [

    ("D0_a",  "GSM3972657_D0.dge.txt.gz",                    "D0",  0,  "timecourse", True),
    ("D1_a",  "GSM3972658_D1.dge.txt.gz",                    "D1",  1,  "timecourse", True),
    ("D2_a",  "GSM3972659_D2.dge.txt.gz",                    "D2",  2,  "timecourse", True),
    ("D4_a",  "GSM3972660_D4.dge.txt.gz",                    "D4",  4,  "timecourse", True),
    ("D9_a",  "GSM3972661_D9.dge.txt.gz",                    "D9",  9,  "timecourse", True),
    ("D11_a", "GSM3972662_D11.dge.txt.gz",                   "D11", 11, "timecourse", True),
    ("D0_b",  "GSM4494347_D0C_expression_matrix.txt.gz",     "D0",  0,  "replicate",  False),
    ("D11_b", "GSM4494348_D11C_expression_matrix.txt.gz",    "D11", 11, "replicate",  False),
    ("D11_c", "GSM4494349_D11_ERL_expression_matrix.txt.gz", "D11", 11, "replicate",  False),
    ("D3_a",  "GSM4494350_PC9D3_ERL1_expression_matrix.txt.gz", "D3", 3, "replicate", False),
    ("D3_b",  "GSM4494351_PC9D3_ERL2_expression_matrix.txt.gz", "D3", 3, "replicate", False),
]

ads, qc_rows = [], []
for sid, fn, tp, day, series, prim in SAMPLES:
    df = pd.read_csv(f"{RAW149}/{fn}", sep="\t", index_col=0)
    df = df[~df.index.duplicated()]
    a = ad.AnnData(X=sparse.csr_matrix(df.T.values.astype(np.float32)),
                   obs=pd.DataFrame(index=[f"{sid}_{bc}" for bc in df.columns]),
                   var=pd.DataFrame(index=df.index.astype(str)))
    a.obs["sample"] = sid; a.obs["timepoint"] = tp; a.obs["day"] = day
    a.obs["series"] = series; a.obs["in_primary"] = prim
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    n0 = a.n_obs
    keep = (a.obs.n_genes_by_counts >= 200) & (a.obs.n_genes_by_counts <= 8000) & (a.obs.pct_counts_mt < 20)
    n_lowgene = int((a.obs.n_genes_by_counts < 200).sum())
    a = a[keep].copy()
    qc_rows.append(dict(sample=sid, file=fn, timepoint=tp, series=series, barcodes_deposited=n0,
                        barcodes_lt200_genes=n_lowgene, cells_passing_QC=a.n_obs,
                        pct_retained=round(100 * a.n_obs / n0, 1),
                        median_genes=float(a.obs.n_genes_by_counts.median()),
                        median_UMI=float(a.obs.total_counts.median()),
                        median_pct_mt=float(a.obs.pct_counts_mt.median())))
    print(qc_rows[-1], flush=True)
    ads.append(a)

adata = ad.concat(ads, join="outer", fill_value=0)
adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()
qc = pd.DataFrame(qc_rows)
qc.to_csv(f"{TBL}/S_QC_per_sample.csv", index=False)

pb = pd.DataFrame({s: np.asarray(adata[adata.obs["sample"] == s].layers["counts"].sum(axis=0)).ravel()
                   for s in qc["sample"]}, index=adata.var_names)
pb.to_csv(f"{CACHE}/pseudobulk_counts.csv")
meta = qc[["sample", "timepoint", "series", "cells_passing_QC"]].copy()
meta["group"] = meta.timepoint.map(lambda t: "naive" if t == "D0" else ("DTP" if t in ("D9", "D11") else "other"))
meta.to_csv(f"{CACHE}/pseudobulk_meta.csv", index=False)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad(f"{CACHE}/pc9_extended.h5ad")
print(adata)
print(qc.to_string())
