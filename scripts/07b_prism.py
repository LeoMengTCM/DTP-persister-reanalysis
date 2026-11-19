#!/usr/bin/env python3
"""Consensus signature vs PRISM primary screen."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
sc = pd.read_csv(f"{TBL}/S_depmap_signature_scores.csv", index_col=0)
ti = pd.read_csv(f"{CACHE}/prism/primary-screen-replicate-collapsed-treatment-info.csv")
lfc = pd.read_csv(f"{CACHE}/prism/primary-screen-replicate-collapsed-logfold-change.csv", index_col=0)
lfc = lfc[~lfc.index.str.contains("FAILED")]
common = lfc.index.intersection(sc.index); print("PRISM lines with signature score:", len(common), "| compounds:", lfc.shape[1])
s = sc.loc[common]; L = lfc.loc[common]
lin = s.lineage.where(s.lineage.map(s.lineage.value_counts()) >= 5, "Other_small")

X = pd.get_dummies(lin, drop_first=True).astype(float); X.insert(0, "const", 1.0); Xm = X.values
P = Xm @ np.linalg.pinv(Xm)
rows = []
sv = s["consensus"].values
for col in L.columns:
    y = L[col].values; ok = ~np.isnan(y)
    if ok.sum() < 200: continue
    Xo = Xm[ok]; Po = Xo @ np.linalg.pinv(Xo)
    sr = sv[ok] - Po @ sv[ok]; yr = y[ok] - Po @ y[ok]
    df = ok.sum() - np.linalg.matrix_rank(Xo) - 1
    b = (sr @ yr) / (sr @ sr); se = np.sqrt(((yr - b * sr) ** 2).sum() / df / (sr @ sr)); t = b / se
    rows.append((col, ok.sum(), t / np.sqrt(t * t + df), t))
R = pd.DataFrame(rows, columns=["column_name", "n", "partial_r", "t"]).merge(ti[["column_name", "name", "moa", "target"]], on="column_name")
from scipy import stats
R["P"] = 2 * stats.t.sf(R.t.abs(), R.n - X.shape[1] - 1); R["FDR"] = multipletests(R.P, method="fdr_bh")[1]
R["rank_most_negative"] = R.t.rank(method="min")
R.sort_values("t").to_csv(f"{TBL}/S_prism_consensus_signature_all_compounds.csv", index=False)
print("compounds tested:", len(R), "| FDR<0.05 more sensitive in signature-high:", int(((R.FDR < 0.05) & (R.t < 0)).sum()), "| less sensitive:", int(((R.FDR < 0.05) & (R.t > 0)).sum()))
print(R[R.name.isin(["erastin", "buthionine-sulfoximine", "sulfasalazine"])].round(4).to_string())
print(R.sort_values("t").head(25)[["name", "moa", "partial_r", "t", "FDR"]].round(4).to_string())
