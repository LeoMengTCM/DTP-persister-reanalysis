#!/usr/bin/env python3
"""Consensus signature vs sensitivity to ferroptosis inducers."""
import sys, os, json, re
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd
from scipy import stats
sc_ = pd.read_csv(f"{TBL}/S_depmap_signature_scores.csv", index_col=0)
model = pd.read_csv(f"{DEPMAP}/Model.csv").set_index("ModelID")
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())
name2id = {}
for mid, r in model.iterrows():
    for n in (r.get("StrippedCellLineName"), r.get("CellLineName")):
        if isinstance(n, str): name2id.setdefault(norm(n), mid)
def partial(y, s, covs):
    X = np.column_stack([np.ones(len(y))] + covs)
    P = X @ np.linalg.pinv(X); yr = y - P @ y; sr = s - P @ s
    df = len(y) - np.linalg.matrix_rank(X) - 1
    r = (yr @ sr) / np.sqrt((yr @ yr) * (sr @ sr)); t = r * np.sqrt(df / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), df)
rows = []
for name, dataset in [("ML210", "CTRPv2"), ("ML162", "CTRPv2"), ("RSL3", "CTRPv2"), ("Erastin", "CTRPv2"), ("Piperlongumine", "CTRPv2"),
                      ("Erastin", "gCSI"), ("Piperlongumine", "GDSC1")]:
    ex = json.load(open(f"{CACHE}/pharmacodb/experiments_{name}.json"))
    d = pd.DataFrame([(e["cell_line"]["name"], e["profile"]["AAC"]) for e in ex if e["dataset"]["name"] == dataset and e["profile"] and e["profile"]["AAC"] is not None],
                     columns=["cell", "AAC"])
    d["ModelID"] = d.cell.map(lambda c: name2id.get(norm(c)))
    n_all = d.cell.nunique()
    d = d.dropna(subset=["ModelID"]).groupby("ModelID").AAC.mean()
    idx = d.index.intersection(sc_.index)
    y = d.loc[idx].values.astype(float); lin = sc_.loc[idx, "lineage"]
    vc = lin.value_counts(); lin = lin.where(lin.map(vc) >= 5, "Other_small")
    D = pd.get_dummies(lin, drop_first=True).astype(float); covs = [D[c].values for c in D.columns]
    for sig in ("consensus", "pc9_erl", "orig6", "hallmark_EMT"):
        s = sc_.loc[idx, sig].values.astype(float)
        r0, p0 = stats.pearsonr(s, y); r1, p1 = partial(y, s, covs)
        r2, p2 = partial(y, s, covs + [sc_.loc[idx, "hallmark_EMT"].values.astype(float)]) if sig != "hallmark_EMT" else (np.nan, np.nan)
        ns = sc_.loc[idx, "NSCLC"].values.astype(bool)
        rn, pn = stats.pearsonr(s[ns], y[ns]) if ns.sum() > 10 else (np.nan, np.nan)
        rows.append(dict(compound=name, dataset=dataset, signature=sig, cell_lines_in_dataset=n_all, n_lines=len(idx), r_unadjusted=r0, P_unadjusted=p0,
                         partial_r_lineage=r1, P_lineage=p1, partial_r_lineage_EMT=r2, P_lineage_EMT=p2, n_NSCLC=int(ns.sum()), r_NSCLC=rn, P_NSCLC=pn))
R = pd.DataFrame(rows); R.to_csv(f"{TBL}/S_ctrp_signature_drug_sensitivity.csv", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 30)
print(R.round(4).to_string())
