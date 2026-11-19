#!/usr/bin/env python3
"""TF activity (CollecTRI, DoRothEA) and Hallmark scores from edgeR statistics."""
import sys, os, glob
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, decoupler as dc
from statsmodels.stats.multitest import multipletests

cdir = f"{TBL}/contrasts"
T = {}
for fn in sorted(glob.glob(f"{cdir}/DE_*.csv")):
    d = pd.read_csv(fn).set_index("gene")
    T[os.path.basename(fn)[3:-4]] = d["t"]
T = pd.DataFrame(T).T
print("contrast matrix", T.shape)

def get_net(name):
    path = f"{CACHE}/net_{name}.csv"
    if os.path.exists(path):
        return pd.read_csv(path)
    net = {"collectri": lambda: dc.op.collectri(organism="human"),
           "dorothea": lambda: dc.op.dorothea(organism="human", levels=["A", "B", "C"]),
           "hallmark": lambda: dc.op.hallmark(organism="human")}[name]()
    net.to_csv(path, index=False); return net

def run(net, tag, tmin=5):
    out = []
    for c in T.index:
        x = T.loc[[c]].dropna(axis=1)
        es, pv = dc.mt.ulm(data=x, net=net, tmin=tmin)

        df = pd.DataFrame({"source": es.columns, "score": es.iloc[0].values, "FDR": pv.iloc[0].values})
        nt = net[net.target.isin(x.columns)].groupby("source").size()
        df["n_targets"] = df.source.map(nt).values
        df["contrast"] = c
        out.append(df)
    out = pd.concat(out); out.to_csv(f"{TBL}/Table_{tag}_all_contrasts.csv", index=False)
    return out

collectri = get_net("collectri"); print("CollecTRI:", len(collectri), "interactions;", collectri.source.nunique(), "TFs")
tf = run(collectri, "TF_activity_CollecTRI")
doro = get_net("dorothea"); print("DoRothEA A-C:", len(doro), "interactions;", doro.source.nunique(), "TFs")
tfd = run(doro, "TF_activity_DoRothEA")
hm = get_net("hallmark"); hm = hm.drop_duplicates(["source", "target"]); print("Hallmark:", hm.source.nunique(), "sets")
if "weight" not in hm.columns: hm["weight"] = 1.0
hall = run(hm, "Hallmark_ULM")

pd.set_option("display.width", 330); pd.set_option("display.max_columns", 40)
short = lambda c: (c.replace("A_PC9_erlotinib_", "A:").replace("B_PC9_gef_osi_9d_vs_DMSO", "B:9d").replace("C_", "")
                   .replace("_vs_DMSO", "").replace("pooled4lines", "POOL").replace("_vs_naive", "").replace("longwash", "wash"))
def show(df, sources, title):
    d = df[df.source.isin(sources)].copy()
    d["cell"] = d.apply(lambda r: f"{r.score:+.1f}" + ("***" if r.FDR < 0.001 else "**" if r.FDR < 0.01 else "*" if r.FDR < 0.05 else ""), axis=1)
    d["c"] = d.contrast.map(short)
    cols = [c for c in d.c.unique() if "vs_acute" not in c and "vs_D3" not in c]
    print("\n" + title); print(d.pivot(index="source", columns="c", values="cell").reindex(sources)[cols].to_string())

named = ["SREBF1", "SREBF2", "NFE2L2", "ATF3", "FOXO3", "ID4", "DLX5", "TFE3", "MITF", "STAT3", "PPARG", "PPARD", "CREBBP",
         "JUN", "SP1", "TP53", "NFKB1", "NFKB", "AP1", "ATF4", "HIF1A", "MYC", "E2F1", "E2F4", "TEAD1", "YAP1", "SMAD3", "KDM5A", "BACH1", "CLOCK", "AR", "CREB1"]
show(tf, named, "CollecTRI ULM score (signed) per contrast")
show(tfd, [s for s in named if s in set(tfd.source)], "DoRothEA A-C ULM score per contrast")

dtp_contrasts = ["A_PC9_erlotinib_late_vs_naive", "B_PC9_gef_osi_9d_vs_DMSO", "C_PC9_DTP_vs_DMSO", "C_HCC827_DTP_vs_DMSO",
                 "C_HCC2935_DTP_vs_DMSO", "C_H1975_DTP_vs_DMSO"]
def consistency(df, tag):
    d = df[df.contrast.isin(dtp_contrasts)]
    sc_ = d.pivot(index="source", columns="contrast", values="score")[dtp_contrasts]
    fd = d.pivot(index="source", columns="contrast", values="FDR")[dtp_contrasts]
    res = pd.DataFrame({"n_contrasts": sc_.notna().sum(1), "n_up_FDR05": ((sc_ > 0) & (fd < 0.05)).sum(1),
                        "n_down_FDR05": ((sc_ < 0) & (fd < 0.05)).sum(1), "n_up": (sc_ > 0).sum(1), "n_down": (sc_ < 0).sum(1),
                        "mean_score": sc_.mean(1)})
    res = pd.concat([res, sc_.add_prefix("score_"), fd.add_prefix("FDR_")], axis=1).sort_values("mean_score")
    res.to_csv(f"{TBL}/Table_{tag}_consistency.csv"); return res
ctf = consistency(tf, "TF_activity_CollecTRI")
print("\nTFs consistently repressed (FDR<0.05 in >=5 of 6 DTP contrasts):", ctf[ctf.n_down_FDR05 >= 5].index.tolist())
print("TFs consistently activated (FDR<0.05 in >=5 of 6 DTP contrasts):", ctf[ctf.n_up_FDR05 >= 5].index.tolist())
print("significant in PC9 erlotinib late-vs-naive (FDR<0.05):", int((tf[tf.contrast == dtp_contrasts[0]].FDR < 0.05).sum()), "of", int((tf.contrast == dtp_contrasts[0]).sum()))
print(ctf.head(15)[["n_down_FDR05", "n_up_FDR05", "mean_score"]].round(2).to_string()); print(ctf.tail(15)[["n_down_FDR05", "n_up_FDR05", "mean_score"]].round(2).to_string())
ch = consistency(hall, "Hallmark_ULM")
show(hall, ch.index.tolist(), "Hallmark ULM score per contrast (sorted by mean across DTP contrasts)")

targets = {"SREBP_targets": ["LDLR", "INSIG1", "HMGCS1", "HMGCR", "MVD", "MVK", "FDPS", "FDFT1", "SQLE", "LSS", "DHCR7", "DHCR24",
                             "IDI1", "SC5D", "FASN", "ACACA", "ACLY", "SCD", "ELOVL6", "FADS1", "FADS2", "ACSS2", "SREBF1", "SREBF2"],
           "NRF2_targets": ["NQO1", "HMOX1", "GCLC", "GCLM", "TXNRD1", "SRXN1", "SLC7A11", "G6PD", "PGD", "TALDO1", "ME1",
                            "AKR1C1", "AKR1C2", "AKR1C3", "AKR1B10", "GPX2", "ABCC2", "FTH1", "FTL", "OSGIN1", "NFE2L2", "KEAP1"],
           "GPX4_axis": ["GPX4", "AIFM2", "DHODH", "ACSL4", "LPCAT3", "SLC3A2", "TFRC", "ALDH1A1", "ALDH1A3", "ALDH3A1", "ALDH2"]}
rows = []
for fn in sorted(glob.glob(f"{cdir}/DE_*.csv")):
    c = os.path.basename(fn)[3:-4]; d = pd.read_csv(fn).set_index("gene")
    for k, gl in targets.items():
        for g in gl:
            if g in d.index:
                rows.append(dict(group=k, gene=g, contrast=c, logFC=d.loc[g, "logFC"], PValue=d.loc[g, "PValue"], FDR=d.loc[g, "FDR"]))
gt = pd.DataFrame(rows); gt.to_csv(f"{TBL}/Table_target_gene_logFC.csv", index=False)
for k in targets:
    d = gt[gt.group == k].copy(); d["c"] = d.contrast.map(short)
    d["cell"] = d.apply(lambda r: f"{r.logFC:+.1f}" + ("*" if r.FDR < 0.05 else ""), axis=1)
    cols = [c for c in d.c.unique() if "vs_acute" not in c and "vs_D3" not in c and "wash" not in c and "POOL" not in c]
    print("\n" + k); print(d.pivot(index="gene", columns="c", values="cell").reindex(targets[k])[cols].to_string())
    dd = d[d.contrast.isin(dtp_contrasts)]
    print("  genes down (FDR<0.05) in >=4/6 DTP contrasts:", dd[(dd.FDR < 0.05) & (dd.logFC < 0)].groupby("gene").size().loc[lambda s: s >= 4].index.tolist())
    print("  genes up (FDR<0.05) in >=4/6 DTP contrasts:", dd[(dd.FDR < 0.05) & (dd.logFC > 0)].groupby("gene").size().loc[lambda s: s >= 4].index.tolist())
print("done")
