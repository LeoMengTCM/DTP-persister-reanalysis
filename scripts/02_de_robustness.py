#!/usr/bin/env python3
"""Single-cell DE with cell-number and depth matching, compared with pseudobulk."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
set_pub_style()
sc.settings.verbosity = 0

ad_ = sc.read_h5ad(f"{PROC}/PC9_DTP_timecourse_final.h5ad")
ad_.obs["timepoint"] = pd.Categorical(ad_.obs["timepoint"], categories=TP_ORDER, ordered=True)

s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','MLF1IP',
           'HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51',
           'RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8']
g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B','MKI67','TMPO','CENPF',
             'TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B',
             'HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2',
             'KIF23','HMMR','AURKA','PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']
s_ok = [g for g in s_genes if g in ad_.raw.var_names]; g2m_ok = [g for g in g2m_genes if g in ad_.raw.var_names]
sc.tl.score_genes_cell_cycle(ad_, s_genes=s_ok, g2m_genes=g2m_ok, use_raw=True)
print(f"cell-cycle genes used: S {len(s_ok)}/{len(s_genes)}, G2M {len(g2m_ok)}/{len(g2m_genes)}")

ct = pd.crosstab(ad_.obs["leiden"], ad_.obs["timepoint"])
ct["n_cells"] = ct.sum(axis=1)
frac = ct[TP_ORDER].div(ct["n_cells"], axis=0).round(3).add_prefix("frac_")
phase = pd.crosstab(ad_.obs["leiden"], ad_.obs["phase"], normalize="index").round(3).add_prefix("phase_")
qc = ad_.obs.groupby("leiden", observed=True)[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].median().round(1).add_prefix("median_")
comp = pd.concat([ct, frac, phase, qc], axis=1)
comp.to_csv(f"{TBL}/S_cluster_composition.csv")
print(comp.to_string())

tp_share = pd.crosstab(ad_.obs["timepoint"], ad_.obs["leiden"], normalize="index").round(3)
tp_share.to_csv(f"{TBL}/S_timepoint_cluster_share.csv")
print(tp_share.to_string())
cc_tp = pd.crosstab(ad_.obs["timepoint"], ad_.obs["phase"], normalize="index").round(3)
cc_tp["n"] = ad_.obs["timepoint"].value_counts().reindex(TP_ORDER).values
cc_tp.to_csv(f"{TBL}/S_cellcycle_by_timepoint.csv")
print(cc_tp.to_string())

marker_groups = {
    "Proliferation": ["MKI67", "TOP2A", "PCNA"],
    "Stress response": ["ATF3", "DDIT3", "JUN"],
    "Persister / EMT": ["AXL", "IGFBP3", "VIM", "FN1", "SERPINE1", "CYP1B1"],
    "Iron / redox": ["FTL", "FTH1", "ALDH1A3", "GPX4"],
    "Down in DTP": ["S100A6", "HMGA1"],
}
marker_groups = {k: [g for g in v if g in ad_.raw.var_names] for k, v in marker_groups.items()}
rows = []
X = ad_.raw.to_adata()
for grp, genes in marker_groups.items():
    for g in genes:
        x = np.asarray(X[:, g].X.todense()).ravel() if hasattr(X[:, g].X, "todense") else np.asarray(X[:, g].X).ravel()
        for cl in ad_.obs["leiden"].cat.categories:
            m = (ad_.obs["leiden"] == cl).values
            rows.append(dict(group=grp, gene=g, cluster=cl, mean_lognorm=x[m].mean(), pct_expressing=100 * (x[m] > 0).mean(), n_cells=int(m.sum())))
        for tp in TP_ORDER:
            m = (ad_.obs["timepoint"] == tp).values
            rows.append(dict(group=grp, gene=g, cluster=f"tp_{tp}", mean_lognorm=x[m].mean(), pct_expressing=100 * (x[m] > 0).mean(), n_cells=int(m.sum())))
mk = pd.DataFrame(rows)
mk.to_csv(f"{TBL}/S_marker_expression_by_cluster.csv", index=False)
piv_pct = mk[~mk.cluster.str.startswith("tp_")].pivot(index="gene", columns="cluster", values="pct_expressing").round(1)
print("pct expressing per cluster\n", piv_pct.loc[[g for v in marker_groups.values() for g in v]].to_string())
piv_tp = mk[mk.cluster.str.startswith("tp_")].pivot(index="gene", columns="cluster", values="pct_expressing").round(1)[[f"tp_{t}" for t in TP_ORDER]]
print("pct expressing per time point\n", piv_tp.loc[[g for v in marker_groups.values() for g in v]].to_string())

for groupby, fn, w in [("leiden", "FigS_marker_dotplot_cluster", 12.5), ("timepoint", "FigS_marker_dotplot_timepoint", 12.5)]:
    dp = sc.pl.dotplot(ad_, marker_groups, groupby=groupby, use_raw=True, standard_scale="var",
                       dot_max=1.0, figsize=(w, 5.2), show=False, return_fig=True,
                       colorbar_title="Scaled mean\nexpression", size_title="% cells\nexpressing")
    dp.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=0.4)
    dp.savefig(f"{FIG}/{fn}.png", dpi=300, bbox_inches="tight"); dp.savefig(f"{FIG}/{fn}.pdf", bbox_inches="tight")
    plt.close("all")

sc.tl.rank_genes_groups(ad_, "leiden", method="wilcoxon", use_raw=True, pts=True)
top = []
for cl in ad_.obs["leiden"].cat.categories:
    d = sc.get.rank_genes_groups_df(ad_, cl).head(15); d["cluster"] = cl; top.append(d)
top = pd.concat(top); top.to_csv(f"{TBL}/S_cluster_top_markers.csv", index=False)
print(top.groupby("cluster")["names"].apply(lambda s: ", ".join(s.head(10))).to_string())

deg = pd.read_csv(f"{PROC}/DEG_DTP_vs_naive.csv").set_index("names")
def count(d, lfc=1.0, q=0.05):
    return int(((d.logfoldchanges > lfc) & (d.pvals_adj < q)).sum()), int(((d.logfoldchanges < -lfc) & (d.pvals_adj < q)).sum())
grid = []
for lfc in (0.25, 0.5, 1.0, 1.5, 2.0):
    for q in (0.05, 0.01, 0.001):
        u, dn = count(deg, lfc, q); grid.append(dict(method="single-cell Wilcoxon", abs_logFC=lfc, FDR=q, up=u, down=dn))
pb = pd.read_csv(f"{TBL}/contrasts/DE_A_PC9_erlotinib_late_vs_naive.csv").set_index("gene")
for lfc in (0.25, 0.5, 1.0, 1.5, 2.0):
    for q in (0.05, 0.01, 0.001):
        grid.append(dict(method="pseudobulk edgeR", abs_logFC=lfc, FDR=q,
                         up=int(((pb.logFC > lfc) & (pb.FDR < q)).sum()), down=int(((pb.logFC < -lfc) & (pb.FDR < q)).sum())))
grid = pd.DataFrame(grid); grid["up_to_down"] = (grid.up / grid.down.replace(0, np.nan)).round(2)
grid.to_csv(f"{TBL}/S_DE_threshold_grid.csv", index=False)
print(grid.to_string())

sub = ad_[ad_.obs["timepoint"].isin(["D0", "D9", "D11"])].copy()
sub.obs["grp"] = np.where(sub.obs["timepoint"] == "D0", "naive", "DTP")
print(sub.obs.groupby("grp")[["total_counts", "n_genes_by_counts"]].median())

rng = np.random.default_rng(1)
naive_idx = np.where(sub.obs["grp"] == "naive")[0]; dtp_idx = np.where(sub.obs["grp"] == "DTP")[0]
ds = []
for i in range(50):
    pick = np.concatenate([rng.choice(naive_idx, size=len(dtp_idx), replace=False), dtp_idx])
    s = sub[pick].copy()
    sc.tl.rank_genes_groups(s, "grp", groups=["DTP"], reference="naive", method="wilcoxon", use_raw=True)
    d = sc.get.rank_genes_groups_df(s, "DTP")
    ds.append(count(d))
ds = pd.DataFrame(ds, columns=["up", "down"])
print("down-sampled (n=685 vs 685, 50 draws): up median %d [%d-%d], down median %d [%d-%d]" % (
    ds.up.median(), ds.up.min(), ds.up.max(), ds.down.median(), ds.down.min(), ds.down.max()))
ds.to_csv(f"{TBL}/S_DE_downsampling.csv", index=False)

import scipy.sparse as sp
cnt = sub.layers["counts"]; cnt = cnt.toarray() if sp.issparse(cnt) else np.asarray(cnt)
target = 700
keep = cnt.sum(1) >= target
cnt = cnt[keep].astype(int); obs = sub.obs[keep].copy()
rng = np.random.default_rng(2)
thin = np.vstack([rng.multivariate_hypergeometric(c, target) for c in cnt]).astype(np.float32)
import anndata as adm
t = adm.AnnData(thin, obs=obs, var=pd.DataFrame(index=sub.var_names))
sc.pp.normalize_total(t, target_sum=1e4); sc.pp.log1p(t)
sc.tl.rank_genes_groups(t, "grp", groups=["DTP"], reference="naive", method="wilcoxon")
dthin = sc.get.rank_genes_groups_df(t, "DTP")
print("depth-matched (every cell thinned to %d UMIs; n naive=%d, n DTP=%d): up/down =" % (target, (obs.grp == "naive").sum(), (obs.grp == "DTP").sum()), count(dthin))
dthin.to_csv(f"{TBL}/S_DE_depth_matched.csv", index=False)

c = pb.index.intersection(deg.index)
rho = stats.spearmanr(pb.loc[c, "logFC"], deg.loc[c, "logfoldchanges"]).statistic
up = deg[(deg.logfoldchanges > 1) & (deg.pvals_adj < 0.05)].index; dn = deg[(deg.logfoldchanges < -1) & (deg.pvals_adj < 0.05)].index
summ = []
for nm, s, sg in (("single-cell up (759)", up, 1), ("single-cell down (67)", dn, -1)):
    tt = pb.reindex(s).dropna(subset=["logFC"])
    summ.append(dict(set=nm, tested_in_pseudobulk=len(tt), same_direction=int((tt.logFC * sg > 0).sum()),
                     pseudobulk_FDR05_same_direction=int(((tt.FDR < 0.05) & (tt.logFC * sg > 0)).sum()),
                     pseudobulk_P05_same_direction=int(((tt.PValue < 0.05) & (tt.logFC * sg > 0)).sum())))
summ = pd.DataFrame(summ); summ["spearman_logFC_all_genes"] = round(rho, 3)
summ.to_csv(f"{TBL}/S_DE_singlecell_vs_pseudobulk.csv", index=False)
print(summ.to_string())

ad_.obs[["S_score", "G2M_score", "phase"]].to_csv(f"{CACHE}/primary_cellcycle.csv")
print("done")
