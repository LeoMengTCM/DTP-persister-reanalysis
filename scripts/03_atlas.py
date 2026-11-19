#!/usr/bin/env python3
"""Embedding, clustering, cell cycle and marker tables for the erlotinib time course."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import silhouette_score, adjusted_rand_score
set_pub_style(); sc.settings.verbosity = 0

a0 = sc.read_h5ad(f"{PROC}/PC9_DTP_timecourse_final.h5ad")
orig_stats = dict(
    silhouette_timepoint_harmonyPCs=silhouette_score(a0.obsm["X_pca_harmony"][:, :30], a0.obs.timepoint, sample_size=3000, random_state=0),
    ARI_origLeiden_timepoint=adjusted_rand_score(a0.obs.leiden, a0.obs.timepoint),
    rho_pseudotime_day=stats.spearmanr(a0.obs.dpt_pseudotime, a0.obs.day).statistic)

keep_obs = ["sample", "timepoint", "day", "stage", "n_genes_by_counts", "total_counts", "log1p_total_counts", "pct_counts_mt"]
a = sc.AnnData(a0.layers["counts"].copy(), obs=a0.obs[keep_obs].copy(), var=pd.DataFrame(index=a0.var_names))
a.obs["leiden_original"] = a0.obs["leiden"].values
a.obs["timepoint"] = pd.Categorical(a.obs["timepoint"], categories=TP_ORDER, ordered=True)
cc = pd.read_csv(f"{CACHE}/primary_cellcycle.csv", index_col=0)
for c_ in ("S_score", "G2M_score"):
    a.obs[c_] = cc.loc[a.obs_names, c_].astype(float).values
a.obs["phase"] = pd.Categorical(cc.loc[a.obs_names, "phase"].astype(str).values, categories=["G1", "S", "G2M"])
a.layers["counts"] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a); a.raw = a
sc.pp.highly_variable_genes(a, n_top_genes=3000, batch_key="timepoint")
h = a[:, a.var.highly_variable].copy()
sc.pp.regress_out(h, ["log1p_total_counts"])
sc.pp.scale(h, max_value=10); sc.tl.pca(h, n_comps=50, random_state=0)
a.obsm["X_pca"] = h.obsm["X_pca"]
sc.pp.neighbors(a, n_neighbors=20, n_pcs=30, random_state=0)
sc.tl.umap(a, min_dist=0.3, random_state=0)
sc.tl.leiden(a, resolution=0.8, random_state=0, flavor="igraph", n_iterations=2, directed=False, key_added="cluster")

order = a.obs.groupby("cluster", observed=True)["day"].mean().sort_values().index.tolist()
a.obs["cluster"] = pd.Categorical(a.obs["cluster"].map({c: f"C{i+1}" for i, c in enumerate(order)}).astype(str),
                                  categories=[f"C{i+1}" for i in range(len(order))])

new_stats = dict(
    silhouette_timepoint_PCs=silhouette_score(a.obsm["X_pca"][:, :30], a.obs.timepoint, sample_size=3000, random_state=0),
    ARI_cluster_timepoint=adjusted_rand_score(a.obs.cluster, a.obs.timepoint),
    max_abs_rho_PC1to5_depth=max(abs(stats.spearmanr(a.obsm["X_pca"][:, i], a.obs.log1p_total_counts).statistic) for i in range(5)),
    max_abs_rho_PC1to5_depth_original=max(abs(stats.spearmanr(a0.obsm["X_pca_harmony"][:, i], a0.obs.log1p_total_counts).statistic) for i in range(5)))
pd.Series({**orig_stats, **new_stats}).round(3).to_csv(f"{TBL}/S_embedding_diagnostics.csv", header=["value"])
print(pd.Series({**orig_stats, **new_stats}).round(3).to_string())

ct = pd.crosstab(a.obs["cluster"], a.obs["timepoint"]); ct["n_cells"] = ct.sum(1)
frac = ct[TP_ORDER].div(ct["n_cells"], axis=0).round(3).add_prefix("frac_")
phase = pd.crosstab(a.obs["cluster"], a.obs["phase"], normalize="index").round(3).add_prefix("phase_")
qc = a.obs.groupby("cluster", observed=True)[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].median().round(1).add_prefix("median_")
comp = pd.concat([ct, frac, phase, qc], axis=1)
comp.to_csv(f"{TBL}/Table_cluster_composition.csv"); print(comp.to_string())
share = pd.crosstab(a.obs["timepoint"], a.obs["cluster"], normalize="index").round(3)
share.to_csv(f"{TBL}/Table_timepoint_cluster_share.csv"); print(share.to_string())

sc.tl.rank_genes_groups(a, "cluster", method="wilcoxon", use_raw=True, pts=True)
top = pd.concat([sc.get.rank_genes_groups_df(a, c).head(25).assign(cluster=c) for c in a.obs["cluster"].cat.categories])
top.to_csv(f"{TBL}/S_cluster_top_markers.csv", index=False)
print(top.groupby("cluster")["names"].apply(lambda s: ", ".join(s.head(12))).to_string())

marker_groups = {
    "Proliferation": ["MKI67", "TOP2A", "PCNA"],
    "Stress response": ["ATF3", "DDIT3", "JUN"],
    "Persister / EMT": ["AXL", "IGFBP3", "VIM", "FN1", "SERPINE1", "CYP1B1"],
    "Iron / redox": ["FTL", "FTH1", "ALDH1A3", "GPX4"],
    "Down in DTP": ["S100A6", "HMGA1"],
}
rows = []
R = a.raw.to_adata()
for grp, genes in marker_groups.items():
    for g in genes:
        x = R[:, g].X; x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
        for key in ("cluster", "timepoint"):
            for lev in a.obs[key].cat.categories:
                m = (a.obs[key] == lev).values
                rows.append(dict(marker_group=grp, gene=g, grouping=key, level=lev, n_cells=int(m.sum()),
                                 mean_lognorm=round(float(x[m].mean()), 4), pct_expressing=round(100 * float((x[m] > 0).mean()), 1)))
mk = pd.DataFrame(rows); mk.to_csv(f"{TBL}/Table_marker_expression.csv", index=False)
for key in ("cluster", "timepoint"):
    p = mk[mk.grouping == key].pivot(index="gene", columns="level", values="pct_expressing")
    p = p.loc[[g for v in marker_groups.values() for g in v], list(a.obs[key].cat.categories)]
    print(f"% expressing by {key}\n", p.to_string())

for key, fn in (("cluster", "panel_marker_dotplot_cluster"), ("timepoint", "panel_marker_dotplot_timepoint")):
    dp = sc.pl.dotplot(a, marker_groups, groupby=key, use_raw=True, standard_scale="var", dot_max=1.0,
                       figsize=(13, 4.6), show=False, return_fig=True,
                       colorbar_title="Scaled mean\nexpression", size_title="% cells\nexpressing")
    dp.style(cmap="Reds", dot_edge_color="black", dot_edge_lw=0.4)
    dp.savefig(f"{FIG}/{fn}.png", dpi=300, bbox_inches="tight"); dp.savefig(f"{FIG}/{fn}.pdf", bbox_inches="tight"); plt.close("all")

a.write_h5ad(f"{CACHE}/primary_atlas.h5ad")

pb = pd.read_csv(f"{CACHE}/pseudobulk_counts.csv", index_col=0)
cpm = np.log2(pb / pb.sum() * 1e6 + 1)
expr = cpm[(pb > 0).all(axis=1) & (cpm.mean(axis=1) > np.log2(11))]
hv = expr.loc[expr.var(axis=1).sort_values(ascending=False).index[:2000]]
corr = hv.corr(method="spearman").round(3)
corr.to_csv(f"{TBL}/S_library_correlation.csv"); print(corr.to_string())
meta = pd.read_csv(f"{CACHE}/pseudobulk_meta.csv").set_index("sample")
pairs = []
libs = list(corr.index)
for i in range(len(libs)):
    for j in range(i + 1, len(libs)):
        t1, t2 = meta.loc[libs[i], "timepoint"], meta.loc[libs[j], "timepoint"]
        g = lambda t: "naive" if t == "D0" else ("late" if t in ("D9", "D11") else "early")
        pairs.append(dict(a=libs[i], b=libs[j], same_timepoint=t1 == t2, same_stage=g(t1) == g(t2), rho=corr.iloc[i, j]))
pairs = pd.DataFrame(pairs); pairs.to_csv(f"{TBL}/S_library_correlation_pairs.csv", index=False)
print(pairs.groupby("same_timepoint").rho.agg(["mean", "min", "max", "size"]).round(3))
print(pairs.groupby("same_stage").rho.agg(["mean", "min", "max", "size"]).round(3))
import scipy.cluster.hierarchy as hc
Z = hc.linkage(hv.T.values, "average", metric="correlation")
print("hierarchical clustering order:", [libs[i] for i in hc.leaves_list(Z)])
print("done")
