#!/usr/bin/env python3
"""Pathway scores per treatment day, raw and depth-matched."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests
sc.settings.verbosity = 0
rng = np.random.default_rng(7)

a = sc.read_h5ad(f"{CACHE}/primary_atlas.h5ad")
a.obs["timepoint"] = pd.Categorical(a.obs["timepoint"], categories=TP_ORDER, ordered=True)
SETS = dict(METABOLIC_PATHWAYS); SETS["Defense_score_9genes"] = DEFENSE_GENES

def score_all(adata, use_raw):
    names = adata.raw.var_names if use_raw else adata.var_names
    used = {}
    for k, g in SETS.items():
        gg = resolve(g, names); used[k] = gg
        sc.tl.score_genes(adata, gg, score_name=f"s_{k}", use_raw=use_raw, random_state=0)
    return used
used = score_all(a, True)
pd.DataFrame([(k, len(SETS[k]), len(v), ", ".join(v), ", ".join(sorted(set(SETS[k]) - set(v) - {x for x in SETS[k] if ALIASES.get(x) in v})))
              for k, v in used.items()], columns=["set", "n_curated", "n_detected", "genes_scored", "genes_not_detected"]).to_csv(f"{TBL}/Table_gene_sets_scored.csv", index=False)

cnt = a.layers["counts"]; cnt = cnt.toarray() if sp.issparse(cnt) else np.asarray(cnt)
TARGET = 700
ok = cnt.sum(1) >= TARGET
thin = np.vstack([rng.multivariate_hypergeometric(c, TARGET) for c in cnt[ok].astype(int)]).astype(np.float32)
t = ad.AnnData(thin, obs=a.obs[ok][["timepoint", "day", "cluster", "S_score", "G2M_score"]].copy(), var=pd.DataFrame(index=a.var_names))
sc.pp.normalize_total(t, target_sum=1e4); sc.pp.log1p(t); score_all(t, False)
n_tab = pd.DataFrame({"all_cells": a.obs.timepoint.value_counts().reindex(TP_ORDER), "depth_matched_cells": t.obs.timepoint.value_counts().reindex(TP_ORDER),
                      "median_UMI": a.obs.groupby("timepoint", observed=True).total_counts.median().reindex(TP_ORDER)})
print(n_tab.to_string()); n_tab.to_csv(f"{TBL}/S_cells_per_timepoint.csv")

def boot_ci(x, y, n=1000):
    ds = [cohens_d(rng.choice(x, len(x)), rng.choice(y, len(y))) for _ in range(n)]
    return np.percentile(ds, [2.5, 97.5])

rows = []
for k in SETS:
    col = f"s_{k}"
    kw = stats.kruskal(*[a.obs.loc[a.obs.timepoint == tp, col].values for tp in TP_ORDER])
    kw_t = stats.kruskal(*[t.obs.loc[t.obs.timepoint == tp, col].values for tp in TP_ORDER])

    X = np.column_stack([np.ones(a.n_obs), a.obs.S_score.values, a.obs.G2M_score.values])
    res = a.obs[col].values - X @ np.linalg.lstsq(X, a.obs[col].values, rcond=None)[0]
    d0, d0t, r0 = a.obs.loc[a.obs.timepoint == "D0", col].values, t.obs.loc[t.obs.timepoint == "D0", col].values, res[(a.obs.timepoint == "D0").values]
    groups = {tp: [tp] for tp in TP_ORDER[1:]}; groups["D9+D11"] = ["D9", "D11"]
    for lab, tps in groups.items():
        m = a.obs.timepoint.isin(tps).values; mt = t.obs.timepoint.isin(tps).values
        x, xt, xr = a.obs.loc[m, col].values, t.obs.loc[mt, col].values, res[m]
        lo, hi = boot_ci(x, d0)

        sub = [cohens_d(rng.choice(x, min(200, len(x)), replace=False), rng.choice(d0, 200, replace=False)) for _ in range(200)]
        rows.append(dict(set=k, comparison=f"{lab} vs D0", n_group=int(m.sum()), n_D0=len(d0), mean_group=x.mean(), mean_D0=d0.mean(),
                         cohens_d=cohens_d(x, d0), ci_low=lo, ci_high=hi, MWU_P=stats.mannwhitneyu(x, d0).pvalue,
                         d_depth_matched=cohens_d(xt, d0t), MWU_P_depth_matched=stats.mannwhitneyu(xt, d0t).pvalue,
                         n_group_depth_matched=int(mt.sum()), n_D0_depth_matched=len(d0t),
                         d_cellcycle_adjusted=cohens_d(xr, r0), d_n200_p2_5=np.percentile(sub, 2.5), d_n200_p97_5=np.percentile(sub, 97.5),
                         KW_H=kw.statistic, KW_P=kw.pvalue, KW_P_depth_matched=kw_t.pvalue))
R = pd.DataFrame(rows)
R["MWU_FDR"] = multipletests(R.MWU_P, method="fdr_bh")[1]; R["MWU_FDR_depth_matched"] = multipletests(R.MWU_P_depth_matched, method="fdr_bh")[1]
R.to_csv(f"{TBL}/Table_pathway_timepoint_singlecell.csv", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 30)
for what in ("cohens_d", "d_depth_matched", "d_cellcycle_adjusted"):
    print(f"\n{what} vs D0"); print(R.pivot(index="set", columns="comparison", values=what)[[f"{tp} vs D0" for tp in TP_ORDER[1:]] + ["D9+D11 vs D0"]].round(2).to_string())
sig = R.assign(sig=(R.MWU_FDR < 0.05) & (R.MWU_FDR_depth_matched < 0.05) & (np.sign(R.cohens_d) == np.sign(R.d_depth_matched)))
print("\nsignificant (FDR<0.05) in both raw and depth-matched, same sign:"); print(sig.pivot(index="set", columns="comparison", values="sig")[[f"{tp} vs D0" for tp in TP_ORDER[1:]] + ["D9+D11 vs D0"]].replace({True: "yes", False: "."}).to_string())
print("\nKruskal-Wallis across six time points: raw P<0.05 in", int((R.drop_duplicates('set').query("set != 'Defense_score_9genes'").KW_P < 0.05).sum()),
      "of 15; depth-matched P<0.05 in", int((R.drop_duplicates('set').query("set != 'Defense_score_9genes'").KW_P_depth_matched < 0.05).sum()), "of 15")

ds = pd.concat([a.obs.groupby("timepoint", observed=True)["s_Defense_score_9genes"].agg(["mean", "std", "size"]).add_prefix("raw_"),
                t.obs.groupby("timepoint", observed=True)["s_Defense_score_9genes"].agg(["mean", "std", "size"]).add_prefix("depthmatched_")], axis=1)
print("\nDefense score by time point\n", ds.round(3).to_string()); ds.to_csv(f"{TBL}/Table_defense_score_timepoint.csv")
gr = []
RA = a.raw.to_adata()
for g in resolve(DEFENSE_GENES, a.raw.var_names) + ["FTH1", "FTL"]:
    x = np.asarray(RA[:, g].X.todense()).ravel() if sp.issparse(RA[:, g].X) else np.asarray(RA[:, g].X).ravel()
    xt = np.asarray(t[:, g].X).ravel()
    for tp in TP_ORDER:
        m = (a.obs.timepoint == tp).values; mt = (t.obs.timepoint == tp).values
        gr.append(dict(gene=g, timepoint=tp, mean_raw=x[m].mean(), pct_raw=100 * (x[m] > 0).mean(), mean_depth_matched=xt[mt].mean(), pct_depth_matched=100 * (xt[mt] > 0).mean()))
gr = pd.DataFrame(gr); gr.to_csv(f"{TBL}/Table_defense_genes_timepoint.csv", index=False)
print("\n% cells expressing, raw\n", gr.pivot(index="gene", columns="timepoint", values="pct_raw")[TP_ORDER].round(1).to_string())
print("\n% cells expressing, depth-matched (700 UMIs/cell)\n", gr.pivot(index="gene", columns="timepoint", values="pct_depth_matched")[TP_ORDER].round(1).to_string())

lc = pd.read_csv(f"{CACHE}/pseudobulk_logCPM_all_libraries.csv", index_col=0)
ref = lc[["D0_a", "D0_b"]].mean(1)
lib = pd.DataFrame({k: (lc.loc[resolve(g, lc.index)].sub(ref.loc[resolve(g, lc.index)], axis=0)).mean() for k, g in SETS.items()}).T
lib = lib[["D0_a", "D0_b", "D1_a", "D2_a", "D3_a", "D3_b", "D4_a", "D9_a", "D11_a", "D11_b", "D11_c"]]
lib.to_csv(f"{TBL}/Table_pathway_library_level.csv"); print("\nlibrary-level mean log2FC vs mean of D0 libraries\n", lib.round(2).to_string())

late = a.obs.cluster.isin(["C4", "C5"])
print("\nD0 cells in clusters C4/C5:", int((late & (a.obs.timepoint == "D0")).sum()), "of", int((a.obs.timepoint == "D0").sum()))
sigg = pd.read_csv(f"{TBL}/Table_persister_signatures.csv"); pg = [g for g in sigg[sigg.signature == "consensus"].gene if g in t.var_names]
sc.tl.score_genes(t, pg, score_name="s_persister", use_raw=False, random_state=0)
thr = t.obs.loc[t.obs.timepoint.isin(["D9", "D11"]), "s_persister"].median()
frac = t.obs.groupby("timepoint", observed=True)["s_persister"].apply(lambda s: 100 * (s >= thr).mean()).round(1)
print(f"consensus persister-up score ({len(pg)} genes, depth-matched cells): % of cells at or above the D9/D11 median\n", frac.to_string())
frac.to_csv(f"{TBL}/S_persister_like_fraction_by_timepoint.csv")
a.obs.to_csv(f"{CACHE}/primary_obs_scores.csv"); t.obs.to_csv(f"{CACHE}/primary_obs_scores_depth_matched.csv")
print("done")
