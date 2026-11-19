#!/usr/bin/env python3
"""Persister signatures vs CRISPR gene effect in DepMap, lineage-adjusted."""
import sys, os, glob
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

model = pd.read_csv(f"{DEPMAP}/Model.csv").set_index("ModelID")
crispr = pd.read_csv(f"{DEPMAP}/CRISPRGeneEffect.csv", index_col=0)
crispr.columns = [c.split(" (")[0] for c in crispr.columns]
crispr = crispr.loc[:, ~crispr.columns.duplicated()]
expr = pd.read_csv(f"{DEPMAP}/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", index_col=0)
expr = expr[expr["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")
expr = expr.drop(columns=["SequencingID", "ModelConditionID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"])
expr.columns = [c.split(" (")[0] for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()]
common = crispr.index.intersection(expr.index).intersection(model.index)
print(f"CRISPR lines {crispr.shape[0]} x {crispr.shape[1]} genes | expression models {expr.shape[0]} | both {len(common)}")
crispr, expr, meta = crispr.loc[common], expr.loc[common], model.loc[common]
lin = meta["OncotreeLineage"].fillna("Unknown")
vc = lin.value_counts(); lin = lin.where(lin.map(vc) >= 5, "Other_small")
print("lineages used as fixed effects:", lin.nunique()); print(vc.head(12).to_string())
nsclc = (meta["OncotreeLineage"] == "Lung") & meta["OncotreePrimaryDisease"].fillna("").str.contains("Non-Small Cell")
print("NSCLC lines with both data types:", int(nsclc.sum()))

dtp_contrasts = ["A_PC9_erlotinib_late_vs_naive", "B_PC9_gef_osi_9d_vs_DMSO", "C_PC9_DTP_vs_DMSO", "C_HCC827_DTP_vs_DMSO",
                 "C_HCC2935_DTP_vs_DMSO", "C_H1975_DTP_vs_DMSO"]
DE = {c: pd.read_csv(f"{TBL}/contrasts/DE_{c}.csv").set_index("gene") for c in dtp_contrasts}
up = pd.DataFrame({c: (d.FDR < 0.05) & (d.logFC >= 1) for c, d in DE.items()}).fillna(False)
dn = pd.DataFrame({c: (d.FDR < 0.05) & (d.logFC <= -1) for c, d in DE.items()}).fillna(False)
n_up, n_dn = up.sum(1), dn.sum(1)
print("genes up in >=k of 6 contrasts:", {k: int((n_up >= k).sum()) for k in range(1, 7)})
print("genes down in >=k of 6 contrasts:", {k: int((n_dn >= k).sum()) for k in range(1, 7)})
excl = all_metabolic_genes()
cons_all = n_up[n_up >= 4].index.tolist()
cons = [g for g in cons_all if g not in excl and g in expr.columns]
print(f"consensus-up genes (>=4/6): {len(cons_all)}; removed as metabolic/candidate: {sorted(set(cons_all) & excl)}; scored in DepMap: {len(cons)}")
pc9_all = DE[dtp_contrasts[0]].query("FDR < 0.05 and logFC >= 1").index.tolist()
pc9 = [g for g in pc9_all if g not in excl and g in expr.columns]
print(f"PC9-erlotinib-up genes: {len(pc9_all)}; removed: {sorted(set(pc9_all) & excl)}; scored: {len(pc9)}")
cons_dn = [g for g in n_dn[n_dn >= 4].index if g not in excl and g in expr.columns]
hall = pd.read_csv(f"{CACHE}/net_hallmark.csv")
emt = [g for g in hall[hall.source.str.contains("EPITHELIAL_MESENCHYMAL")].target.unique() if g in expr.columns]
sigs = {"consensus": cons, "pc9_erl": pc9, "orig9": ORIG_DTP_SIGNATURE, "orig6": [g for g in ORIG_DTP_SIGNATURE if g not in ("GPX4", "SLC7A11", "GCLC")],
        "hallmark_EMT": emt, "consensus_down": cons_dn}
pd.DataFrame([(k, g, int(n_up.get(g, 0))) for k in ("consensus", "pc9_erl") for g in sigs[k]], columns=["signature", "gene", "n_contrasts_up"]).to_csv(f"{TBL}/Table_persister_signatures.csv", index=False)
print("consensus signature genes:", ", ".join(sorted(cons)))
print("overlap consensus x Hallmark EMT:", len(set(cons) & set(emt)), "of", len(cons))

z = (expr - expr.mean()) / expr.std(ddof=0).replace(0, np.nan)
score = pd.DataFrame({k: z[[g for g in v if g in z.columns]].mean(1) for k, v in sigs.items()})
score["consensus_up_minus_down"] = score["consensus"] - score["consensus_down"]
score["orig9_meanTPM"] = expr[ORIG_DTP_SIGNATURE].mean(1)
score.assign(lineage=lin, NSCLC=nsclc).to_csv(f"{TBL}/S_depmap_signature_scores.csv")
print(score.corr(method="spearman").round(2).to_string())
print("signature score by lineage (consensus, median):\n", score.groupby(lin)["consensus"].median().sort_values().round(2).to_string())

def assoc(sig, covar=None, lines=None, fixed=True):
    idx = score.index if lines is None else lines
    y = crispr.loc[idx]; s = score.loc[idx, sig].values.astype(float)
    cols = [np.ones(len(idx))]
    if fixed:
        d = pd.get_dummies(lin.loc[idx], drop_first=True).astype(float)
        d = d.loc[:, d.sum() > 0]; cols += [d[c].values for c in d.columns]
    if covar is not None:
        for c in covar: cols.append(score.loc[idx, c].values.astype(float))
    X = np.column_stack(cols)
    Y = y.values.astype(float); nanmask = np.isnan(Y)
    Y = np.where(nanmask, np.nanmean(Y, axis=0), Y)
    P = X @ np.linalg.pinv(X)
    sr = s - P @ s; Yr = Y - P @ Y
    dfree = len(idx) - np.linalg.matrix_rank(X) - 1
    beta = (sr @ Yr) / (sr @ sr)
    resid = Yr - np.outer(sr, beta)
    se = np.sqrt((resid ** 2).sum(0) / dfree / (sr @ sr))
    t = beta / se; p = 2 * stats.t.sf(np.abs(t), dfree)
    r = t / np.sqrt(t ** 2 + dfree)
    out = pd.DataFrame({"beta_per_SD": beta * s.std(), "partial_r": r, "t": t, "P": p}, index=y.columns)
    out["FDR_genomewide"] = multipletests(out.P.fillna(1), method="fdr_bh")[1]
    out["rank_most_negative"] = out.t.rank(method="min")
    out["pct_rank"] = 100 * out["rank_most_negative"] / len(out)
    out["n_lines"] = len(idx); out["median_effect"] = np.nanmedian(y.values, axis=0)
    return out

cand = [g for g in CANDIDATES if g in crispr.columns]
runs = {
    "consensus | lineage-adjusted": assoc("consensus"),
    "consensus | unadjusted": assoc("consensus", fixed=False),
    "consensus | lineage + Hallmark EMT adjusted": assoc("consensus", covar=["hallmark_EMT"]),
    "consensus | NSCLC only": assoc("consensus", lines=score.index[nsclc.values], fixed=False),
    "consensus | solid tumours, lineage-adjusted": assoc("consensus", lines=score.index[~lin.isin(["Lymphoid", "Myeloid"]).values]),
    "pc9_erl | lineage-adjusted": assoc("pc9_erl"),
    "consensus_up_minus_down | lineage-adjusted": assoc("consensus_up_minus_down"),
    "hallmark_EMT | lineage-adjusted": assoc("hallmark_EMT"),
    "orig6 (no candidate genes) | lineage-adjusted": assoc("orig6"),
    "orig9 (original, with GPX4/SLC7A11/GCLC) | lineage-adjusted": assoc("orig9"),
    "orig9 (original) | unadjusted": assoc("orig9", fixed=False),
}
rows = []
for name, r in runs.items():
    d = r.loc[cand].copy(); d["analysis"] = name; d["gene"] = d.index; rows.append(d)
    n_sig_neg = int(((r.FDR_genomewide < 0.05) & (r.t < 0)).sum()); n_sig_pos = int(((r.FDR_genomewide < 0.05) & (r.t > 0)).sum())
    print(f"\n== {name} (n={int(r.n_lines.iloc[0])} lines, {len(r)} genes; FDR<0.05: {n_sig_neg} more essential / {n_sig_pos} less essential in signature-high lines)")
    print(d.sort_values("t")[["partial_r", "t", "P", "FDR_genomewide", "rank_most_negative", "pct_rank", "median_effect"]].head(8).round(4).to_string())
    print("  genome-wide top 15 (more essential with higher signature):", ", ".join(r.sort_values("t").head(15).index))
res = pd.concat(rows); res.to_csv(f"{TBL}/Table_depmap_candidate_associations.csv", index=False)
runs["consensus | lineage-adjusted"].sort_values("t").to_csv(f"{TBL}/S_depmap_genomewide_consensus_lineage_adjusted.csv")
runs["consensus | lineage + Hallmark EMT adjusted"].sort_values("t").to_csv(f"{TBL}/S_depmap_genomewide_consensus_lineage_EMT_adjusted.csv")

key = ["GPX4", "NFE2L2", "PHGDH", "CPT1A", "TXNRD1", "SOD2", "DHFR", "SLC7A11", "GCLC", "NQO1"]
piv_r = res[res.gene.isin(key)].pivot(index="gene", columns="analysis", values="partial_r").loc[key, list(runs)]
piv_q = res[res.gene.isin(key)].pivot(index="gene", columns="analysis", values="FDR_genomewide").loc[key, list(runs)]
pd.set_option("display.width", 330); pd.set_option("display.max_columns", 30)
print("\npartial r\n", piv_r.round(3).T.to_string()); print("\ngenome-wide FDR\n", piv_q.applymap(lambda v: f"{v:.1e}").T.to_string())

base = runs["consensus | lineage-adjusted"].loc["GPX4", "t"]
loo = []
for g in cons:
    score["_loo"] = z[[x for x in cons if x != g]].mean(1)
    loo.append(assoc("_loo").loc["GPX4", "t"])
print(f"\nGPX4 t with full consensus signature {base:.2f}; leave-one-gene-out range {min(loo):.2f} to {max(loo):.2f}")
rng = np.random.default_rng(0)
expressed = expr.columns[(expr.mean() > 1)]
rand_t = []
for i in range(500):
    score["_rnd"] = z[rng.choice(expressed, size=len(cons), replace=False)].mean(1)
    rand_t.append(assoc("_rnd").loc["GPX4", "t"])
rand_t = np.array(rand_t)
print(f"random size-matched signatures (500): GPX4 t mean {rand_t.mean():.2f}, sd {rand_t.std():.2f}; empirical P(t <= observed) = {(np.sum(rand_t <= base) + 1) / 501:.4f}")
pd.DataFrame({"random_signature_GPX4_t": rand_t}).to_csv(f"{TBL}/S_depmap_random_signature_GPX4.csv", index=False)

med = score["orig9_meanTPM"].median()
hi = score.index[score["orig9_meanTPM"] >= med]; lo = score.index[score["orig9_meanTPM"] < med]
rep = []
for g in ["GPX4", "NFE2L2", "PHGDH", "CPT1A", "NQO1", "SLC7A11", "TXNRD1", "GCLC", "SOD2", "DHFR", "PSAT1", "ALDH1A3"]:
    a, b = crispr.loc[hi, g].dropna(), crispr.loc[lo, g].dropna()
    rep.append(dict(gene=g, n_high=len(a), n_low=len(b), median_high=a.median(), median_low=b.median(), MWU_P=stats.mannwhitneyu(a, b).pvalue))
rep = pd.DataFrame(rep); rep.to_csv(f"{TBL}/S_depmap_original_median_split_reproduced.csv", index=False)
print("\noriginal 9-gene median split reproduced:\n", rep.round(4).to_string())

med = score["consensus"].median(); hi = score.index[score["consensus"] >= med]; lo = score.index[score["consensus"] < med]
a, b = crispr.loc[hi, "GPX4"].dropna(), crispr.loc[lo, "GPX4"].dropna()
print(f"consensus signature median split, GPX4: n {len(a)}/{len(b)}, median {a.median():.3f} vs {b.median():.3f}, MWU P={stats.mannwhitneyu(a, b).pvalue:.2e}")
print("done")
