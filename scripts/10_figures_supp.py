#!/usr/bin/env python3
"""Figures S1-S3."""
import sys, os, importlib.util
sys.path.insert(0, os.path.dirname(__file__))
spec = importlib.util.spec_from_file_location("figmain", os.path.join(os.path.dirname(__file__), "09_figures_main.py"))
src = open(spec.origin).read().split("if __name__")[0]

exec(src.split("def figure1():")[0])
import json, re
_save0 = save
def save(fig, name):

    fig.canvas.draw()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                if not (lo - 1e-9 <= tk.get_loc() <= hi + 1e-9): tk.label1.set_visible(False); tk.label2.set_visible(False)
    r = fig.canvas.get_renderer()
    for ax in fig.axes:
        for lab in ax.get_xticklabels() + ax.get_yticklabels():
            bb = lab.get_window_extent(r)
            if bb.x1 > fig.bbox.x1 + 1 or bb.x0 < -1 or bb.y1 > fig.bbox.y1 + 1 or bb.y0 < -1:
                print("  hidden off-canvas tick label:", repr(lab.get_text()), "axis ylabel:", ax.get_ylabel()[:30], "ylim:", ax.get_ylim()); lab.set_visible(False)
    _save0(fig, name)

def bars(ax, labels, vals, colors, ylabel, title=None, fmt="{:.0f}", rot=0):
    x = np.arange(len(labels)); ax.bar(x, vals, color=colors, width=0.7)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=rot); ax.set_ylabel(ylabel)
    if title: ax.set_title(title)
    if fmt:
        for xi, v in zip(x, vals): ax.text(xi, v, fmt.format(v), ha="center", va="bottom", fontsize=7)
    ax.margins(y=0.15)

def figS1():
    fig = plt.figure(figsize=(W, 7.6), layout="constrained")
    gs = fig.add_gridspec(3, 2)
    qc = pd.read_csv(f"{TBL}/S_QC_per_sample.csv")
    col = [TP_COLORS.get(t, "#8c6bb1") for t in qc.timepoint]
    ax = fig.add_subplot(gs[0, 0]); bars(ax, qc["sample"], qc.cells_passing_QC, col, "Cells passing QC", rot=90); panel(ax, "a")
    ax = fig.add_subplot(gs[0, 1]); bars(ax, qc["sample"], qc.median_UMI, col, "Median UMIs per cell", rot=90); panel(ax, "b")
    cc = pd.read_csv(f"{TBL}/S_cellcycle_by_timepoint.csv").set_index("timepoint").loc[TP_ORDER]
    ax = fig.add_subplot(gs[1, 0]); bottom = np.zeros(len(cc))
    for ph, c in (("G1", "#bdbdbd"), ("S", "#fdae61"), ("G2M", "#2c7bb6")):
        ax.bar(cc.index, 100 * cc[ph], bottom=bottom, color=c, label=ph, width=0.7); bottom += 100 * cc[ph].values
    ax.set_ylabel("Cells (%)"); ax.set_xlabel("Erlotinib treatment day"); ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0)); panel(ax, "c")
    ed = pd.read_csv(f"{TBL}/S_embedding_diagnostics.csv", index_col=0)["value"]
    ax = fig.add_subplot(gs[1, 1]); x = np.arange(3); w = 0.36
    o = [ed["ARI_origLeiden_timepoint"], ed["silhouette_timepoint_harmonyPCs"], ed["max_abs_rho_PC1to5_depth_original"]]
    n = [ed["ARI_cluster_timepoint"], ed["silhouette_timepoint_PCs"], ed["max_abs_rho_PC1to5_depth"]]
    ax.bar(x - w / 2, o, w, color="#bdbdbd", label="Batch-corrected by time point"); ax.bar(x + w / 2, n, w, color="#3b6fb6", label="No batch correction")
    for xi, v in zip(x - w / 2, o): ax.text(xi, max(v, 0), f"{v:.2f}", ha="center", va="bottom", fontsize=7)
    for xi, v in zip(x + w / 2, n): ax.text(xi, max(v, 0), f"{v:.2f}", ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color="k", lw=0.6); ax.set_xticks(x)
    ax.set_xticklabels(["ARI\n(clusters vs day)", "Silhouette\n(day)", "max |rho|\n(PC1-5 vs depth)"]); ax.set_ylabel("Value"); ax.set_ylim(-0.1, 1.25)
    ax.legend(loc="upper left"); panel(ax, "d")
    lc = pd.read_csv(f"{TBL}/S_library_correlation.csv", index_col=0)
    ax = fig.add_subplot(gs[2, 0]); im = ax.imshow(lc.values, cmap="viridis", vmin=0.2, vmax=1)
    ax.set_xticks(range(len(lc))); ax.set_yticks(range(len(lc))); ax.set_xticklabels(lc.columns, rotation=90, fontsize=7); ax.set_yticklabels(lc.index, fontsize=7)
    cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label("Spearman rho"); panel(ax, "e")
    fr = pd.read_csv(f"{TBL}/S_persister_like_fraction_by_timepoint.csv", index_col=0).iloc[:, 0].reindex(TP_ORDER)
    ax = fig.add_subplot(gs[2, 1]); bars(ax, fr.index, fr.values, [TP_COLORS[t] for t in fr.index], "Cells with persister-like\nscore (%)", fmt="{:.0f}")
    ax.set_xlabel("Erlotinib treatment day"); panel(ax, "f")
    save(fig, "Figure_S1")

def figS2():
    fig = plt.figure(figsize=(W, 5.6), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    g = pd.read_csv(f"{TBL}/S_DE_threshold_grid.csv")
    ax = fig.add_subplot(gs[0, 0])
    for m, c in zip(g.method.unique(), ("#d73027", "#3b6fb6", "#7f7f7f", "#f4a259")):
        d = g[(g.method == m) & (g.FDR == 0.05)]
        ax.plot(d.abs_logFC, np.log2(d.up / d.down.replace(0, np.nan)), "o-", color=c, label=m, ms=4)
    ax.axhline(0, color="k", lw=0.6, ls=":"); ax.set_xlabel("|log2 fold change| threshold"); ax.set_ylabel("log2 (up : down), FDR < 0.05")
    ax.legend(loc="upper right", fontsize=7); panel(ax, "a")
    ds = pd.read_csv(f"{TBL}/S_DE_downsampling.csv")
    ax = fig.add_subplot(gs[0, 1]); ax.boxplot([ds.up, ds.down], widths=0.5, showfliers=False)
    for i, k in enumerate(("up", "down")): ax.scatter(np.random.default_rng(1).normal(i + 1, 0.05, len(ds)), ds[k], s=6, color="0.4", zorder=3)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Up", "Down"]); ax.set_ylabel("Genes (685 vs 685 cells, 50 draws)"); ax.set_ylim(0, 800); ax.set_yticks([0, 200, 400, 600, 800]); panel(ax, "b")
    scde = pd.read_csv(f"{PROC}/DEG_DTP_vs_naive.csv").set_index("names")["logfoldchanges"]
    pb = pd.read_csv(f"{TBL}/contrasts/DE_A_PC9_erlotinib_late_vs_naive.csv").set_index("gene")["logFC"]
    d = pd.concat([scde.rename("sc"), pb.rename("pb")], axis=1).dropna(); d = d[np.isfinite(d).all(1)]
    ax = fig.add_subplot(gs[1, 0])
    ax.scatter(d.sc.clip(-8, 8), d.pb.clip(-8, 8), s=2, color="0.45", alpha=0.4, rasterized=True)
    rho = stats.spearmanr(d.sc, d.pb)[0]; ax.text(0.04, 0.97, f"Spearman rho = {rho:.2f}\nn = {len(d):,} genes", transform=ax.transAxes, va="top", fontsize=7.5)
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("log2FC, single-cell Wilcoxon"); ax.set_ylabel("log2FC, pseudobulk edgeR"); panel(ax, "c")
    cs = pd.read_csv(f"{TBL}/Table_contrast_summary.csv")
    keep = ["A_PC9_erlotinib_late_vs_naive", "A_PC9_erlotinib_D3_vs_naive", "A_PC9_erlotinib_late_vs_D3", "B_PC9_gef_osi_9d_vs_DMSO", "C_PC9_DTP_vs_DMSO",
            "C_HCC827_DTP_vs_DMSO", "C_HCC2935_DTP_vs_DMSO", "C_H1975_DTP_vs_DMSO"]
    lab = ["PC9 erl D9-11\nvs D0", "PC9 erl D3\nvs D0", "PC9 erl D9-11\nvs D3", "PC9 gef/osi\n9 d", "PC9 osi\nDTP", "HCC827 osi\nDTP", "HCC2935 osi\nDTP", "H1975 osi\nDTP"]
    d = cs.set_index("contrast").loc[keep]
    ax = fig.add_subplot(gs[1, 1]); x = np.arange(len(keep))
    ax.bar(x - 0.2, d.up_FDR05_lfc1, 0.4, color="#d73027", label="Up"); ax.bar(x + 0.2, d.down_FDR05_lfc1, 0.4, color="#3b6fb6", label="Down")
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=90, fontsize=7); ax.set_ylabel("Genes (FDR < 0.05, |log2FC| >= 1)"); ax.legend(loc="upper left"); panel(ax, "d")
    save(fig, "Figure_S2")

def figS3():
    fig = plt.figure(figsize=(W, 6.0), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    rt = pd.read_csv(f"{TBL}/S_depmap_random_signature_GPX4.csv").iloc[:, 0]
    gw = pd.read_csv(f"{TBL}/S_depmap_genomewide_consensus_lineage_adjusted.csv", index_col=0)
    obs = gw.loc["GPX4", "t"]
    ax = fig.add_subplot(gs[0, 0]); ax.hist(rt, bins=30, color="#bdbdbd"); ax.axvline(obs, color="#d73027", lw=1.5)
    ax.text(obs + 0.1, ax.get_ylim()[1] * 0.9, f"observed\nt = {obs:.2f}", color="#d73027", fontsize=7.5, va="top")
    ax.set_xlabel("GPX4 t statistic, 500 random\nsize-matched signatures"); ax.set_ylabel("Signatures"); panel(ax, "a")
    ca = pd.read_csv(f"{TBL}/Table_depmap_candidate_associations.csv"); d = ca[ca.gene == "GPX4"].set_index("analysis")
    order = ["consensus | unadjusted", "consensus | lineage-adjusted", "consensus | lineage + Hallmark EMT adjusted", "consensus | solid tumours, lineage-adjusted",
             "consensus | NSCLC only", "pc9_erl | lineage-adjusted", "orig6 (no candidate genes) | lineage-adjusted", "orig9 (original, with GPX4/SLC7A11/GCLC) | lineage-adjusted"]
    lab = ["Consensus, unadjusted", "Consensus, lineage", "Consensus, lineage + EMT", "Consensus, solid tumours", "Consensus, NSCLC only", "PC9-erlotinib signature",
           "Literature, 6 genes", "Literature, 9 genes"]
    ax = fig.add_subplot(gs[0, 1]); y = np.arange(len(order))[::-1]; r = d.loc[order, "partial_r"].values; n = d.loc[order, "n_lines"].values
    se = np.sqrt((1 - r ** 2) / (n - 2)); ax.errorbar(r, y, xerr=1.96 * se, fmt="o", color="#3b6fb6", ms=4, lw=1)
    ax.axvline(0, color="k", lw=0.6, ls=":"); ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=7); ax.set_xlabel("GPX4 dependency vs signature\n(partial r, 95% CI)"); panel(ax, "b")
    ct = pd.read_csv(f"{TBL}/S_ctrp_signature_drug_sensitivity.csv"); c = ct[(ct.signature == "consensus") & (ct.dataset == "CTRPv2")].set_index("compound")
    ax = fig.add_subplot(gs[1, 0]); x = np.arange(len(c)); w = 0.26
    ax.bar(x - w, c.partial_r_lineage, w, color="#3b6fb6", label="All lineages, lineage-adjusted"); ax.bar(x, c.partial_r_lineage_EMT, w, color="#bdbdbd", label="+ EMT-adjusted")
    ax.bar(x + w, c.r_NSCLC, w, color="#d73027", label="NSCLC lines only"); ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(c.index, rotation=90); ax.set_ylabel("Correlation of signature with\ndrug sensitivity (AAC)"); ax.set_ylim(-0.1, 0.48)
    ax.legend(loc="upper right", fontsize=7); panel(ax, "c")
    pr = pd.read_csv(f"{TBL}/S_prism_consensus_signature_all_compounds.csv")
    rc = [k for k in pr.columns if k.lower() in ("partial_r", "r", "partial_r_lineage")][0]
    ax = fig.add_subplot(gs[1, 1]); ax.hist(pr[rc].dropna(), bins=60, color="#bdbdbd")
    ax.set_xlabel("Partial r, signature vs PRISM log fold change\n(all compounds, lineage-adjusted)"); ax.set_ylabel("Compounds"); panel(ax, "d")
    save(fig, "Figure_S3")

figS1(); figS2(); figS3()
