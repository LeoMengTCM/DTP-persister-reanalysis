#!/usr/bin/env python3
"""Figures 1-4."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np, pandas as pd, scanpy as sc, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

mpl.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 9, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "legend.fontsize": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
                     "legend.frameon": False, "figure.dpi": 100})
W = 7.2
CL_COLORS = {"C1": "#7f7f7f", "C2": "#1b1b1b", "C3": "#d73027", "C4": "#f4a259", "C5": "#3b6fb6"}

_PANELS = []
def panel(ax, letter, dx=None, dy=None):
    _PANELS.append((ax, letter))

def place_panels(fig):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    for ax, letter in _PANELS:
        if ax.figure is not fig: continue
        tb = ax.get_tightbbox(r); x0 = tb.x0 / fig.bbox.width; y1 = max(ax.get_position().y1, tb.y1 / fig.bbox.height)
        fig.text(max(x0 - 0.002, 0.002), min(y1 + 0.002, 0.976), letter, fontsize=12, fontweight="bold", va="bottom", ha="left")
    _PANELS[:] = [(a, l) for a, l in _PANELS if a.figure is not fig]

def check_layout(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fb = fig.bbox
    items = []
    for t in fig.findobj(mpl.text.Text):
        if not t.get_visible() or not t.get_text().strip(): continue
        try: bb = t.get_window_extent(r)
        except Exception: continue
        if bb.width <= 0 or bb.height <= 0: continue
        items.append((t.get_text().replace("\n", " ")[:28], bb, t.get_fontsize()))
    out = [n for n, bb, _ in items if bb.x0 < fb.x0 - 1 or bb.x1 > fb.x1 + 1 or bb.y0 < fb.y0 - 1 or bb.y1 > fb.y1 + 1]
    ov = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            w = min(a.x1, b.x1) - max(a.x0, b.x0); h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if w > 1.5 and h > 1.5 and (w * h) > 0.08 * min(a.width * a.height, b.width * b.height):
                ov.append((items[i][0], items[j][0]))
    axes_bb = [(ax, ax.get_tightbbox(r)) for ax in fig.axes if ax.get_visible()]
    axov = []
    for i in range(len(axes_bb)):
        for j in range(i + 1, len(axes_bb)):
            a, b = axes_bb[i][1], axes_bb[j][1]
            if a is None or b is None: continue
            w = min(a.x1, b.x1) - max(a.x0, b.x0); h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if w > 3 and h > 3: axov.append((i, j, round(w), round(h)))
    print(f"[layout {name}] texts={len(items)} min_font={min(f for *_, f in items):.1f}pt | outside canvas: {out[:8]} | "
          f"text overlaps: {len(ov)} {ov[:10]} | axes overlaps: {axov[:8]}")

def save(fig, name):
    eng = fig.get_layout_engine()
    if eng is not None:
        eng.set(rect=(0, 0, 1, 0.972), h_pad=0.11, w_pad=0.06)
    place_panels(fig)
    check_layout(fig, name)
    fig.savefig(f"{FIG}/{name}.png", dpi=600); fig.savefig(f"{FIG}/{name}.pdf"); plt.close(fig)
    print("saved", name)

def heat(ax, M, vmax, cmap="RdBu_r", annot=None, fs=7, stars=None, xlabels=True, ylabels=True):
    im = ax.imshow(M.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(M.shape[1])); ax.set_yticks(range(M.shape[0]))
    ax.set_xticklabels(M.columns if xlabels else [], rotation=90); ax.set_yticklabels(M.index if ylabels else [])
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            if np.isnan(v):
                ax.text(j, i, "-", ha="center", va="center", fontsize=fs, color="0.5"); continue
            txt = ""
            if annot is not None: txt = annot.format(v)
            if stars is not None and bool(stars.values[i, j]): txt += "*"
            if txt:
                ax.text(j, i, txt, ha="center", va="center", fontsize=fs, color="white" if abs(v) > 0.6 * vmax else "black")
    return im

def figure1():
    a = sc.read_h5ad(f"{CACHE}/primary_atlas.h5ad")
    um = a.obsm["X_umap"]; obs = a.obs
    fig = plt.figure(figsize=(W, 9.2))
    g_top = fig.add_gridspec(1, 3, left=0.06, right=0.985, top=0.955, bottom=0.745, wspace=0.40)
    g_mid = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.585, bottom=0.405, width_ratios=[10, 1.5], wspace=0.02)
    g_bot = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.285, bottom=0.085, wspace=0.28)
    order = np.random.default_rng(0).permutation(a.n_obs)
    ax = fig.add_subplot(g_top[0, 0])
    ax.scatter(um[order, 0], um[order, 1], c=[TP_COLORS[t] for t in obs.timepoint.values[order]], s=2.5, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_title("Treatment day")
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=TP_COLORS[t], label=t, ms=4.5) for t in TP_ORDER], ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), columnspacing=0.8, handletextpad=0.1)
    panel(ax, "a")
    ax = fig.add_subplot(g_top[0, 1])
    ax.scatter(um[order, 0], um[order, 1], c=[CL_COLORS[c] for c in obs.cluster.values[order]], s=2.5, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_title("Cluster")
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=CL_COLORS[c], label=c, ms=4.5) for c in CL_COLORS], ncol=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), columnspacing=0.35, handletextpad=0.0)
    panel(ax, "b")
    ax = fig.add_subplot(g_top[0, 2])
    share = pd.crosstab(obs.timepoint, obs.cluster, normalize="index").reindex(TP_ORDER)
    n = obs.timepoint.value_counts().reindex(TP_ORDER)
    bottom = np.zeros(len(TP_ORDER))
    for c in CL_COLORS:
        ax.bar(range(len(TP_ORDER)), share[c].values, bottom=bottom, color=CL_COLORS[c], width=0.8, label=c); bottom += share[c].values
    ax.set_xticks(range(len(TP_ORDER))); ax.set_xticklabels([f"{t}\n{n[t]}" for t in TP_ORDER], fontsize=7.5)
    ax.set_xlabel("Treatment day (number of cells)"); ax.set_ylabel("Fraction of cells"); ax.set_ylim(0, 1); ax.set_title("Cluster composition")
    panel(ax, "c", dx=-0.12)

    mk = pd.read_csv(f"{TBL}/Table_marker_expression.csv")
    genes = list(dict.fromkeys(mk.gene)); groups = mk.drop_duplicates("gene").set_index("gene").marker_group
    ax = fig.add_subplot(g_mid[0, 0])
    rows = [("cluster", c) for c in ["C1", "C2", "C3", "C4", "C5"]] + [("timepoint", t) for t in TP_ORDER]
    ypos = {}; y = 0
    for i, r in enumerate(rows):
        if i == 5: y += 0.8
        ypos[r] = y; y += 1
    for key in ("cluster", "timepoint"):
        d = mk[mk.grouping == key]
        mean = d.pivot(index="level", columns="gene", values="mean_lognorm")[genes]; pct = d.pivot(index="level", columns="gene", values="pct_expressing")[genes]
        scaled = (mean - mean.min()) / (mean.max() - mean.min())
        for lev in mean.index:
            sca = ax.scatter(range(len(genes)), [ypos[(key, lev)]] * len(genes), s=pct.loc[lev].values * 1.15 + 2, c=scaled.loc[lev].values,
                             cmap="Reds", vmin=0, vmax=1, edgecolors="black", linewidths=0.4)
    ax.set_yticks(list(ypos.values())); ax.set_yticklabels([r[1] for r in rows]); ax.invert_yaxis()
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontstyle="italic")
    ax.set_xlim(-0.7, len(genes) - 0.3); ax.set_ylim(max(ypos.values()) + 0.7, -1.6)
    start = 0
    for g in dict.fromkeys(groups.values):
        k = int((groups.values == g).sum())
        ax.plot([start - 0.35, start + k - 0.65], [-0.85, -0.85], color="black", lw=0.8, clip_on=False)
        ax.text(start + (k - 1) / 2, -1.0, g, ha="center", va="bottom", fontsize=7.5)
        start += k
    lax = fig.add_subplot(g_mid[0, 1]); lax.axis("off")
    hs = [lax.scatter([], [], s=p * 1.15 + 2, c="0.6", edgecolors="black", linewidths=0.4, label=f"{p}%") for p in (20, 60, 100)]
    lax.legend(handles=hs, title="% cells\nexpressing", loc="upper center", bbox_to_anchor=(0.5, 1.0), labelspacing=0.9, title_fontsize=7.5, fontsize=7.5)
    cax = lax.inset_axes([0.38, 0.02, 0.14, 0.42])
    cb = fig.colorbar(sca, cax=cax); cb.set_label("Scaled mean\nexpression", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    panel(ax, "d", dx=-0.03, dy=1.0)

    d = pd.read_csv(f"{TBL}/contrasts/DE_A_PC9_erlotinib_late_vs_naive.csv").set_index("gene")
    ax = fig.add_subplot(g_bot[0, 0])
    y = -np.log10(d.FDR.clip(lower=1e-12)); up = (d.FDR < 0.05) & (d.logFC >= 1); dn = (d.FDR < 0.05) & (d.logFC <= -1)
    ax.scatter(d.logFC[~(up | dn)], y[~(up | dn)], s=2, c="0.75", linewidths=0, rasterized=True)
    ax.scatter(d.logFC[up], y[up], s=5, c="#c0392b", linewidths=0, rasterized=True, label=f"Up (n={int(up.sum())})")
    ax.scatter(d.logFC[dn], y[dn], s=5, c="#2c6fbb", linewidths=0, rasterized=True, label=f"Down (n={int(dn.sum())})")
    ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.6); ax.axvline(1, ls="--", c="0.5", lw=0.6); ax.axvline(-1, ls="--", c="0.5", lw=0.6)
    lab = ["FN1", "SERPINE1", "IGFBP3", "INHBA", "CYP1B1", "ALDH1A3", "TAGLN", "KANK4", "FASN", "S100A6", "HMGA1", "DUSP6", "ETV4", "EMP1", "AGR2", "C3", "HMGA2", "ITGA2"]
    from adjustText import adjust_text
    txts = [ax.text(d.logFC[g], y[g], g, fontsize=7, fontstyle="italic") for g in lab if g in d.index]
    adjust_text(txts, ax=ax, arrowprops=dict(arrowstyle="-", color="0.3", lw=0.4))
    ax.set_xlabel("log$_2$ fold change, D9-D11 vs D0 (pseudobulk, edgeR)"); ax.set_ylabel("-log$_{10}$ FDR")
    ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.0), handletextpad=0.1, columnspacing=1.0, markerscale=1.6)
    panel(ax, "e", dx=-0.1)

    ax = fig.add_subplot(g_bot[0, 1])
    grid = pd.read_csv(f"{TBL}/S_DE_threshold_grid.csv"); dsamp = pd.read_csv(f"{TBL}/S_DE_downsampling.csv"); dm = pd.read_csv(f"{TBL}/S_DE_depth_matched.csv")
    g1 = grid.query("method == 'single-cell Wilcoxon' and abs_logFC == 1 and FDR == 0.05").iloc[0]
    vals = [("Single-cell\nWilcoxon\n(685 vs 1,588\ncells)", g1.up, g1.down),
            ("Single-cell,\ncell-number\nmatched", dsamp.up.median(), dsamp.down.median()),
            ("Single-cell,\ndepth\nmatched", int(((dm.logfoldchanges > 1) & (dm.pvals_adj < 0.05)).sum()), int(((dm.logfoldchanges < -1) & (dm.pvals_adj < 0.05)).sum())),
            ("Pseudobulk\nedgeR\n(4 vs 2\nlibraries)", int(up.sum()), int(dn.sum()))]
    x = np.arange(len(vals)); wd = 0.38
    b1 = ax.bar(x - wd / 2, [v[1] for v in vals], wd, color="#c0392b", label="Up"); b2 = ax.bar(x + wd / 2, [v[2] for v in vals], wd, color="#2c6fbb", label="Down")
    for b in list(b1) + list(b2): ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 8, f"{int(b.get_height())}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x); ax.set_xticklabels([v[0] for v in vals], fontsize=7.5); ax.set_ylabel("Number of genes\n(FDR < 0.05, |log$_2$FC| $\\geq$ 1)", fontsize=8)
    ax.legend(loc="upper right"); ax.set_ylim(0, 880)
    panel(ax, "f", dx=-0.12)
    save(fig, "Figure_1")

def figure2():
    R = pd.read_csv(f"{TBL}/Table_pathway_timepoint_singlecell.csv")
    P = [k for k in METABOLIC_PATHWAYS]
    comps = [f"{t} vs D0" for t in TP_ORDER[1:]]
    D = R.pivot(index="set", columns="comparison", values="cohens_d").loc[P, comps]
    sig = (R.assign(s=(R.MWU_FDR < 0.05) & (R.MWU_FDR_depth_matched < 0.05) & (np.sign(R.cohens_d) == np.sign(R.d_depth_matched)))
           .pivot(index="set", columns="comparison", values="s").loc[P, comps])
    order = D.mean(1).sort_values().index.tolist()
    D, sig = D.loc[order], sig.loc[order]
    lib = pd.read_csv(f"{TBL}/Table_pathway_library_level.csv", index_col=0).loc[order]
    fig = plt.figure(figsize=(W, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[0.78, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    Dd = D.copy(); Dd.index = [PATHWAY_LABELS[i] for i in Dd.index]; Dd.columns = TP_ORDER[1:]
    im = heat(ax, Dd, 1.3, annot="{:+.2f}", fs=6.8, stars=sig)
    ax.set_xticklabels(Dd.columns, rotation=0); ax.set_xlabel("Treatment day (vs D0)"); ax.set_title("Single cells: Cohen's $d$")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02, aspect=22, shrink=0.75, location="bottom"); cb.set_label("Cohen's $d$ vs D0", fontsize=8); cb.ax.tick_params(labelsize=7)
    panel(ax, "a", dx=-0.75, dy=1.01)
    ax = fig.add_subplot(gs[0, 1])
    L = lib.copy(); L.columns = [c.replace("_a", " a").replace("_b", " b").replace("_c", " c") for c in L.columns]
    im = heat(ax, L, 1.5, annot="{:+.1f}", fs=6.3, ylabels=False)
    ax.set_title("Independent libraries: mean log$_2$FC vs D0"); ax.set_xlabel("Drop-seq library (treatment day, replicate)")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02, aspect=26, shrink=0.75, location="bottom"); cb.set_label("Mean log$_2$ fold change vs mean of D0 libraries", fontsize=8); cb.ax.tick_params(labelsize=7)
    panel(ax, "b", dx=-0.02, dy=1.01)

    obs = pd.read_csv(f"{CACHE}/primary_obs_scores.csv", index_col=0); obst = pd.read_csv(f"{CACHE}/primary_obs_scores_depth_matched.csv", index_col=0)
    days = {"D0": 0, "D1": 1, "D2": 2, "D4": 4, "D9": 9, "D11": 11}
    sub = gs[1, :].subgridspec(2, 4, wspace=0.08, hspace=0.25)
    show = ["Fatty_acid_synthesis", "Cholesterol_metabolism", "Glycolysis", "Serine_glycine_1C", "ALDH_family", "Autophagy_core", "Ferroptosis_regulators", "Defense_score_9genes"]
    lab = dict(PATHWAY_LABELS); lab["Defense_score_9genes"] = "Nine-gene defense score"
    for i, k in enumerate(show):
        ax = fig.add_subplot(sub[i // 4, i % 4]); col = f"s_{k}"
        for o, ls, c, name in ((obs, "-", "#222222", "All cells"), (obst, "--", "#c0392b", "Depth-matched")):
            g = o.groupby("timepoint")[col]; m = g.mean().reindex(TP_ORDER); se = (g.std() / np.sqrt(g.size())).reindex(TP_ORDER)
            m0 = m["D0"]; xs = [days[t] for t in TP_ORDER]
            ax.errorbar(xs, m - m0, yerr=1.96 * se, fmt="o" + ls, color=c, ms=3, lw=1.0, capsize=1.5, label=name)
        ax.axhline(0, color="0.6", lw=0.6, ls=":"); ax.set_title(lab[k], fontsize=8.5); ax.set_xticks([0, 1, 2, 4, 9, 11]); ax.tick_params(labelsize=7)
        if i >= 4: ax.set_xlabel("Treatment day")
        if i % 4 == 0: ax.set_ylabel("$\\Delta$ score vs D0")
        if i == 0: ax.legend(fontsize=7, loc="lower right", handlelength=2.2); panel(ax, "c", dx=-0.3, dy=1.05)
    save(fig, "Figure_2")

CONTR = [("A_PC9_erlotinib_D3_vs_naive", "PC9 D3"), ("A_PC9_erlotinib_late_vs_naive", "PC9 D9-11"), ("B_PC9_gef_osi_9d_vs_DMSO", "PC9 9 d"),
         ("C_PC9_acute_vs_DMSO", "PC9"), ("C_H1975_acute_vs_DMSO", "H1975"), ("C_HCC827_acute_vs_DMSO", "HCC827"), ("C_HCC2935_acute_vs_DMSO", "HCC2935"),
         ("C_PC9_DTP_vs_DMSO", "PC9"), ("C_H1975_DTP_vs_DMSO", "H1975"), ("C_HCC827_DTP_vs_DMSO", "HCC827"), ("C_HCC2935_DTP_vs_DMSO", "HCC2935"),
         ("C_PC9_longwash_vs_DMSO", "PC9"), ("C_H1975_longwash_vs_DMSO", "H1975"), ("C_HCC827_longwash_vs_DMSO", "HCC827"), ("C_HCC2935_longwash_vs_DMSO", "HCC2935")]
GROUPS = [(0, 1, "Erlotinib\n(Drop-seq)"), (2, 2, "Gef/Osi\n(sc-rDSeq)"), (3, 6, "Osimertinib 24 h\n(bulk)"), (7, 10, "Osimertinib DTP 21 d\n(bulk)"), (11, 14, "Washout, regrowth\n(bulk)")]

def group_header(ax, n_rows):
    for a0, a1, name in GROUPS:
        ax.plot([a0 - 0.4, a1 + 0.4], [-0.75, -0.75], color="black", lw=0.8, clip_on=False)
        ax.text((a0 + a1) / 2, -0.9, name, ha="center", va="bottom", fontsize=7)
    for a0, a1, _ in GROUPS[1:]:
        ax.axvline(a0 - 0.5, color="white", lw=2.0)

def figure3():
    f = pd.read_csv(f"{TBL}/Table_pathway_fry_all_contrasts.csv")
    cols = [c for c, _ in CONTR]; names = [n for _, n in CONTR]
    P = list(METABOLIC_PATHWAYS)
    M = f.pivot(index="set", columns="contrast", values="mean_logFC").loc[P, cols]; S = (f.pivot(index="set", columns="contrast", values="FDR").loc[P, cols] < 0.05)
    order = M.iloc[:, 7:11].mean(1).sort_values().index; M, S = M.loc[order], S.loc[order]
    M.index = [PATHWAY_LABELS[i] for i in M.index]; M.columns = names
    tf = pd.read_csv(f"{TBL}/Table_TF_activity_CollecTRI_all_contrasts.csv")
    tfs = ["E2F4", "E2F1", "MYC", "SREBF1", "SREBF2", "STAT3", "PPARG", "NFE2L2", "ATF4", "ATF3", "JUN", "NFKB1", "DLX5", "ID4", "TFE3", "TP53", "SMAD3", "FOXO1", "FOXO3"]
    T = tf.pivot(index="source", columns="contrast", values="score").loc[tfs, cols]; TS = (tf.pivot(index="source", columns="contrast", values="FDR").loc[tfs, cols] < 0.05)
    T.columns = names
    g = pd.read_csv(f"{TBL}/Table_target_gene_logFC.csv"); g = g[g.group == "SREBP_targets"]
    genes = ["LDLR", "INSIG1", "HMGCS1", "HMGCR", "MVK", "MVD", "FDPS", "FDFT1", "SQLE", "DHCR7", "DHCR24", "IDI1", "FASN", "ACACA", "ACLY", "SCD", "ELOVL6"]
    G = g.pivot(index="gene", columns="contrast", values="logFC").reindex(index=genes, columns=cols); GS = (g.pivot(index="gene", columns="contrast", values="FDR").reindex(index=genes, columns=cols) < 0.05)
    G.columns = names
    fig = plt.figure(figsize=(W, 9.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[15, 19, 17])
    ax = fig.add_subplot(gs[0]); im = heat(ax, M, 1.2, annot="{:+.1f}", fs=6.3, stars=S, xlabels=False); group_header(ax, len(M))
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, aspect=18); cb.set_label("Mean log$_2$FC\nof pathway genes", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    panel(ax, "a", dx=-0.3, dy=1.06)
    ax = fig.add_subplot(gs[1]); im = heat(ax, T, 10, annot="{:+.0f}", fs=6.3, stars=TS, xlabels=False)
    for a0, a1, _ in GROUPS[1:]: ax.axvline(a0 - 0.5, color="white", lw=2.0)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, aspect=18); cb.set_label("TF activity\n(ULM score)", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    panel(ax, "b", dx=-0.3, dy=1.0)
    ax = fig.add_subplot(gs[2]); im = heat(ax, G, 2.5, annot="{:+.1f}", fs=6.3, stars=GS)
    for a0, a1, _ in GROUPS[1:]: ax.axvline(a0 - 0.5, color="white", lw=2.0)
    ax.set_yticklabels(G.index, fontstyle="italic")
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, aspect=18); cb.set_label("log$_2$FC\nvs control", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    panel(ax, "c", dx=-0.3, dy=1.0)
    save(fig, "Figure_3")

def figure4():
    gw = pd.read_csv(f"{TBL}/S_depmap_genomewide_consensus_lineage_adjusted.csv", index_col=0).sort_values("t")
    res = pd.read_csv(f"{TBL}/Table_depmap_candidate_associations.csv")
    sc_ = pd.read_csv(f"{TBL}/S_depmap_signature_scores.csv", index_col=0)
    hdr = pd.read_csv(f"{DEPMAP}/CRISPRGeneEffect.csv", nrows=0).columns
    gcol = [c for c in hdr if c.split(" (")[0] == "GPX4"][0]
    gp = pd.read_csv(f"{DEPMAP}/CRISPRGeneEffect.csv", usecols=[hdr[0], gcol], index_col=0)[gcol].reindex(sc_.index)
    ok = gp.notna(); lin = sc_.lineage.where(sc_.lineage.map(sc_.lineage.value_counts()) >= 5, "Other_small")
    X = pd.get_dummies(lin[ok], drop_first=True).astype(float); X.insert(0, "c", 1.0); Pm = X.values @ np.linalg.pinv(X.values)
    sr = sc_.consensus[ok].values - Pm @ sc_.consensus[ok].values; yr = gp[ok].values - Pm @ gp[ok].values
    fig = plt.figure(figsize=(W, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(np.arange(1, len(gw) + 1), gw.t.values, s=1.5, c="0.75", linewidths=0, rasterized=True)
    cand = [g for g in CANDIDATES if g in gw.index]; rk = pd.Series(np.arange(1, len(gw) + 1), index=gw.index)
    ax.scatter(rk[cand], gw.t[cand], s=10, c="#2c6fbb", linewidths=0, zorder=3, label="27 candidate metabolic genes")
    ax.scatter(rk["GPX4"], gw.t["GPX4"], s=28, c="#c0392b", zorder=4)
    ax.annotate(f"GPX4 (rank {int(rk['GPX4'])})", (rk["GPX4"], gw.t["GPX4"]), xytext=(2500, gw.t["GPX4"] - 3.2), fontsize=7.5, fontstyle="italic",
                arrowprops=dict(arrowstyle="-", lw=0.5), color="#c0392b")
    top = ["FERMT2", "ITGAV", "TLN1", "PTK2", "YAP1"]
    ax.scatter(rk[top], gw.t[top], s=10, c="#444444", linewidths=0, zorder=3)
    ax.text(0.97, 0.04, "Top-ranked (rank):\n" + "\n".join(f"{g} ({int(rk[g])})" for g in top), transform=ax.transAxes, ha="right", va="bottom", fontsize=7, fontstyle="italic", color="0.25")
    ax.set_xlim(-600, len(gw) + 600); ax.set_ylim(-12.5, 12.5)
    ax.axhline(0, color="0.5", lw=0.6, ls=":"); ax.set_xlabel("Gene rank (18,531 genes)"); ax.set_ylabel("Lineage-adjusted $t$ statistic\n(negative: more essential in signature-high lines)", fontsize=8)
    ax.legend(loc="upper left", handletextpad=0.1, fontsize=7.5)
    panel(ax, "a", dx=-0.16)
    ax = fig.add_subplot(gs[0, 1])
    ns = sc_.NSCLC[ok].values
    ax.scatter(sr[~ns], yr[~ns], s=4, c="0.7", linewidths=0, rasterized=True, label=f"Other lines (n={int((~ns).sum())})")
    ax.scatter(sr[ns], yr[ns], s=9, c="#c0392b", linewidths=0, label=f"NSCLC (n={int(ns.sum())})")
    b = np.polyfit(sr, yr, 1); xx = np.linspace(sr.min(), sr.max(), 50); ax.plot(xx, np.polyval(b, xx), color="black", lw=1)
    r_all = gw.loc["GPX4"]
    ax.set_xlabel("Persister-like signature score\n(lineage-adjusted residual)"); ax.set_ylabel("GPX4 CRISPR gene effect\n(lineage-adjusted residual)")
    ax.text(0.97, 0.97, f"partial $r$ = {r_all.partial_r:.2f}\n$P$ = {r_all.P:.1e}\nn = {int(ok.sum()):,} lines", transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
    ax.legend(loc="lower left", handletextpad=0.1, fontsize=7.5, markerscale=1.3)
    panel(ax, "b", dx=-0.18)

    ana = [("consensus | lineage-adjusted", "Consensus signature,\nlineage-adjusted"), ("consensus | lineage + Hallmark EMT adjusted", "Consensus,\nlineage + EMT-adjusted"),
           ("consensus | NSCLC only", "Consensus,\nNSCLC lines only"), ("pc9_erl | lineage-adjusted", "PC9-erlotinib signature,\nlineage-adjusted"),
           ("orig6 (no candidate genes) | lineage-adjusted", "Literature signature,\nno candidate genes"),
           ("orig9 (original, with GPX4/SLC7A11/GCLC) | lineage-adjusted", "Literature 9-gene\nsignature")]
    M = res.pivot(index="analysis", columns="gene", values="partial_r").loc[[a for a, _ in ana]]
    Q = res.pivot(index="analysis", columns="gene", values="FDR_genomewide").loc[[a for a, _ in ana]] < 0.05
    order = M.iloc[0].sort_values().index; M, Q = M[order], Q[order]; M.index = [n for _, n in ana]
    ax = fig.add_subplot(gs[1, :])
    im = heat(ax, M, 0.25, annot=None, fs=8, stars=Q)
    ax.set_xticklabels(M.columns, rotation=90, fontstyle="italic"); ax.set_yticklabels(M.index, fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01, aspect=14); cb.set_label("Partial $r$ with\ngene effect", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    ax.set_title("Candidate metabolic genes (* genome-wide FDR < 0.05; blue = more essential in signature-high lines)", fontsize=8.5)
    panel(ax, "c", dx=-0.2, dy=1.02)
    save(fig, "Figure_4")

if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4"]
    for w in which:
        {"1": figure1, "2": figure2, "3": figure3, "4": figure4}[w]()
