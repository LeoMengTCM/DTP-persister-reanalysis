"""Paths and gene sets."""
import os

PROJECT = os.environ.get("DTP_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REV = PROJECT
RAW149 = f"{PROJECT}/data/raw/scRNA/GSE149383"
RAW282 = f"{PROJECT}/data/raw/scRNA/GSE282074"
PROC = f"{PROJECT}/data/processed"
DEPMAP = f"{PROJECT}/data/raw/depmap"
CACHE = f"{REV}/cache"
TBL = f"{REV}/results/tables"
FIG = f"{REV}/results/figures"
for d in (CACHE, TBL, FIG):
    os.makedirs(d, exist_ok=True)

TP_ORDER = ["D0", "D1", "D2", "D4", "D9", "D11"]
TP_COLORS = {"D0": "#4d4d4d", "D1": "#d73027", "D2": "#f46d43", "D4": "#fdae61",
             "D9": "#74add1", "D11": "#4575b4"}

METABOLIC_PATHWAYS = {
    "Glycolysis": ["HK1", "HK2", "GPI", "PFKP", "PFKL", "PFKM", "ALDOA", "ALDOB",
                   "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB"],
    "TCA_cycle": ["CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "OGDH", "SUCLA2",
                  "SUCLG1", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2"],
    "OXPHOS_CI": ["NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFB1", "NDUFB2",
                  "NDUFB3", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS6", "NDUFS7", "NDUFS8"],
    "OXPHOS_CII_V": ["SDHA", "SDHB", "UQCRC1", "UQCRC2", "UQCRFS1", "CYC1", "CYCS",
                     "COX4I1", "COX5A", "COX5B", "COX6A1", "COX7A2", "ATP5F1A", "ATP5F1B", "ATP5F1C"],
    "Fatty_acid_oxidation": ["CPT1A", "CPT1B", "CPT1C", "CPT2", "ACOX1", "ACOX2", "ACADM",
                             "ACADL", "ACADVL", "HADHA", "HADHB", "ECHS1", "ECI1", "DECR1"],
    "Fatty_acid_synthesis": ["FASN", "ACACA", "ACACB", "ACLY", "SCD", "FADS1", "FADS2",
                             "ELOVL1", "ELOVL5", "ELOVL6", "MCAT", "ACSL1", "ACSL3", "ACSL4", "ACSL5"],
    "Cholesterol_metabolism": ["HMGCR", "HMGCS1", "HMGCS2", "MVK", "MVD", "FDPS", "SQLE",
                               "DHCR7", "DHCR24", "LDLR", "ABCA1", "NPC1", "NPC2", "SREBF1", "SREBF2"],
    "Glutamine_metabolism": ["GLS", "GLS2", "GLUD1", "GLUD2", "GLUL", "SLC1A5", "SLC7A5",
                             "SLC38A1", "SLC38A2", "GOT1", "GOT2", "GPT", "GPT2", "ASNS"],
    "Serine_glycine_1C": ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD1", "MTHFD2",
                          "MTHFD1L", "ALDH1L1", "ALDH1L2", "TYMS", "DHFR", "MTR"],
    "Tryptophan_IDO": ["IDO1", "IDO2", "TDO2", "KMO", "KYNU", "HAAO", "QPRT", "AFMID",
                       "WARS1", "SLC7A5", "AHR"],
    "Glutathione_redox": ["GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2", "GPX3", "GPX4",
                          "GSTP1", "GSTA1", "MGST1", "MGST2", "TXNRD1", "TXN", "PRDX1"],
    "Ferroptosis_regulators": ["GPX4", "SLC7A11", "SLC3A2", "NFE2L2", "KEAP1", "ACSL4", "LPCAT3",
                               "FTH1", "FTL", "TFRC", "IREB2", "HMOX1", "NQO1", "DHODH", "FSP1"],
    "ROS_detox": ["SOD1", "SOD2", "SOD3", "CAT", "PRDX1", "PRDX2", "PRDX3", "PRDX4",
                  "PRDX5", "PRDX6", "TXN", "TXN2", "TXNRD1", "TXNRD2", "GLRX"],
    "ALDH_family": ["ALDH1A1", "ALDH1A2", "ALDH1A3", "ALDH1B1", "ALDH2", "ALDH3A1",
                    "ALDH3A2", "ALDH3B1", "ALDH3B2", "ALDH4A1", "ALDH5A1", "ALDH6A1",
                    "ALDH7A1", "ALDH8A1", "ALDH9A1"],
    "Autophagy_core": ["ATG5", "ATG7", "ATG12", "BECN1", "MAP1LC3B", "SQSTM1", "ULK1",
                       "ULK2", "ATG13", "ATG14", "PIK3C3", "LAMP1", "LAMP2", "BNIP3", "BNIP3L"],
}

ALIASES = {"FSP1": "AIFM2", "WARS1": "WARS", "ATP5F1A": "ATP5A1", "ATP5F1B": "ATP5B", "ATP5F1C": "ATP5C1"}

PATHWAY_LABELS = {
    "Glycolysis": "Glycolysis", "TCA_cycle": "TCA cycle", "OXPHOS_CI": "OXPHOS complex I",
    "OXPHOS_CII_V": "OXPHOS complexes II-V", "Fatty_acid_oxidation": "Fatty acid oxidation",
    "Fatty_acid_synthesis": "Fatty acid synthesis", "Cholesterol_metabolism": "Cholesterol metabolism",
    "Glutamine_metabolism": "Glutamine metabolism", "Serine_glycine_1C": "Serine-glycine-1C",
    "Tryptophan_IDO": "Tryptophan/kynurenine", "Glutathione_redox": "Glutathione redox",
    "Ferroptosis_regulators": "Ferroptosis regulators", "ROS_detox": "ROS detoxification",
    "ALDH_family": "ALDH family", "Autophagy_core": "Autophagy core",
}

DEFENSE_GENES = ["GPX4", "SLC7A11", "GCLC", "GCLM", "ALDH1A1", "ALDH3A1", "CPT1A", "SOD2", "NFE2L2"]

ORIG_DTP_SIGNATURE = ["AXL", "IGFBP3", "VIM", "FN1", "SERPINE1", "ALDH1A1", "GPX4", "SLC7A11", "GCLC"]

CANDIDATES = ["GPX4", "SLC7A11", "GCLC", "GCLM", "GSS", "GSR",
              "ALDH1A1", "ALDH1A3", "ALDH3A1", "ALDH2",
              "CPT1A", "CPT2", "ACOX1", "HADHA", "HADHB",
              "PHGDH", "PSAT1", "MTHFD2", "SHMT2", "DHFR",
              "SOD2", "NQO1", "TXNRD1", "NFE2L2",
              "ATG5", "ATG7", "BECN1"]

def all_metabolic_genes():
    s = set(CANDIDATES) | set(DEFENSE_GENES)
    for g in METABOLIC_PATHWAYS.values():
        s |= set(g)
    s |= set(ALIASES.values())
    return s

def resolve(genes, var_names):
    vs = set(var_names)
    rev = {v: k for k, v in ALIASES.items()}
    out = []
    for g in genes:
        if g in vs:
            out.append(g)
        elif g in ALIASES and ALIASES[g] in vs:
            out.append(ALIASES[g])
        elif g in rev and rev[g] in vs:
            out.append(rev[g])
    return list(dict.fromkeys(out))

def cohens_d(a, b):
    import numpy as np
    a = np.asarray(a, float); b = np.asarray(b, float)
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / s if s > 0 else np.nan

def set_pub_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11,
        "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
        "savefig.dpi": 300, "figure.dpi": 100, "axes.linewidth": 1.0,
    })
