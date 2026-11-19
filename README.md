# DTP re-analysis (EGFR-mutant NSCLC)

Code for the analysis of drug-tolerant persister transcriptomes in EGFR-mutant
NSCLC cell lines and their relation to DepMap dependencies.

## Data

Download into `data/raw/` (not included here):

- GSE149383 (PC9 + erlotinib, Drop-seq) -> `data/raw/scRNA/GSE149383/`
- GSE282074 (PC9 + gefitinib/osimertinib, sc-rDSeq) -> `data/raw/scRNA/GSE282074/`
- GSE193258 (four lines + osimertinib, bulk RNA-seq counts) -> `cache/GSE193258/`
- DepMap 26Q1: `CRISPRGeneEffect.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `Model.csv` -> `data/raw/depmap/`
- PRISM primary repurposing screen -> `cache/prism/`

## Running

Scripts are numbered in the order they are run. Paths are relative to the
repository root; set `DTP_ROOT` to use another location.

    python scripts/01_qc_pseudobulk.py
    python scripts/03_atlas.py
    python scripts/04_gse282074.py
    Rscript scripts/05_contrasts_edgeR.R
    python scripts/02_de_robustness.py
    python scripts/06_tf_hallmark.py
    python scripts/07_depmap.py
    python scripts/08_timecourse_scores.py
    python scripts/09_figures_main.py
    python scripts/10_figures_supp.py

`03_atlas.py` also reads `data/processed/PC9_DTP_timecourse_final.h5ad`, a Harmony-integrated
object used only for the embedding comparison in Fig. S1d.

## Requirements

Python 3.11: scanpy 1.11, decoupler 2.1, anndata, statsmodels, scikit-learn,
adjustText, matplotlib. R 4.6: edgeR 4.10, limma, jsonlite.
