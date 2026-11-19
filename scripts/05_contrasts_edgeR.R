#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(edgeR); library(jsonlite)})
rev <- Sys.getenv("DTP_ROOT", ".")
gs <- fromJSON(file.path(rev, "cache/gene_sets.json"))
alias <- unlist(gs$aliases); rev_alias <- setNames(names(alias), alias)
resolve <- function(genes, universe) {
  out <- c()
  for (g in genes) {
    if (g %in% universe) out <- c(out, g)
    else if (!is.na(alias[g]) && alias[g] %in% universe) out <- c(out, alias[g])
    else if (!is.na(rev_alias[g]) && rev_alias[g] %in% universe) out <- c(out, rev_alias[g])
  }
  unique(unname(out))
}
outdir <- file.path(rev, "results/tables/contrasts"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
all_fry <- list(); summ <- list()

fit_model <- function(counts, design) {
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  y <- estimateDisp(y, design, robust = TRUE)
  list(y = y, fit = glmQLFit(y, design, robust = TRUE), design = design)
}
do_contrast <- function(m, contrast, tag, dataset, label) {
  qlf <- glmQLFTest(m$fit, contrast = contrast)
  res <- topTags(qlf, n = Inf, sort.by = "none")$table
  res$gene <- rownames(res); res$t <- sign(res$logFC) * sqrt(res$F)
  write.csv(res[, c("gene", "logFC", "logCPM", "F", "t", "PValue", "FDR")], file.path(outdir, paste0("DE_", tag, ".csv")), row.names = FALSE)
  idx <- lapply(gs$sets, function(g) which(rownames(m$y) %in% resolve(g, rownames(m$y))))
  fr <- fry(m$y, index = idx, design = m$design, contrast = contrast, sort = "none")
  fr$set <- rownames(fr)
  fr$mean_logFC <- sapply(idx, function(i) mean(res$logFC[i]))
  fr$median_logFC <- sapply(idx, function(i) median(res$logFC[i]))
  fr$contrast <- tag; fr$dataset <- dataset; fr$label <- label
  all_fry[[tag]] <<- fr
  summ[[tag]] <<- data.frame(contrast = tag, dataset = dataset, label = label, n_genes = nrow(res),
                             BCV = round(sqrt(m$y$common.dispersion), 3),
                             up_FDR05_lfc1 = sum(res$FDR < 0.05 & res$logFC > 1), down_FDR05_lfc1 = sum(res$FDR < 0.05 & res$logFC < -1),
                             up_FDR05 = sum(res$FDR < 0.05 & res$logFC > 0), down_FDR05 = sum(res$FDR < 0.05 & res$logFC < 0))
  cat(tag, "done\n")
}

cnt <- as.matrix(read.csv(file.path(rev, "cache/pseudobulk_counts.csv"), row.names = 1, check.names = FALSE))
sA <- c("D0_a", "D0_b", "D3_a", "D3_b", "D9_a", "D11_a", "D11_b", "D11_c")
gA <- factor(c("naive", "naive", "D3", "D3", "late", "late", "late", "late"), levels = c("naive", "D3", "late"))
dA <- model.matrix(~0 + gA); colnames(dA) <- levels(gA)
mA <- fit_model(cnt[, sA], dA)
do_contrast(mA, makeContrasts(late - naive, levels = dA), "A_PC9_erlotinib_late_vs_naive", "GSE149383", "PC9 erlotinib D9-11 vs D0")
do_contrast(mA, makeContrasts(D3 - naive, levels = dA), "A_PC9_erlotinib_D3_vs_naive", "GSE149383", "PC9 erlotinib D3 vs D0")
do_contrast(mA, makeContrasts(late - D3, levels = dA), "A_PC9_erlotinib_late_vs_D3", "GSE149383", "PC9 erlotinib D9-11 vs D3")

cb <- as.matrix(read.csv(file.path(rev, "cache/GSE282074_pseudobulk_counts.csv"), row.names = 1, check.names = FALSE))
gB <- factor(c("DMSO", "DMSO", "TKI", "TKI", "TKI"), levels = c("DMSO", "TKI"))
dB <- model.matrix(~0 + gB); colnames(dB) <- levels(gB)
mB <- fit_model(cb, dB)
do_contrast(mB, makeContrasts(TKI - DMSO, levels = dB), "B_PC9_gef_osi_9d_vs_DMSO", "GSE282074", "PC9 gefitinib/osimertinib 9 d vs DMSO")

cc <- read.delim(gzfile(file.path(rev, "cache/GSE193258/GSE193258_RNAseq_estimated_counts.tsv.gz")), row.names = 1, check.names = FALSE)
cc <- round(as.matrix(cc))
sm <- data.frame(sample = colnames(cc))
sm$line <- sub("_.*", "", sm$sample)
sm$cond <- sub("_[0-9]+$", "", sub("^[^_]+_", "", sm$sample))
sm$cond <- gsub("osi_", "", sm$cond)
print(table(sm$line, sm$cond))
write.csv(sm, file.path(rev, "cache/GSE193258/sample_sheet.csv"), row.names = FALSE)
for (ln in unique(sm$line)) {
  i <- sm$line == ln
  g <- factor(sm$cond[i], levels = c("DMSO", "acute", "DTP", "short_wash", "long_wash"))
  d <- model.matrix(~0 + g); colnames(d) <- levels(g)
  m <- fit_model(cc[, i], d)
  do_contrast(m, makeContrasts(DTP - DMSO, levels = d), paste0("C_", ln, "_DTP_vs_DMSO"), "GSE193258", paste(ln, "osimertinib DTP 21 d vs DMSO"))
  do_contrast(m, makeContrasts(acute - DMSO, levels = d), paste0("C_", ln, "_acute_vs_DMSO"), "GSE193258", paste(ln, "osimertinib 24 h vs DMSO"))
  do_contrast(m, makeContrasts(DTP - acute, levels = d), paste0("C_", ln, "_DTP_vs_acute"), "GSE193258", paste(ln, "osimertinib DTP 21 d vs 24 h"))
  do_contrast(m, makeContrasts(long_wash - DMSO, levels = d), paste0("C_", ln, "_longwash_vs_DMSO"), "GSE193258", paste(ln, "drug washout (regrowth) vs DMSO"))
  lc <- cpm(m$y, log = TRUE, prior.count = 2)
  write.csv(lc, file.path(rev, sprintf("cache/GSE193258/logCPM_%s.csv", ln)))
}
g <- factor(sm$cond, levels = c("DMSO", "acute", "DTP", "short_wash", "long_wash")); ln <- factor(sm$line)
d <- model.matrix(~0 + g + ln); colnames(d) <- sub("^g", "", colnames(d))
m <- fit_model(cc, d)
do_contrast(m, makeContrasts(DTP - DMSO, levels = d), "C_pooled4lines_DTP_vs_DMSO", "GSE193258", "4 EGFR-mutant lines pooled, DTP vs DMSO")
do_contrast(m, makeContrasts(acute - DMSO, levels = d), "C_pooled4lines_acute_vs_DMSO", "GSE193258", "4 EGFR-mutant lines pooled, 24 h vs DMSO")
do_contrast(m, makeContrasts(long_wash - DMSO, levels = d), "C_pooled4lines_longwash_vs_DMSO", "GSE193258", "4 EGFR-mutant lines pooled, washout vs DMSO")

fr <- do.call(rbind, all_fry)
write.csv(fr, file.path(rev, "results/tables/Table_pathway_fry_all_contrasts.csv"), row.names = FALSE)
sm2 <- do.call(rbind, summ)
write.csv(sm2, file.path(rev, "results/tables/Table_contrast_summary.csv"), row.names = FALSE)
print(sm2, row.names = FALSE)
yA <- normLibSizes(DGEList(cnt[rowSums(cnt > 0) >= 6, ]))
write.csv(cpm(yA, log = TRUE, prior.count = 2), file.path(rev, "cache/pseudobulk_logCPM_all_libraries.csv"))
