library(limma)
library(ggplot2)
library(ggrepel)

DATA <- "../data/GSE63885"
OUT  <- "../results/limma"
FIG  <- "../figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

exp <- read.csv(file.path(DATA, "expression_for_limma2.csv"),
                header = TRUE, row.names = "Gene.Symbol")
ann <- read.csv(file.path(DATA, "annotation_for_limma.csv"),
                header = TRUE, row.names = "X")

stopifnot(all(rownames(ann) == colnames(exp)))

col_status <- "clinical.status.post.1st.line.chemotherapy..cr...complete.response..pr...partial.response..sd...stable.disease..p...progression..ch1"
slope <- factor(ann[[col_status]], levels = c("pCR", "pNC"), labels = c(1, 0))
pCR <- as.integer(as.vector(slope))

keep <- !is.na(pCR)
exp <- exp[, keep]
ann <- ann[keep, ]
pCR <- pCR[keep]

cat("Samples after filtering: pCR =", sum(pCR == 1), ", pNC =", sum(pCR == 0), "\n")

design <- cbind(npCR = rep(1, length(pCR)), pCR = pCR)
fit <- lmFit(exp, design)
fit <- eBayes(fit)
top <- topTable(fit, coef = "pCR", adjust = "BH", n = Inf)

top$gene_symbol <- rownames(top)
pval_thr <- 0.05
top$significant_logFC1 <- ifelse(top$P.Value < pval_thr & abs(top$logFC) > 1, "Significant", "Not significant")
top$significant_logFC2 <- ifelse(top$P.Value < pval_thr & abs(top$logFC) > 2, "Significant", "Not significant")
top$significant_logFC3 <- ifelse(top$P.Value < pval_thr & abs(top$logFC) > 3, "Significant", "Not significant")

write.csv(top, file.path(OUT, "limma_results.csv"), row.names = FALSE)

counts <- data.frame(
  logFC_threshold = c(1, 2, 3),
  significant_genes = c(
    sum(top$significant_logFC1 == "Significant"),
    sum(top$significant_logFC2 == "Significant"),
    sum(top$significant_logFC3 == "Significant")
  ),
  up_in_pCR = c(
    sum(top$significant_logFC1 == "Significant" & top$logFC > 0),
    sum(top$significant_logFC2 == "Significant" & top$logFC > 0),
    sum(top$significant_logFC3 == "Significant" & top$logFC > 0)
  ),
  up_in_pNC = c(
    sum(top$significant_logFC1 == "Significant" & top$logFC < 0),
    sum(top$significant_logFC2 == "Significant" & top$logFC < 0),
    sum(top$significant_logFC3 == "Significant" & top$logFC < 0)
  )
)
write.csv(counts, file.path(OUT, "limma_significance_counts.csv"), row.names = FALSE)
print(counts)

for (thr in c(1, 2, 3)) {
  col <- paste0("significant_logFC", thr)
  sig <- subset(top, top[[col]] == "Significant")
  top_labels <- head(sig[order(sig$adj.P.Val), ], 15)

  p <- ggplot(top, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = .data[[col]]), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
    geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(pval_thr), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_labels, aes(label = gene_symbol),
                    size = 3, max.overlaps = 25) +
    labs(title = paste0("Volcano plot: pCR vs pNC (|logFC| > ", thr, ", p < 0.05)"),
         x = "log Fold Change",
         y = "-log10(p-value)",
         color = "") +
    theme_minimal()

  ggsave(file.path(FIG, paste0("volcano_limma_logFC", thr, ".png")),
         p, width = 8, height = 6, dpi = 150)
}

cat("\n=== Done ===\n")
