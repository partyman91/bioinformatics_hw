library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

DATA <- "../data"
OUT  <- "../results/deseq2"
FIG  <- "../figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

expr_raw <- read.csv(file.path(DATA, "raw_counts_ici_samples.tsv"), sep = "\t", row.names = 1)

meta <- read.csv(file.path(DATA, "meta_responses.tsv"), sep = "\t", row.names = 1)
meta <- meta %>% filter(X0 %in% c("R", "NR"))
rownames(meta) <- gsub("-", ".", rownames(meta))

expr_raw <- expr_raw[, rownames(meta)]
expr_raw <- expr_raw[rowSums(expr_raw) >= 10, ]
expr_raw <- round(expr_raw)

dds <- DESeqDataSetFromMatrix(countData = expr_raw, colData = meta, design = ~ X0)
dds$X0 <- relevel(dds$X0, ref = "NR")
dds <- DESeq(dds)

vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "X0", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_data$short_name <- gsub("_.*", "", pca_data$name)
center_x <- median(pca_data$PC1)
center_y <- median(pca_data$PC2)
dist_to_center <- sqrt((pca_data$PC1 - center_x)^2 + (pca_data$PC2 - center_y)^2)
outlier_thr <- quantile(dist_to_center, 0.85)
pca_data$is_outlier <- dist_to_center > outlier_thr

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = X0)) +
  stat_ellipse(aes(group = X0, fill = X0), geom = "polygon",
               alpha = 0.15, level = 0.9, linetype = "dashed") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(data = subset(pca_data, is_outlier),
                  aes(label = short_name),
                  size = 3.2, max.overlaps = 30,
                  box.padding = 0.5, segment.color = "grey50") +
  scale_color_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
  scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA (VST-normalized counts), R vs NR") +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(FIG, "pca.png"), p_pca, width = 8, height = 6, dpi = 150)

write.csv(pca_data, file.path(OUT, "pca_data.csv"), row.names = FALSE)

extreme <- pca_data$PC1 > 80
if (any(extreme)) {
  kept <- rownames(pca_data)[!extreme]
  vsd_clean <- vsd[, kept]
  pca_clean <- plotPCA(vsd_clean, intgroup = "X0", returnData = TRUE)
  var_clean <- round(100 * attr(pca_clean, "percentVar"))
  pca_clean$short_name <- gsub("_.*", "", pca_clean$name)

  p_pca_clean <- ggplot(pca_clean, aes(PC1, PC2, color = X0)) +
    stat_ellipse(aes(group = X0, fill = X0), geom = "polygon",
                 alpha = 0.15, level = 0.9, linetype = "dashed") +
    geom_point(size = 4, alpha = 0.85) +
    scale_color_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    xlab(paste0("PC1: ", var_clean[1], "% variance")) +
    ylab(paste0("PC2: ", var_clean[2], "% variance")) +
    ggtitle("PCA после удаления выброса (PC1 > 80)") +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "pca_no_outlier.png"), p_pca_clean, width = 8, height = 6, dpi = 150)
}

sample_dist <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dist)
annot_col <- data.frame(Response = meta$X0)
rownames(annot_col) <- rownames(meta)

png(file.path(FIG, "clustermap.png"), width = 1200, height = 1000, res = 150)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_col = annot_col,
         main = "Sample-to-sample distances (VST)",
         fontsize = 7)
dev.off()

res <- results(dds, contrast = c("X0", "R", "NR"))
summary(res)

hgnc <- read.csv(file.path(DATA, "hgnc_complete_set.txt"), row.names = 1, sep = "\t")
symbol_map <- setNames(hgnc$symbol, hgnc$ensembl_gene_id)
ens_ids <- rownames(res)
gene_symbols <- ifelse(ens_ids %in% names(symbol_map) & !is.na(symbol_map[ens_ids]),
                       symbol_map[ens_ids], ens_ids)

res_df <- as.data.frame(res)
res_df$ensembl_id <- ens_ids
res_df$gene_symbol <- gene_symbols

pvalue_threshold <- 0.05
log2fc_threshold <- 1

res_df$significant <- ifelse(!is.na(res_df$padj) &
                             res_df$padj < pvalue_threshold &
                             abs(res_df$log2FoldChange) > log2fc_threshold,
                             "Significant", "Not significant")

res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(OUT, "deseq2_results.csv"), row.names = FALSE)

sig_genes <- subset(res_df, significant == "Significant")
write.csv(sig_genes, file.path(OUT, "deseq2_significant.csv"), row.names = FALSE)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = head(sig_genes, 20),
                  aes(label = gene_symbol), size = 3, max.overlaps = 30) +
  labs(title = "Volcano plot: R vs NR (LuC, ICI therapy)",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  xlim(-10, 10) +
  theme_minimal() +
  theme(legend.title = element_blank())
ggsave(file.path(FIG, "volcano_deseq2.png"), p_volcano, width = 8, height = 6, dpi = 150)

png(file.path(FIG, "ma_plot.png"), width = 1000, height = 800, res = 150)
plotMA(res, main = "MA plot: R vs NR", ylim = c(-5, 5))
dev.off()

cat("\n=== Summary ===\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant (padj<0.05, |log2FC|>1):", sum(res_df$significant == "Significant"), "\n")
cat("  up in R:  ", sum(res_df$significant == "Significant" & res_df$log2FoldChange > 0), "\n")
cat("  up in NR: ", sum(res_df$significant == "Significant" & res_df$log2FoldChange < 0), "\n")
