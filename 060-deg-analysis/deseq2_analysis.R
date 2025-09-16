#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
})

project_dir <- "../"
counts_file <- file.path(project_dir, "050-quantification", "counts_matrix.tsv")
meta_file   <- file.path(project_dir, "010-reference", "RNAseq_metadata.tsv")
out_dir     <- file.path(project_dir, "060-deg-analysis")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[INFO] Lendo matrizes e metadata...\n")

counts <- read.delim(counts_file, header = TRUE, sep = "\t", row.names = 1)
counts <- round(counts)
meta   <- read.delim(meta_file, header = TRUE, sep = "\t")

samples <- meta[[3]]
rownames(meta) <- samples
counts <- counts[, samples]  # garantir mesma ordem


cat("[INFO] Criando DESeqDataSet...\n")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Location + Sex)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

cat("[INFO] Exportando resultados DEG...\n")

res <- results(dds, alpha = 0.05)
resOrdered <- res[order(res$padj), ]
res_df <- as.data.frame(resOrdered) %>% rownames_to_column("gene_id")

write_tsv(res_df, file.path(out_dir, "DEGs_results.tsv"))

# =============================
# PCA plot
# =============================
cat("[INFO] Gerando PCA plot...\n")

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("Location", "Sex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = Location, shape = Sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  theme_minimal()

ggsave(file.path(out_dir, "PCA_plot.png"), p_pca, width = 6, height = 5)

# =============================
# Heatmap top genes
# =============================
cat("[INFO] Gerando heatmap...\n")

top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

pheatmap(mat,
         annotation_col = as.data.frame(colData(vsd)[, c("Location", "Sex")]),
         show_rownames = FALSE,
         fontsize_col = 8,
         filename = file.path(out_dir, "heatmap_top50.png"))

# =============================
# Volcano plot
# =============================
cat("[INFO] Gerando volcano plot...\n")

res_df <- res_df %>%
  mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "NotSig"))

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot", x = "log2 Fold Change", y = "-log10 adjusted p-value")

ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 6, height = 5)

cat("[OK] Análise DEG concluída. Resultados em 060-DEG-analysis/\n")
