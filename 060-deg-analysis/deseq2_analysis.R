#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  library(rtracklayer)
})

# ================================
# CONFIGURAÇÕES
# ================================
project_dir <- "../"
counts_file <- file.path(project_dir, "050-quantification", "counts_matrix.tsv")
meta_file   <- file.path(project_dir, "010-reference", "RNAseq_metadata.tsv")
gff_file    <- file.path(project_dir, "010-reference/data/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3")
out_dir     <- file.path(project_dir, "060-deg-analysis")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ================================
# LOG
# ================================
log_info <- function(msg) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")

# ================================
# VERIFICAÇÃO DE ARQUIVOS
# ================================
if (!file.exists(counts_file)) stop("Counts matrix não encontrada!")
if (!file.exists(meta_file)) stop("Metadata não encontrada!")

# ================================
# LEITURA DE DADOS
# ================================
log_info("Lendo matrizes e metadata...")
counts <- read.delim(counts_file, header = TRUE, sep = "\t", row.names = 1) %>% round()
meta <- read.delim(meta_file, header = TRUE, sep = "\t")
samples <- meta[[3]]
rownames(meta) <- samples

# Limpeza de fatores
meta <- meta %>%
  mutate(across(c(Location, Sex), ~as.factor(gsub("[^a-zA-Z0-9_.]", "_", .x))))

# Subset counts para amostras presentes no metadata
counts <- counts[, samples]

# Checagem de integridade
zero_genes <- sum(rowSums(counts) == 0)
if(zero_genes > 0) log_info(paste(zero_genes, "genes com contagem zero em todas as amostras"))

# ================================
# FUNÇÃO DE PROCESSAMENTO DO GFF3
# ================================
process_gff3 <- function(gff_file, counts_genes) {
  if(!file.exists(gff_file)) {
    warning("Arquivo GFF3 não encontrado, usando apenas IDs de gene")
    return(data.frame(
      gene_id = counts_genes,
      gene_name = counts_genes,
      biotype = "Unknown",
      stringsAsFactors = FALSE
    ))
  }
  
  log_info("Lendo arquivo de anotações GFF3...")
  gff <- import(gff_file)
  genes <- gff[gff$type == "gene"]
  
  # Garantir que attributes seja character
  attributes_char <- as.character(mcols(genes)$attributes)
  
  get_attr <- function(attr_string, key) {
    sapply(attr_string, function(x) {
      if(is.na(x) || x=="") return(NA)
      kv_pairs <- strsplit(x, ";")[[1]]
      kv_list <- strsplit(kv_pairs, "=")
      kv_list <- kv_list[sapply(kv_list, length)==2]
      val <- kv_list[sapply(kv_list, `[`, 1) == key]
      if(length(val)==0) return(NA)
      kv_list[[which(sapply(kv_list, `[`, 1) == key)]][2]
    }, USE.NAMES = FALSE)
  }
  
  gene_id   <- get_attr(attributes_char, "ID")
  gene_name <- get_attr(attributes_char, "Name")
  biotype   <- get_attr(attributes_char, "biotype")
  
  n <- length(counts_genes)
  
  # Fallback seguro: se não há valores do GFF3, preencher com counts_genes
  if(length(gene_id) == 0)   gene_id   <- counts_genes
  if(length(gene_name) == 0) gene_name <- gene_id
  if(length(biotype) == 0)   biotype   <- rep("Unknown", length(gene_id))
  
  # Garantir que todos têm mesmo tamanho
  gene_id   <- rep_len(gene_id, n)
  gene_name <- rep_len(gene_name, n)
  biotype   <- rep_len(biotype, n)
  
  annotations <- data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    biotype = biotype,
    stringsAsFactors = FALSE
  )
  annotations <- annotations[!duplicated(annotations$gene_id), ]
  return(annotations)
}


annotations <- process_gff3(gff_file, rownames(counts))

# ================================
# ANÁLISE DE EXPRESSÃO DIFERENCIAL
# ================================
log_info("Criando DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
  countData = counts[rowSums(counts) > 10, ],
  colData   = meta,
  design    = ~ Location + Sex
)

dds <- DESeq(dds)
saveRDS(dds, file.path(out_dir, "dds_object.rds"))

res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res[order(res$padj), ]) %>% rownames_to_column("gene_id")

# ================================
# ANOTAÇÃO DE RESULTADOS
# ================================
annotated_res <- res_df %>%
  left_join(annotations, by="gene_id") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
    biotype   = ifelse(is.na(biotype), "Unknown", biotype)
  )

write_tsv(annotated_res, file.path(out_dir, "DEGs_annotated_results.tsv"))

# ================================
# COUNTS NORMALIZADOS
# ================================
norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file.path(out_dir, "normalized_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)

# ================================
# ENRIQUECIMENTO POR BIOTIPO
# ================================
log_info("Realizando análise de enriquecimento por biotipo...")
sig_genes <- annotated_res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_id)
if(length(sig_genes) > 0){
  all_biotypes <- table(annotations$biotype)
  sig_biotypes <- table(annotations$biotype[annotations$gene_id %in% sig_genes])
  enrichment_df <- data.frame(
    Biotype = names(all_biotypes),
    All_Genes = as.numeric(all_biotypes),
    DEGs = as.numeric(sig_biotypes[names(all_biotypes)]),
    stringsAsFactors = FALSE
  )
  enrichment_df$DEGs[is.na(enrichment_df$DEGs)] <- 0
  enrichment_df$Percent_DEGs <- round((enrichment_df$DEGs / enrichment_df$All_Genes)*100, 2)
  enrichment_df$Enrichment_Score <- round((enrichment_df$DEGs/sum(enrichment_df$DEGs)) / 
                                            (enrichment_df$All_Genes/sum(enrichment_df$All_Genes)), 2)
  write_tsv(enrichment_df, file.path(out_dir, "biotype_enrichment.tsv"))
} else {
  log_info("Nenhum gene significativo para enriquecimento")
}

# ================================
# VISUALIZAÇÕES
# ================================
vsd <- vst(dds, blind=FALSE)

# PCA
pcaData <- plotPCA(vsd, intgroup=c("Location","Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_pca <- ggplot(pcaData, aes(PC1, PC2, color=Location, shape=Sex)) +
  geom_point(size=3) + stat_ellipse(level=0.95) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  theme_minimal()
ggsave(file.path(out_dir,"PCA_plot.png"), p_pca, width=8, height=6)

# Heatmap top 100 genes
top_genes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 100)
mat <- t(scale(t(assay(vsd)[top_genes, ])))
pheatmap(mat,
         annotation_col=as.data.frame(colData(vsd)[,c("Location","Sex")]),
         show_rownames=FALSE,
         clustering_method="ward.D2",
         fontsize_col=8,
         filename=file.path(out_dir,"heatmap_top100.png"))

# Volcano plot
volcano_plot <- function(df, title){
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  top_labels <- df %>% filter(padj<0.05 & abs(log2FoldChange)>1) %>% arrange(padj) %>% head(10)
  ggplot(df, aes(x=log2FoldChange, y=-log10(padj), color=(padj<0.05 & abs(log2FoldChange)>1))) +
    geom_point(alpha=0.6) +
    geom_text_repel(data=top_labels, aes(label=gene_id), max.overlaps=20) +
    scale_color_manual(values=c("gray","red")) +
    theme_minimal() +
    labs(title=title, x="log2 Fold Change", y="-log10 padj")
}
p_volcano <- volcano_plot(annotated_res, "Volcano plot - DEG Analysis")
ggsave(file.path(out_dir,"volcano_plot.png"), p_volcano, width=8, height=6)

# MA plot
ma_plot <- function(df, title){
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  ggplot(df, aes(x=baseMean, y=log2FoldChange, color=(padj<0.05 & abs(log2FoldChange)>1))) +
    geom_point(alpha=0.6) + scale_x_log10() +
    scale_color_manual(values=c("gray","red")) +
    theme_minimal() +
    labs(title=title, x="Mean expression", y="log2 Fold Change")
}
p_ma <- ma_plot(annotated_res, "MA plot - DEG Analysis")
ggsave(file.path(out_dir,"MA_plot.png"), p_ma, width=8, height=6)

# ================================
# RELATÓRIO FINAL
# ================================
summary_file <- file.path(out_dir, "analysis_summary.txt")
sink(summary_file)
cat("Análise de Expressão Diferencial - Schistosoma mansoni\n")
cat("=====================================================\n\n")
cat("Data da análise: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Estatísticas gerais:\n")
cat("- Total de genes na matriz: ", nrow(counts), "\n")
cat("- Genes após filtragem (rowSums>10): ", nrow(dds), "\n")
cat("- Genes diferencialmente expressos (padj<0.05 & |log2FC|>1): ", 
    sum(annotated_res$padj<0.05 & abs(annotated_res$log2FoldChange)>1, na.rm=TRUE), "\n\n")
cat("Arquivos gerados:\n")
cat("- Resultados DEG anotados: DEGs_annotated_results.tsv\n")
cat("- Counts normalizados: normalized_counts.tsv\n")
cat("- Objeto DESeq2: dds_object.rds\n")
cat("- Visualizações: PCA_plot.png, heatmap_top100.png, volcano_plot.png, MA_plot.png\n")
sink()
log_info(paste("[OK] Análise DEG concluída. Resultados em", out_dir))
