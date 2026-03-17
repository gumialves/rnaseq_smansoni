#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(viridis)
})

project_dir <- "/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"

tpm_file  <- file.path(project_dir,"050-quantification","tpm_matrix.tsv")
meta_file <- file.path(project_dir,"100-integrative-analysis","01-prepare-data","metadata_clean.tsv")
genes_file <- file.path(project_dir,"100-integrative-analysis","genes.txt")
results_dir <- file.path(project_dir,"100-integrative-analysis","results")

# Ordem biológica correta
stage_order <- c("Eggs",
                 "Miracidia",
                 "1d_Sporocysts",
                 "5d_Sporocysts",
                 "32d_Sporocysts",
                 "Cercariae",
		 "2d_Somules",
                 "26d_Juveniles")

# ==========================================================
# Função para ler genes.txt
# ==========================================================

parse_gene_file <- function(file){
  lines <- readLines(file)
  groups <- list()
  for(line in lines){
    if(nchar(line)==0) next
    parts <- strsplit(line, ":")[[1]]
    group_name <- trimws(parts[1])
    genes <- trimws(unlist(strsplit(parts[2], ",")))
    groups[[group_name]] <- genes
  }
  groups
}

# ==========================================================
# Carregar dados base
# ==========================================================

tpm  <- read_tsv(tpm_file, show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

gene_col <- colnames(tpm)[1]

tpm_mat <- tpm %>%
  column_to_rownames(gene_col) %>%
  as.matrix()

tpm_mat <- tpm_mat[, meta$SampleName]
log_tpm <- log2(tpm_mat + 1)

meta$Stage <- factor(meta$Stage, levels = stage_order, ordered = TRUE)
meta <- meta %>% arrange(Stage)

log_tpm <- log_tpm[, meta$SampleName]

annotation_col <- meta %>%
  select(SampleName, Stage, Sex) %>%
  column_to_rownames("SampleName")

groups <- parse_gene_file(genes_file)

# ==========================================================
# Loop por grupo
# ==========================================================

for(group in names(groups)){

  cat("Processando:", group,"\n")

  genes_input   <- groups[[group]]
  genes_present <- genes_input[genes_input %in% rownames(log_tpm)]
  genes_missing <- setdiff(genes_input, genes_present)

  group_dir <- file.path(results_dir, group)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

  summary_file <- file.path(group_dir,"summary.txt")
  sink(summary_file)

  cat("Resumo do Grupo:",group,"\n")
  cat("====================================\n\n")

  cat("Genes fornecidos:",length(genes_input),"\n")
  cat("Genes encontrados:",length(genes_present),"\n")
  cat("Genes ausentes:",length(genes_missing),"\n\n")

  if(length(genes_missing)>0){
    cat("Lista genes ausentes:\n")
    cat(paste(genes_missing,collapse=", "),"\n\n")
  }

  if(length(genes_present)==0){
    cat("Nenhum gene válido encontrado.\n")
    sink()
    next
  }

  # ==========================================================
  # Expressão
  # ==========================================================

  mat_group <- log_tpm[genes_present,]

  # Métricas básicas
  mean_expression <- mean(mat_group)
  mean_sd_gene <- mean(apply(mat_group,1,sd))

  expressed_genes <- sum(apply(mat_group,1,max) > 1)

  cat("Genes com log2(TPM+1) > 1 em pelo menos 1 estágio:",
      expressed_genes,"\n\n")

  cat("Média global log2(TPM+1):",round(mean_expression,3),"\n")
  cat("Desvio padrão médio por gene:",round(mean_sd_gene,3),"\n\n")

  # ==========================================================
  # Heatmap replicatas
  # ==========================================================

  mat_z <- t(scale(t(mat_group)))

  pheatmap(mat_z,
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           color = viridis(100),
           filename = file.path(group_dir,"heatmap_logTPM.png"),
           width=12,height=8)

  # ==========================================================
  # Heatmap média por estágio
  # ==========================================================

  df_long <- as.data.frame(mat_group) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene,
                 names_to="SampleName",
                 values_to="logTPM") %>%
    left_join(meta,by="SampleName")

  stage_means <- df_long %>%
    group_by(gene,Stage) %>%
    summarise(mean_logTPM=mean(logTPM),.groups="drop")

  stage_matrix <- stage_means %>%
    pivot_wider(names_from=Stage,
                values_from=mean_logTPM) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  stage_matrix <- stage_matrix[,stage_order]

  stage_z <- t(scale(t(stage_matrix)))

  pheatmap(stage_z,
           cluster_rows=TRUE,
           cluster_cols=FALSE,
           color=viridis(100),
           filename=file.path(group_dir,"heatmap_stage_means.png"),
           width=10,height=8)

  # ==========================================================
  # Boxplot ordenado biologicamente
  # ==========================================================

  df_long$Stage <- factor(df_long$Stage,
                          levels=stage_order,
                          ordered=TRUE)

  boxplot_file <- file.path(group_dir,"boxplot_stage_means.png")

  p <- ggplot(df_long,
              aes(x=Stage,y=logTPM)) +
    geom_boxplot(fill="steelblue",alpha=0.7) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    labs(title=paste("Expression across development -",group),
         y="log2(TPM+1)")

  ggsave(boxplot_file,p,width=12,height=6)

}

cat("[OK] Visualização concluída.\n")
