#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(gridExtra)
})

# =====================================================
# CONFIG
# =====================================================

project_dir <- "/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"

dds_cycle_file <- file.path(project_dir,
                            "100-integrative-analysis",
                            "01-prepare-data",
                            "dds_cycle.rds")

results_dir <- file.path(project_dir,
                         "100-integrative-analysis",
                         "results")

log_info <- function(msg){
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n")
}

# =====================================================
# FUNÇÃO TAU
# =====================================================

calc_tau <- function(x){

  if(all(is.na(x))) return(NA)

  x <- as.numeric(x)

  if(max(x, na.rm=TRUE) == 0) return(0)

  x_norm <- x / max(x, na.rm=TRUE)

  tau <- sum(1 - x_norm, na.rm=TRUE) / (length(x) - 1)

  return(tau)
}

# =====================================================
# CARREGAR MATRIZ NORMALIZADA
# =====================================================

log_info("Carregando dds_cycle para cálculo de Tau...")

dds_cycle <- readRDS(dds_cycle_file)

norm_counts <- counts(dds_cycle, normalized=TRUE)

stage <- colData(dds_cycle)$Stage

stage_levels <- unique(stage)

# média por estágio
stage_matrix <- sapply(stage_levels, function(s){
  rowMeans(norm_counts[, stage == s, drop=FALSE])
})

stage_matrix <- as.data.frame(stage_matrix)
stage_matrix$gene_id <- rownames(stage_matrix)

# =====================================================
# CALCULAR TAU
# =====================================================

log_info("Calculando Tau...")

stage_matrix$tau <- apply(stage_matrix[, stage_levels],
                          1,
                          calc_tau)

tau_df <- stage_matrix %>%
  select(gene_id, tau)

# =====================================================
# PROCESSAR GRUPOS
# =====================================================

groups <- list.dirs(results_dir,
                    full.names = FALSE,
                    recursive = FALSE)

summary_global <- list()

for(group in groups){

  log_info(paste("Integrando grupo:", group))

  group_dir <- file.path(results_dir, group)

  life_file <- file.path(group_dir, "life_cycle.tsv")
  sex_file  <- file.path(group_dir, "sex_results.tsv")

  if(!file.exists(life_file) | !file.exists(sex_file)){
    next
  }

  life_df <- read_tsv(life_file, show_col_types=FALSE)
  sex_df  <- read_tsv(sex_file, show_col_types=FALSE)

  df <- full_join(life_df, sex_df, by="gene_id") %>%
        left_join(tau_df, by="gene_id")

  safe_log10 <- function(x){
    x <- ifelse(is.na(x), 1, x)
    -log10(pmax(x, 1e-300))
  }

  df <- df %>%
    mutate(
      cycle_score = safe_log10(padj),
      sex_score   = safe_log10(sex_padj),
      integrated_score = cycle_score + sex_score
    ) %>%
    arrange(desc(integrated_score))

  # salvar tabela integrada
  write_tsv(df,
            file.path(group_dir,
                      "integrated_results.tsv"))

  # summary (movido do 04)
  summary_df <- df %>%
    summarise(
      total_genes = n(),
      mean_tau = mean(tau, na.rm=TRUE),
      high_tau = sum(tau > 0.8, na.rm=TRUE),
      cycle_sig = sum(cycle_score > 2),
      sex_sig = sum(sex_score > 2)
    )

  write_tsv(summary_df,
            file.path(group_dir,
                      "group_summary.tsv"))

  summary_df$group <- group
  summary_global[[group]] <- summary_df

  # =================================================
  # =================================================

pdf(file.path(group_dir, "Figure_integrative_paper_style.pdf"),
    width = 12, height = 10)

# A - Tau distribution
p1 <- ggplot(df, aes(tau)) +
      geom_histogram(bins = 30) +
      theme_bw() +
      ggtitle("A) Tau distribution")
print(p1)

# B - Integrated score
p2 <- ggplot(df, aes(integrated_score)) +
      geom_histogram(bins = 30) +
      theme_bw() +
      ggtitle("B) Integrated score")
print(p2)

# C - Scatter
p3 <- ggplot(df, aes(tau, integrated_score)) +
      geom_point(alpha = 0.4) +
      theme_bw() +
      ggtitle("C) Tau vs Integrated")
print(p3)

# D - Heatmap top 30
top_genes <- df %>%
  slice_max(integrated_score, n = 30) %>%
  pull(gene_id)

heat_df <- stage_matrix[stage_matrix$gene_id %in% top_genes, ]
heat_df <- heat_df[match(top_genes, heat_df$gene_id), ]
heat_mat <- as.matrix(heat_df[, stage_levels])
rownames(heat_mat) <- heat_df$gene_id

pheatmap(heat_mat,
         scale = "row",
         main = "D) Top 30 stage expression")

dev.off()
}
# salvar summary global
bind_rows(summary_global) %>%
  write_tsv(file.path(results_dir,
                      "GLOBAL_group_summary.tsv"))

log_info("[OK] 05 Integrative analysis finalizado.")
