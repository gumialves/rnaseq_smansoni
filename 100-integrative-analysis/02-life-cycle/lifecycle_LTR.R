#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

# =====================================================
# CONFIGURAÇÕES
# =====================================================

project_dir <- "/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"

dds_file <- file.path(project_dir,
                      "100-integrative-analysis",
                      "01-prepare-data",
                      "dds_cycle.rds")

genes_file <- file.path(project_dir,
                        "100-integrative-analysis",
                        "genes.txt")

results_dir <- file.path(project_dir,
                         "100-integrative-analysis",
                         "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

log_info <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n")
}

# =====================================================
# FUNÇÃO PARA LER genes.txt
# =====================================================

parse_gene_file <- function(file){

  if(!file.exists(file)){
    stop("Arquivo genes.txt não encontrado: ", file)
  }

  lines <- readLines(file)

  groups <- list()

  for(line in lines){

    if(trimws(line) == "") next
    if(!grepl(":", line)) next

    parts <- strsplit(line, ":")[[1]]

    group_name <- trimws(parts[1])
    genes <- unlist(strsplit(parts[2], ",")) %>% trimws()

    groups[[group_name]] <- genes
  }

  return(groups)
}

# =====================================================
# CARREGAR DADOS
# =====================================================

log_info("Carregando dds_cycle...")

if(!file.exists(dds_file)){
  stop("Arquivo dds_cycle.rds não encontrado: ", dds_file)
}

dds_cycle <- readRDS(dds_file)

if(!"Stage" %in% colnames(colData(dds_cycle))){
  stop("Variável 'Stage' não encontrada em colData(dds_cycle)")
}

log_info("Executando LRT (~ Stage vs ~1)...")

dds_cycle <- DESeq(dds_cycle,
                   test = "LRT",
                   reduced = ~1)

log_info("Extraindo resultados LRT...")

res_lrt <- results(dds_cycle)

res_df <- as.data.frame(res_lrt) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    padj = as.numeric(padj),
    is_cycle_dynamic = case_when(
      is.na(padj) ~ FALSE,
      padj < 0.05 ~ TRUE,
      TRUE ~ FALSE
    )
  )

log_info(paste("Total de genes testados:", nrow(res_df)))
log_info(paste("Genes dinâmicos (padj < 0.05):",
               sum(res_df$is_cycle_dynamic)))

# =====================================================
# PROCESSAR GRUPOS
# =====================================================

log_info("Carregando grupos de genes...")

groups <- parse_gene_file(genes_file)

log_info(paste("Total de grupos encontrados:", length(groups)))

for(group in names(groups)){

  log_info(paste("Processando grupo:", group))

  genes <- groups[[group]]

  genes_present <- genes[genes %in% res_df$gene_id]

  if(length(genes_present) == 0){
    log_info(paste("Nenhum gene encontrado para grupo", group))
    next
  }

  group_dir <- file.path(results_dir, group)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

  group_res <- res_df %>%
    filter(gene_id %in% genes_present) %>%
    mutate(across(where(is.numeric), as.numeric))

  readr::write_tsv(group_res,
                   file.path(group_dir, "life_cycle.tsv"),
                   na = "")

  log_info(paste("Genes no grupo:", nrow(group_res)))
  log_info(paste("Genes dinâmicos no grupo:",
                 sum(group_res$is_cycle_dynamic)))
}

log_info("[OK] 02-life-cycle concluído com sucesso.")
