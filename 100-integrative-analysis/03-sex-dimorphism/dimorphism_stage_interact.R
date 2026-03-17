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
                      "dds_sex.rds")

genes_file <- file.path(project_dir,
                        "100-integrative-analysis",
                        "genes.txt")

results_dir <- file.path(project_dir,
                         "100-integrative-analysis",
                         "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

log_info <- function(msg){
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

log_info("Carregando dds_sex...")

if(!file.exists(dds_file)){
  stop("Arquivo dds_sex.rds não encontrado: ", dds_file)
}

dds_sex <- readRDS(dds_file)

# Verificações básicas
if(!"Stage" %in% colnames(colData(dds_sex))){
  stop("Variável 'Stage' não encontrada em colData(dds_sex)")
}

if(!"Sex" %in% colnames(colData(dds_sex))){
  stop("Variável 'Sex' não encontrada em colData(dds_sex)")
}

log_info("Executando modelo ~ Stage * Sex...")

dds_sex <- DESeq(dds_sex)

# Salvar objeto ajustado
saveRDS(dds_sex,
        file.path(project_dir,
                  "100-integrative-analysis",
                  "03-sex-dimorphism",
                  "dds_sex_fitted.rds"))

# =====================================================
# IDENTIFICAR COEFICIENTES
# =====================================================

res_names <- resultsNames(dds_sex)

log_info("Coeficientes disponíveis:")
log_info(paste(res_names, collapse = " | "))

# Coeficiente principal de sexo
sex_coef <- res_names[grep("^Sex_", res_names)][1]

if(is.na(sex_coef)){
  stop("Coeficiente principal de sexo não encontrado.")
}

# Coeficientes de interação Stage:Sex
interaction_coefs <- res_names[grep("Stage.*Sex", res_names)]

if(length(interaction_coefs) == 0){
  log_info("Nenhum coeficiente de interação Stage:Sex encontrado.")
}

# =====================================================
# EXTRAIR RESULTADOS PRINCIPAIS (SEXO)
# =====================================================

log_info(paste("Extraindo coeficiente principal:", sex_coef))

res_sex <- results(dds_sex, name = sex_coef)

res_sex_df <- as.data.frame(res_sex) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    sex_log2FC = as.numeric(log2FoldChange),
    sex_padj = as.numeric(padj),
    sex_significant = case_when(
      is.na(padj) ~ FALSE,
      padj < 0.05 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  select(gene_id,
         baseMean,
         sex_log2FC,
         lfcSE,
         stat,
         pvalue,
         sex_padj,
         sex_significant)

log_info(paste("Genes com dimorfismo significativo:",
               sum(res_sex_df$sex_significant)))

# =====================================================
# EXTRAIR INTERAÇÕES
# =====================================================

interaction_list <- list()

if(length(interaction_coefs) > 0){

  for(coef in interaction_coefs){

    log_info(paste("Extraindo interação:", coef))

    tmp <- results(dds_sex, name = coef) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      mutate(padj = as.numeric(padj)) %>%
      select(gene_id, padj) %>%
      rename(!!paste0("padj_", coef) := padj)

    interaction_list[[coef]] <- tmp
  }

  interaction_df <- reduce(interaction_list,
                           full_join,
                           by = "gene_id")

  final_sex_df <- res_sex_df %>%
    left_join(interaction_df, by = "gene_id")

} else {

  final_sex_df <- res_sex_df
}

# =====================================================
# PROCESSAR GRUPOS
# =====================================================

log_info("Carregando grupos de genes...")

groups <- parse_gene_file(genes_file)

log_info(paste("Total de grupos encontrados:", length(groups)))

for(group in names(groups)){

  log_info(paste("Processando grupo:", group))

  genes <- groups[[group]]
  genes_present <- genes[genes %in% final_sex_df$gene_id]

  if(length(genes_present) == 0){
    log_info(paste("Nenhum gene encontrado para grupo", group))
    next
  }

  group_dir <- file.path(results_dir, group)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

  group_res <- final_sex_df %>%
    filter(gene_id %in% genes_present) %>%
    mutate(across(where(is.numeric), as.numeric))

  readr::write_tsv(group_res,
                   file.path(group_dir, "sex_results.tsv"),
                   na = "")

  log_info(paste("Genes no grupo:", nrow(group_res)))
  log_info(paste("Genes com dimorfismo significativo:",
                 sum(group_res$sex_significant)))
}

log_info("[OK] 03-sex-dimorphism concluído com sucesso.")
