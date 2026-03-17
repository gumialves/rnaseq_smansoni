#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

# ================================
# CONFIGURAÇÕES
# ================================

project_dir <- "/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"
counts_file <- file.path(project_dir, "050-quantification", "counts_matrix.tsv")
meta_file   <- file.path(project_dir, "010-reference", "RNAseq_metadata.tsv")

out_dir <- file.path(project_dir, "100-integrative-analysis", "01-prepare-data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ================================
# LOG
# ================================

log_info <- function(msg) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")

log_info("Lendo dados...")

counts <- read.delim(counts_file, header = TRUE, sep = "\t", row.names = 1) %>% round()
meta   <- read.delim(meta_file, header = TRUE, sep = "\t")

# ================================
# IDENTIFICAR COLUNA DE SAMPLE
# ================================

# Ajuste aqui se necessário:
sample_col <- colnames(meta)[3]   # <- CONFIRMAR se está correto
meta$SampleName <- meta[[sample_col]]

rownames(meta) <- meta$SampleName

# ================================
# EXTRAIR SEXO, ESTÁGIO, RÉPLICA
# ================================

log_info("Extraindo Sex, Stage e Replicate...")

meta_clean <- meta %>%
  mutate(
    Sex = sub("^([A-Z]+)_.*", "\\1", SampleName),
    Replicate = sub(".*_(R[0-9]+)$", "\\1", SampleName),
    Stage = sub("^[A-Z]+_(.*)_R[0-9]+$", "\\1", SampleName)
  )

# Garantir fatores ordenados
meta_clean$Stage <- factor(meta_clean$Stage)
meta_clean$Sex   <- factor(meta_clean$Sex)

# ================================
# SUBSET COUNTS
# ================================

counts <- counts[, rownames(meta_clean)]

# Filtragem leve (igual ao 060)
counts_filtered <- counts[rowSums(counts) > 10, ]

# ================================
# DATASET A — CICLO DE VIDA
# ================================

log_info("Criando dds_cycle (~ Stage)...")

dds_cycle <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = meta_clean,
  design    = ~ Stage
)

dds_cycle <- estimateSizeFactors(dds_cycle)

saveRDS(dds_cycle, file.path(out_dir, "dds_cycle.rds"))

# ================================
# DATASET B — DIMORFISMO
# ================================

log_info("Criando dds_sex (~ Stage * Sex)...")

# Manter apenas estágios com M e F
stages_with_sex <- meta_clean %>%
  group_by(Stage) %>%
  summarise(n_sex = n_distinct(Sex)) %>%
  filter(all(c("M","F") %in% meta_clean$Sex[meta_clean$Stage == Stage])) %>%
  pull(Stage)

meta_sex <- meta_clean %>%
  filter(Stage %in% stages_with_sex & Sex %in% c("M","F"))

counts_sex <- counts_filtered[, rownames(meta_sex)]

dds_sex <- DESeqDataSetFromMatrix(
  countData = counts_sex,
  colData   = meta_sex,
  design    = ~ Stage * Sex
)

dds_sex <- estimateSizeFactors(dds_sex)

saveRDS(dds_sex, file.path(out_dir, "dds_sex.rds"))

# ================================
# SALVAR METADATA
# ================================

write_tsv(meta_clean, file.path(out_dir, "metadata_clean.tsv"))
write_tsv(as.data.frame(counts_filtered),
          file.path(out_dir, "counts_filtered.tsv"))

log_info("[OK] 01-prepare-data concluído.")
