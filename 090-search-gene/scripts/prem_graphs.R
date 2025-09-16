#!/usr/bin/env Rscript

# ============================
# Leitura de parâmetros
# ============================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript prem_graphs.R <TPM_file> <gene_target>")
}
tpm_file <- args[1]
gene_target <- args[2]

# ============================
# Pacotes
# ============================
library(tidyverse)
library(pheatmap)

# ============================
# CONFIGURAÇÕES
# ============================
output_dir <- file.path("../results", gene_target)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================
# 1. LEITURA DO TPM MATRIX
# ============================
cat("INFO: Lendo TPM matrix...\n")
tpm_matrix <- read_tsv(tpm_file)

# ============================
# 2. FILTRAR GENE DE INTERESSE E PIVOT LONGER
# ============================
cat("INFO: Filtrando gene de interesse...\n")
gene_data <- tpm_matrix %>%
  filter(gene == gene_target) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "TPM")

# ============================
# 3. EXTRAIR METADADOS (SEX, STAGE, REPLICATA)
# ============================
cat("INFO: Extraindo metadados...\n")
gene_data <- gene_data %>%
  separate(sample, into = c("sex", "stage", "rep"), sep = "_", extra = "merge") %>%
  mutate(
    stage = str_replace_all(stage, "-", "_"),
    rep = str_remove(rep, "^R")
  )

# ============================
# 4. MAPEAR PARA CATEGORIAS DO CICLO DE VIDA
# ============================
cat("INFO: Mapeando estágios para categorias biológicas...\n")
unique_stages <- unique(gene_data$stage)
stage_map <- setNames(unique_stages, unique_stages)

gene_data <- gene_data %>%
  mutate(
    stage_category = recode(stage, !!!stage_map),
    stage_category = factor(stage_category, levels = unique(stage_category))
  )

# ============================
# 5. SALVAR CSV FILTRADO
# ============================
write_csv(gene_data, file.path(output_dir, "TPM_data_long.csv"))

# ============================
# 6. ESTATÍSTICAS POR CATEGORIA
# ============================
stats_cat <- gene_data %>%
  group_by(stage_category) %>%
  summarise(
    mean_tpm = mean(TPM, na.rm = TRUE),
    sd_tpm = sd(TPM, na.rm = TRUE),
    n = n(),
    cv = sd_tpm / mean_tpm,
    .groups = "drop"
  )
write_csv(stats_cat, file.path(output_dir, "stats_by_category.csv"))

kw_test <- kruskal.test(TPM ~ stage_category, data = gene_data)

# ============================
# 7. GRÁFICOS
# ============================
p_box <- ggplot(gene_data, aes(x = stage_category, y = TPM, fill = stage_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Distribuição de TPM para", gene_target),
       x = "Categoria do ciclo de vida", y = "TPM")
ggsave(file.path(output_dir, "boxplot_category.png"), p_box, width = 10, height = 7, dpi = 300)

p_bar <- ggplot(stats_cat, aes(x = stage_category, y = mean_tpm, fill = stage_category)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_tpm - sd_tpm, ymax = mean_tpm + sd_tpm), width = 0.2) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("TPM médio de", gene_target),
       x = "Categoria do ciclo de vida", y = "TPM médio")
ggsave(file.path(output_dir, "barplot_category.png"), p_bar, width = 10, height = 7, dpi = 300)

p_line <- ggplot(gene_data, aes(x = stage_category, y = TPM, group = sex, color = sex)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  theme_bw(base_size = 14) +
  labs(title = paste("Expressão média de", gene_target, "ao longo do ciclo de vida"),
       x = "Categoria do ciclo de vida", y = "TPM", color = "Sexo")
ggsave(file.path(output_dir, "lineplot_category.png"), p_line, width = 10, height = 7, dpi = 300)

p_box_sex <- ggplot(gene_data, aes(x = stage_category, y = TPM, fill = sex)) +
  geom_boxplot() +
  theme_bw(base_size = 14) +
  labs(title = paste("TPM por sexo e categoria -", gene_target),
       x = "Categoria do ciclo de vida", y = "TPM", fill = "Sexo")
ggsave(file.path(output_dir, "boxplot_sex_category.png"), p_box_sex, width = 10, height = 7, dpi = 300)

p_violin <- ggplot(gene_data, aes(x = stage_category, y = TPM, fill = stage_category)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(title = paste("Distribuição de TPM (violino) -", gene_target),
       x = "Categoria do ciclo de vida", y = "TPM")
ggsave(file.path(output_dir, "violin_category.png"), p_violin, width = 10, height = 7, dpi = 300)

heat_data <- gene_data %>%
  unite("sample_full", sex, stage, rep, sep = "_") %>%
  select(sample_full, TPM) %>%
  column_to_rownames("sample_full") %>%
  t()
pheatmap(heat_data, cluster_rows = FALSE, cluster_cols = FALSE,
         main = paste("Heatmap TPM -", gene_target),
         filename = file.path(output_dir, "heatmap_TPM.png"))

# ============================
# 8. RELATÓRIO EM TXT
# ============================
sink(file.path(output_dir, "report.txt"))
cat("Relatório de Expressão -", gene_target, "\n")
cat("Arquivo de entrada:", tpm_file, "\n\n")
cat("Estatísticas por categoria do ciclo de vida:\n")
print(stats_cat)
cat("\nTeste Kruskal-Wallis entre categorias:\n")
print(kw_test)
sink()

cat("Processo concluído! Resultados em:", normalizePath(output_dir), "\n")
