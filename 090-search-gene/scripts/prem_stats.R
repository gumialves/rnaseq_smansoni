#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript prem_stats.R <arquivo_TPM_long.csv> <pasta_saida>")
}

tpm_file <- args[1]
out_dir  <- args[2]

# ===============================
# 1. Ler dados
# ===============================
data <- read_csv(tpm_file)

# ===============================
# 2. Estatísticas descritivas
# ===============================
stats_desc <- data %>%
  group_by(stage_category) %>%
  summarise(
    n         = n(),
    mean_tpm  = mean(TPM),
    median_tpm= median(TPM),
    sd_tpm    = sd(TPM),
    .groups   = "drop"
  )

# ===============================
# 3. Kruskal-Wallis
# ===============================
kruskal_res <- kruskal_test(data, TPM ~ stage_category)

# ===============================
# 4. Pós-hoc Dunn
# ===============================
dunn_res <- dunn_test(data, TPM ~ stage_category, p.adjust.method = "BH") %>%
  mutate(significance_flag = if_else(p.adj < 0.05, "SIGNIFICATIVO", "ns"))

# ===============================
# 5. Salvar CSV com comparações
# ===============================
comparisons_file <- file.path(out_dir, "comparisons.csv")
write_csv(dunn_res, comparisons_file)

# ===============================
# 6. Salvar relatório TXT
# ===============================
report_file <- file.path(out_dir, "stats_report.txt")
sink(report_file)
cat("Relatório Estatístico -", basename(out_dir), "\n\n")
cat("Arquivo analisado:", tpm_file, "\n\n")

cat("==== Estatísticas descritivas ====\n")
print(stats_desc)

cat("\n==== Teste Kruskal-Wallis ====\n")
print(kruskal_res)

cat("\n==== Testes pós-hoc de Dunn (ajuste BH) ====\n")
print(dunn_res)
sink()

cat("Relatório salvo em:", report_file, "\n")
cat("Comparações salvas em:", comparisons_file, "\n")
