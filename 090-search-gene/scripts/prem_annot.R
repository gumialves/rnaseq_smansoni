#!/usr/bin/env Rscript

# ============================
# Leitura de parâmetros
# ============================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript pre_annot.R <TPM_file> <gene_target>")
}
tpm_file <- args[1]
gene_target <- args[2]

# ============================
# CONFIGURAÇÕES
# ============================
library(tidyverse)
library(seqinr)
library(rtracklayer)

gff_file <- "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3" 
protein_fasta <- "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa"
cds_fasta <- "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.CDS_transcripts.fa"

output_base <- file.path("../results", gsub("[:.]", "_", gene_target), "/")

# ============================
# 1. LEITURA DOS DADOS DE EXPRESSÃO (CORRIGIDA)
# ============================
cat("INFO: Lendo TPM matrix...\n")
tpm_matrix <- read_tsv(tpm_file, show_col_types = FALSE)

# Verificar nomes de genes na matriz
cat("INFO: Primeiros genes na matriz TPM:\n")
print(head(unique(tpm_matrix$gene)))

# Criar padrão flexível para busca do gene
gene_pattern <- paste0("(", gsub("[:.]", "[:.]?", gene_target), ")|(", 
                       gsub("^gene:", "", gene_target), ")")

# Filtrar gene de interesse com expressão regular
gene_data <- tpm_matrix %>%
  filter(str_detect(gene, regex(gene_pattern, ignore_case = TRUE)))

# Verificar se encontrou dados
if (nrow(gene_data) == 0) {
  cat("ERRO: Gene não encontrado na matriz TPM. Padrão usado:", gene_pattern, "\n")
  cat("INFO: Exemplos de genes presentes:\n")
  print(sample(unique(tpm_matrix$gene), 10))
  stop("Abortando execução: gene não encontrado na matriz de expressão")
} else {
  cat("INFO: Encontradas", nrow(gene_data), "amostras para o gene\n")
  
  # Atualizar gene_target para o nome real encontrado
  actual_gene_name <- unique(gene_data$gene)
  if (length(actual_gene_name) == 1) {
    gene_target <- actual_gene_name
    cat("INFO: Usando nome do gene:", gene_target, "\n")
  }
}

# Processamento dos dados
gene_data <- gene_data %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  separate(sample, into = c("sex", "stage", "rep"), sep = "_", extra = "merge") %>%
  mutate(
    stage = str_replace(stage, "-", "_"),
    rep = str_remove(rep, "^R")
  )

# Mapeamento de estágios
stage_map <- c(
  "Eggs" = "Ovos",
  "Miracidia" = "Miracídio",
  "Sporocysts_1d" = "Esporocisto_1d",
  "Sporocysts_5d" = "Esporocisto_5d",
  "Sporocysts_32d" = "Esporocisto_32d",
  "Cercariae" = "Cercária",
  "Somules_2d" = "Esquistossômulo_2d",
  "Juveniles" = "Juvenis"
)

gene_data <- gene_data %>%
  mutate(
    stage_category = recode(stage, !!!stage_map),
    stage_category = factor(stage_category,
                           levels = c("Ovos", "Miracídio", "Esporocisto_1d",
                                      "Esporocisto_5d", "Esporocisto_32d", 
                                      "Cercária", "Esquistossômulo_2d", 
                                      "Juvenis"))
  )

# Salvar dados de expressão
write_csv(gene_data, paste0(output_base, "_TPM_data.csv"))

# ============================
# 2. ANÁLISE DE ANOTAÇÃO (VERSÃO CORRIGIDA)
# ============================
cat("\nINFO: Analisando anotação genética...\n")

extract_gene_info <- function(gff, gene_id) {
  # Converter para dataframe para manipulação segura
  gff_df <- as.data.frame(gff)
  
  # Encontrar gene com tolerância a formatos de ID
  patterns <- c(gene_id, 
                paste0("gene:", gene_id),
                gsub("^gene:", "", gene_id))
  
  gene_entry <- gff_df %>%
    filter(type == "gene") %>%
    filter(ID %in% patterns)
  
  if (nrow(gene_entry) == 0) return(NULL)
  
  # Extrair metadados de forma segura
  get_metadata <- function(field, default = NA) {
    if (field %in% colnames(gene_entry)) {
      value <- gene_entry[[field]][1]
      ifelse(is.na(value), default, as.character(value))
    } else {
      default
    }
  }
  
  # Buscar transcritos associados (CORREÇÃO APLICADA)
  mrnas <- gff_df %>%
    filter(type == "mRNA") %>%
    mutate(Parent = as.character(Parent)) %>%  # Forçar conversão para caractere
    filter(str_detect(Parent, fixed(gene_entry$ID[1])))
  
  # Buscar características estruturais
  features <- gff_df %>%
    filter(Parent %in% mrnas$ID)
  
  exons <- features %>% filter(type == "exon")
  cds <- features %>% filter(type == "CDS")
  
  # Calcular comprimentos
  gene_length <- gene_entry$width
  cds_length <- sum(cds$width, na.rm = TRUE)
  
  # Termos GO - abordagem mais robusta
  go_terms <- character(0)
  ontology <- get_metadata("Ontology_term")
  if (!is.na(ontology)) {
    go_terms <- unlist(strsplit(ontology, ",")) %>%
      str_trim() %>%
      keep(~ str_detect(.x, "GO:"))
  }
  
  # Retornar estrutura
  list(
    gene_id = gene_entry$ID[1],
    location = paste0(gene_entry$seqnames[1], ":", gene_entry$start[1], "-", gene_entry$end[1]),
    strand = as.character(gene_entry$strand[1]),
    gene_length = gene_length,
    cds_length = cds_length,
    exon_count = nrow(exons),
    exon_lengths = paste(exons$width, collapse = ", "),
    go_terms = go_terms,
    description = get_metadata("description", "Sem descrição")
  )
}

# ============================
# 3. ANÁLISE DE SEQUÊNCIAS (VERSÃO CORRIGIDA)
# ============================
find_sequence <- function(fasta_file, pattern) {
  seqs <- tryCatch(
    {
      read.fasta(fasta_file, seqtype = "AA", as.string = TRUE)
    },
    error = function(e) {
      cat("ERRO na leitura do arquivo FASTA:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(seqs)) return(NA)
  
  seq_names <- names(seqs)
  
  # Padrões de busca ampliados
  patterns <- c(
    pattern,
    gsub("^gene:", "", pattern),
    paste0(pattern, "\\.\\d+"),
    paste0(gsub("^gene:", "", pattern), "\\.\\d+"),
    gsub("_", "-", pattern),
    gsub(":", "_", pattern)
  )
  
  for (pat in patterns) {
    match <- grep(pat, seq_names, value = TRUE)
    if (length(match) > 0) {
      cat("INFO: Encontrada sequência para padrão '", pat, "'\n")
      return(seqs[[match[1]]][1])
    }
  }
  
  cat("WARN: Sequência não encontrada para", pattern, "em", basename(fasta_file), "\n")
  return(NA)
}

# ============================
# FLUXO PRINCIPAL CORRIGIDO
# ============================

# Leitura do GFF com tratamento de erros
gff <- tryCatch(
  {
    import.gff3(gff_file)
  },
  error = function(e) {
    cat("ERRO na leitura do GFF3, tentando formato genérico...\n")
    import(gff_file, format = "gff3")
  }
)

# Extrair informações do gene
gene_info <- extract_gene_info(gff, gene_target)

# Tentativas alternativas se gene não for encontrado
if (is.null(gene_info)) {
  cat("WARN: Gene não encontrado com ID original, tentando alternativas...\n")
  alt_ids <- c(
    gsub("^gene:", "", gene_target),
    gsub("_", "-", gene_target),
    gsub(":", "_", gene_target),
    paste0("gene:", gene_target)
  )
  
  for (alt_id in alt_ids) {
    gene_info <- extract_gene_info(gff, alt_id)
    if (!is.null(gene_info)) {
      cat("INFO: Gene encontrado com ID alternativo:", alt_id, "\n")
      gene_target <- alt_id
      break
    }
  }
  
  if (is.null(gene_info)) {
    cat("WARN: Gene ainda não encontrado após tentativas alternativas\n")
    gene_info <- list(gene_id = gene_target, description = "Gene não encontrado na anotação")
  }
}

# Busca de sequências
if (!is.null(gene_info$gene_id)) {
  cat("INFO: Buscando sequência proteica para", gene_info$gene_id, "\n")
  protein_seq <- find_sequence(protein_fasta, gene_info$gene_id)
  
  cat("INFO: Buscando sequência CDS para", gene_info$gene_id, "\n")
  cds_seq <- find_sequence(cds_fasta, gene_info$gene_id)
} else {
  protein_seq <- NA
  cds_seq <- NA
}

if (!is.na(protein_seq)) {
  aa_seq <- unlist(strsplit(protein_seq, ""))
  aa_length <- length(aa_seq)
  aa_counts <- table(aa_seq)
  aa_composition <- paste(names(aa_counts), aa_counts, sep = ":", collapse = ", ")
  
  hydrophobic_aa <- c("A", "V", "I", "L", "M", "F", "W", "Y")
  charged_aa <- c("R", "H", "K", "D", "E")
  
  hydrophobic_count <- sum(aa_seq %in% hydrophobic_aa)
  charged_count <- sum(aa_seq %in% charged_aa)
  
  gene_info$protein_length <- aa_length
  gene_info$aa_composition <- aa_composition
  gene_info$hydrophobic_percent <- round(hydrophobic_count / aa_length * 100, 1)
  gene_info$charged_percent <- round(charged_count / aa_length * 100, 1)
}

# ============================
# 4. ESTATÍSTICAS E GRÁFICOS
# ============================
cat("INFO: Calculando estatísticas de expressão...\n")

# Função segura para estatísticas
safe_stat <- function(x, f, default = NA_real_) {
  if (all(is.na(x))) default else f(x, na.rm = TRUE)
}

stats_cat <- gene_data %>%
  group_by(stage_category) %>%
  summarise(
    mean_tpm = safe_stat(TPM, mean),
    sd_tpm = safe_stat(TPM, sd),
    median_tpm = safe_stat(TPM, median),
    min_tpm = safe_stat(TPM, min, Inf),
    max_tpm = safe_stat(TPM, max, -Inf),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    min_tpm = ifelse(is.infinite(min_tpm), NA, min_tpm),
    max_tpm = ifelse(is.infinite(max_tpm), NA, max_tpm)
  )

# CORREÇÃO AQUI: Usando which.max para evitar problemas com slice
if (nrow(stats_cat) > 0 && any(!is.na(stats_cat$mean_tpm))) {
  max_index <- which.max(stats_cat$mean_tpm)
  max_expression <- stats_cat$stage_category[max_index] %>% as.character()
} else {
  max_expression <- "N/A"
}

# Gráficos
p_box <- ggplot(gene_data, aes(x = stage_category, y = TPM, fill = stage_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = paste("Expressão de", gene_target),
       subtitle = gene_info$description,
       x = "Estágio de desenvolvimento",
       y = "TPM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_base, "boxplot.png"), p_box, width = 10, height = 6)

# ============================
# 5. RELATÓRIO COMPLETO
# ============================
cat("INFO: Gerando relatório completo...\n")

sink(paste0(output_base, "full_report.txt"))
cat("RELATÓRIO DE EXPRESSÃO E ANOTAÇÃO GÊNICA\n")
cat("=======================================\n\n")
cat("Gene analisado: ", gene_target, "\n")
cat("Descrição: ", gene_info$description, "\n\n")

if (!is.null(gene_info$location)) {
  cat("INFORMAÇÕES ESTRUTURAIS:\n")
  cat("------------------------\n")
  cat("Localização: ", gene_info$location, "\n")
  cat("Strand: ", gene_info$strand, "\n")
  cat("Comprimento do gene: ", gene_info$gene_length, " bp\n")
  cat("Comprimento do CDS: ", gene_info$cds_length, " bp\n")
  cat("Número de éxons: ", gene_info$exon_count, "\n")
  cat("Comprimentos dos éxons: ", gene_info$exon_lengths, "\n\n")
}

if (!is.null(gene_info$protein_length)) {
  cat("PROPRIEDADES DA PROTEÍNA:\n")
  cat("-------------------------\n")
  cat("Comprimento: ", gene_info$protein_length, " aminoácidos\n")
  cat("Composição de AA: ", gene_info$aa_composition, "\n")
  cat("% Aminoácidos hidrofóbicos: ", gene_info$hydrophobic_percent, "\n")
  cat("% Aminoácidos carregados: ", gene_info$charged_percent, "\n\n")
}

if (length(gene_info$go_terms) > 0) {
  cat("TERMOS GO (Gene Ontology):\n")
  cat("--------------------------\n")
  for (term in gene_info$go_terms) {
    cat("- ", term, "\n")
  }
  cat("\n")
}

cat("PERFIL DE EXPRESSÃO:\n")
cat("--------------------\n")
cat("Estágio com maior expressão: ", max_expression, "\n")
cat("Expressão média por estágio:\n")
print(knitr::kable(stats_cat, format = "simple"))
cat("\n")

cat("AMOSTRAS ANALISADAS:\n")
cat("--------------------\n")
print(knitr::kable(gene_data, format = "simple"))
sink()

# ============================
# 6. SALVAR SEQUÊNCIAS
# ============================
if (!is.na(protein_seq)) {
  write_lines(protein_seq, paste0(output_base, "protein.fasta"))
  cat("INFO: Sequência proteica salva\n")
}

if (!is.na(cds_seq)) {
  write_lines(cds_seq, paste0(output_base, "cds.fasta"))
  cat("INFO: Sequência CDS salva\n")
}

cat("\nANÁLISE CONCLUÍDA! Arquivos salvos com prefixo:", output_base, "\n")
