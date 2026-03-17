#!/bin/bash
# generate_final_dashboard.sh
#SBATCH --job-name=prem_dashboard
#SBATCH --output=logs/prem_dashboard_%j.out
#SBATCH --error=logs/prem_dashboard_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

# =====================
# Configuração
# =====================
source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate r-analysis

PROJECT_DIR="/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"
SCRIPTS_DIR="$PROJECT_DIR/090-search-gene/scripts"
RESULTS_DIR="$PROJECT_DIR/090-search-gene/results"
cd "$SCRIPTS_DIR"

# =====================
# Função Principal
# =====================
main() {
    echo "=== GERANDO DASHBOARD FINAL ==="
    echo "Data: $(date)"
    echo "Diretório de resultados: $RESULTS_DIR"
    
    # Verificar se existem resultados
    if [[ ! -d "$RESULTS_DIR" ]]; then
        echo "ERRO: Diretório de resultados não encontrado: $RESULTS_DIR"
        exit 1
    fi
    
    # Contar genes processados
    total_genes=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name "Smp_*" | wc -l)
    echo "Total de genes processados: $total_genes"
    
    if [[ $total_genes -eq 0 ]]; then
        echo "AVISO: Nenhum gene processado encontrado"
        exit 0
    fi
    
    # =====================
    # Executar Análise Consolidada em R
    # =====================
    echo "Executando análise consolidada..."
    
    Rscript --vanilla - "
    # Carregar bibliotecas
    library(dplyr)
    library(ggplot2)
    library(jsonlite)
    library(tidyr)
    
    cat('=== ANÁLISE CONSOLIDADA DE TODOS OS GENES ===\\n')
    
    # Configurar caminhos
    results_dir <- '$RESULTS_DIR'
    scripts_dir <- '$SCRIPTS_DIR'
    
    # Função para ler summary de um gene
    read_gene_summary <- function(gene_dir) {
        summary_file <- file.path(gene_dir, 'analysis_summary.json')
        if (file.exists(summary_file)) {
            tryCatch({
                summary_data <- fromJSON(summary_file)
                return(data.frame(
                    gene = summary_data\$gene,
                    status = summary_data\$status,
                    analysis_date = as.POSIXct(summary_data\$analysis_date),
                    n_files = length(summary_data\$output_files),
                    mean_expression = if (!is.null(summary_data\$summary_stats\$mean_expression)) 
                                       summary_data\$summary_stats\$mean_expression else NA,
                    n_samples = if (!is.null(summary_data\$summary_stats\$n_samples)) 
                                 summary_data\$summary_stats\$n_samples else NA,
                    stringsAsFactors = FALSE
                ))
            }, error = function(e) {
                cat('Erro ao ler:', summary_file, '-', e\$message, '\\n')
                return(NULL)
            })
        }
        return(NULL)
    }
    
    # Coletar dados de todos os genes
    cat('Coletando dados de', $total_genes, 'genes...\\n')
    gene_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
    gene_dirs <- gene_dirs[grepl('Smp_', basename(gene_dirs))]
    
    all_data <- list()
    for (gene_dir in gene_dirs) {
        gene_data <- read_gene_summary(gene_dir)
        if (!is.null(gene_data)) {
            all_data[[length(all_data) + 1]] <- gene_data
        }
    }
    
    # Combinar todos os dados
    if (length(all_data) > 0) {
        consolidated <- bind_rows(all_data)
        
        cat('\\n=== ESTATÍSTICAS GERAIS ===\\n')
        cat('Total de genes analisados:', nrow(consolidated), '\\n')
        cat('Data da primeira análise:', as.character(min(consolidated\$analysis_date)), '\\n')
        cat('Data da última análise:', as.character(max(consolidated\$analysis_date)), '\\n')
        cat('Média de arquivos por gene:', round(mean(consolidated\$n_files), 1), '\\n')
        cat('Média de expressão (TPM):', round(mean(consolidated\$mean_expression, na.rm = TRUE), 2), '\\n')
        
        # Salvar dados consolidados
        write.csv(consolidated, file.path(results_dir, 'consolidated_analysis.csv'), row.names = FALSE)
        cat('✓ Dados consolidados salvos: consolidated_analysis.csv\\n')
        
        # =====================
        # Gerar Gráficos do Dashboard
        # =====================
        cat('Gerando gráficos do dashboard...\\n')
        
        # 1. Distribuição de Expressão
        p1 <- ggplot(consolidated, aes(x = mean_expression)) +
            geom_histogram(bins = 20, fill = 'steelblue', alpha = 0.7) +
            labs(title = 'Distribuição da Expressão Média dos Genes',
                 x = 'TPM Médio', y = 'Número de Genes') +
            theme_bw()
        
        # 2. Número de Arquivos por Gene
        p2 <- ggplot(consolidated, aes(x = n_files)) +
            geom_histogram(bins = 15, fill = 'darkorange', alpha = 0.7) +
            labs(title = 'Número de Arquivos Gerados por Gene',
                 x = 'Número de Arquivos', y = 'Número de Genes') +
            theme_bw()
        
        # 3. Timeline das Análises
        p3 <- ggplot(consolidated, aes(x = analysis_date)) +
            geom_histogram(bins = 20, fill = 'darkgreen', alpha = 0.7) +
            labs(title = 'Distribuição Temporal das Análises',
                 x = 'Data da Análise', y = 'Número de Genes') +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        # Salvar gráficos
        ggsave(file.path(results_dir, 'dashboard_expression.png'), p1, width = 8, height = 6, dpi = 300)
        ggsave(file.path(results_dir, 'dashboard_files.png'), p2, width = 8, height = 6, dpi = 300)
        ggsave(file.path(results_dir, 'dashboard_timeline.png'), p3, width = 8, height = 6, dpi = 300)
        
        cat('✓ Gráficos do dashboard gerados\\n')
        
        # =====================
        # Gerar Relatório HTML Simples
        # =====================
        cat('Gerando relatório HTML...\\n')
        
        html_content <- paste('
        <!DOCTYPE html>
        <html>
        <head>
            <title>Dashboard PREM - Análise Consolidada</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                .header { background: #f4f4f4; padding: 20px; border-radius: 5px; }
                .stats { display: flex; justify-content: space-around; margin: 20px 0; }
                .stat-box { background: white; padding: 15px; border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); text-align: center; }
                .images { display: flex; flex-direction: column; gap: 20px; }
                img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; }
                table { width: 100%; border-collapse: collapse; margin: 20px 0; }
                th, td { padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }
                th { background-color: #f2f2f2; }
            </style>
        </head>
        <body>
            <div class=\"header\">
                <h1>Dashboard PREM - Análise Consolidada</h1>
                <p>Data de geração: ', Sys.Date(), '</p>
                <p>Total de genes analisados: ', nrow(consolidated), '</p>
            </div>
            
            <div class=\"stats\">
                <div class=\"stat-box\">
                    <h3>', nrow(consolidated), '</h3>
                    <p>Genes Analisados</p>
                </div>
                <div class=\"stat-box\">
                    <h3>', round(mean(consolidated\$mean_expression, na.rm = TRUE), 2), '</h3>
                    <p>TPM Médio</p>
                </div>
                <div class=\"stat-box\">
                    <h3>', round(mean(consolidated\$n_files), 1), '</h3>
                    <p>Arquivos/Gene</p>
                </div>
            </div>
            
            <div class=\"images\">
                <h2>Visualizações</h2>
                <div>
                    <h3>Distribuição da Expressão</h3>
                    <img src=\"dashboard_expression.png\" alt=\"Distribuição da Expressão\">
                </div>
                <div>
                    <h3>Arquivos por Gene</h3>
                    <img src=\"dashboard_files.png\" alt=\"Arquivos por Gene\">
                </div>
                <div>
                    <h3>Timeline das Análises</h3>
                    <img src=\"dashboard_timeline.png\" alt=\"Timeline das Análises\">
                </div>
            </div>
            
            <div>
                <h2>Top 10 Genes com Maior Expressão</h2>
                <table>
                    <tr><th>Gene</th><th>TPM Médio</th><th>Arquivos</th><th>Data</th></tr>'
        ,
        paste(sapply(head(consolidated[order(-consolidated\$mean_expression), ], 10), function(row) {
            paste('<tr><td>', row\$gene, '</td><td>', round(row\$mean_expression, 2), 
                  '</td><td>', row\$n_files, '</td><td>', as.character(row\$analysis_date), '</td></tr>', sep = '')
        }), collapse = ''),
        '
                </table>
            </div>
            
            <div>
                <h2>Lista Completa de Genes</h2>
                <table>
                    <tr><th>Gene</th><th>Status</th><th>TPM Médio</th><th>Arquivos</th></tr>'
        ,
        paste(sapply(consolidated[order(consolidated\$gene), ], function(row) {
            paste('<tr><td>', row\$gene, '</td><td>', row\$status, 
                  '</td><td>', round(row\$mean_expression, 2), '</td><td>', row\$n_files, '</td></tr>', sep = '')
        }), collapse = ''),
        '
                </table>
            </div>
        </body>
        </html>'
        , sep = '')
        
        writeLines(html_content, file.path(results_dir, 'dashboard.html'))
        cat('✓ Dashboard HTML gerado: dashboard.html\\n')
        
    } else {
        cat('AVISO: Nenhum dado válido encontrado para consolidar\\n')
    }
    
    cat('\\\\n=== DASHBOARD FINALIZADO ===\\\\n')
    "

    echo "=== DASHBOARD CONCLUÍDO ==="
    echo "Arquivos gerados em: $RESULTS_DIR/"
    echo "- consolidated_analysis.csv"
    echo "- dashboard.html" 
    echo "- dashboard_expression.png"
    echo "- dashboard_files.png"
    echo "- dashboard_timeline.png"
}

# Executar
main "$@"
