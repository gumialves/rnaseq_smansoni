#!/bin/bash
# run_integrated_pipeline.sh
#SBATCH --job-name=prem_integrated
#SBATCH --output=logs/prem_integrated_%A_%a.out
#SBATCH --error=logs/prem_integrated_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-3

set -euo pipefail

# =====================
# Configuração
# =====================
source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate r-analysis

PROJECT_DIR="/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"
SCRIPTS_DIR="$PROJECT_DIR/090-search-gene/scripts"
cd "$SCRIPTS_DIR"

# =====================
# Variáveis
# =====================
MAIN_PIPELINE="$SCRIPTS_DIR/analysis_pipeline.R"
TPM_FILE="$PROJECT_DIR/050-quantification/tpm_matrix.tsv"
GENES_FILE="genes.txt"

# Obter gene
GENE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENES_FILE" | tr -d '\r' | xargs)

echo "=== PIPELINE INTEGRADO COMPLETO ==="
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Gene: $GENE"
echo "Data: $(date)"
echo "Usando todos os módulos..."

# =====================
# Validar
# =====================
if [[ -z "$GENE" ]]; then
    echo "ERRO: Gene vazio"
    exit 1
fi

# Verificar se todos os módulos existem
MODULES=(
    "analysis_pipeline.R"
    "config_analysis.R" 
    "prem_utils.R"
    "qc_utils.R"
    "logging_utils.R"
    "reporting_utils.R"
    "monitoring_utils.R"
)

for module in "${MODULES[@]}"; do
    if [[ ! -f "$module" ]]; then
        echo "AVISO: Módulo $module não encontrado"
    else
        echo "✓ Módulo: $module"
    fi
done

# =====================
# Executar PIPELINE COMPLETO
# =====================
echo "Executando pipeline integrado..."
Rscript "$MAIN_PIPELINE" "$TPM_FILE" "$GENE" "$PROJECT_DIR"

echo "=== PIPELINE INTEGRADO CONCLUÍDO ==="
echo "Gene: $GENE"
echo "Status: FINALIZADO"
echo "Data: $(date)"
