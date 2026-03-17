#!/bin/bash
# run_anaysis.sh
#SBATCH --job-name=090_CompleteRun
#SBATCH --output=logs/090_complete_run_%A_%a.out
#SBATCH --error=logs/090_complete_run_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-27

set -euo pipefail

source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate r-analysis

PROJECT_DIR="/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"
SCRIPTS_DIR="$PROJECT_DIR/090-search-gene/scripts"
cd "$SCRIPTS_DIR"

PIPELINE_SCRIPT="$SCRIPTS_DIR/analysis_pipeline.R"
TPM_FILE="$PROJECT_DIR/050-quantification/tpm_matrix.tsv"
GENES_FILE="genes.txt"

GENE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENES_FILE" | tr -d '\r' | xargs)

echo "=== EXECUTANDO PIPELINE COMPLETO ==="
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Gene: $GENE"
echo "Data: $(date)"

# ==================================
if [[ -z "$GENE" ]]; then
    echo "ERRO: Gene vazio"
    exit 1
fi

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
    echo "ERRO: Pipeline script não encontrado: $PIPELINE_SCRIPT"
    exit 1
fi

# =================================
echo "Executando pipeline completo..."
Rscript "$PIPELINE_SCRIPT" "$TPM_FILE" "$GENE" "$PROJECT_DIR"

echo "=== PIPELINE CONCLUÍDO ==="
echo "Gene: $GENE"
echo "Status: FINALIZADO"
echo "Data: $(date)"
