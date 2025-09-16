#!/bin/bash
#SBATCH --job-name=prem_analysis
#SBATCH --output=logs/prem/prem_%A.out
#SBATCH --error=logs/prem/prem_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00

set -euo pipefail

# =====================
# Ativar conda
# =====================
source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate r-analysis

# =====================
# Caminhos
# =====================
PREM_GRAPHS="prem_graphs.R"
PREM_ANNOT="prem_annot.R"
PREM_STATS="prem_stats.R"
TPM_FILE="../../050-quantification/tpm_matrix.tsv"
GENES_FILE="genes.txt"

mkdir -p logs/prem ../results

# =====================
# Loop sobre genes
# =====================
while read -r gene; do
    echo "=== Processando gene ${gene} ==="
    mkdir -p "../results/${gene}"

    Rscript "$PREM_GRAPHS" "$TPM_FILE" "$gene"
    Rscript "$PREM_ANNOT" "$TPM_FILE" "$gene"
    Rscript "$PREM_STATS" "../results/${gene}/TPM_data_long.csv" "../results/${gene}"

    echo "=== Concluído: ${gene} ==="
done < "$GENES_FILE"

echo "[OK] Todas as análises concluídas. Resultados em ../results/"
