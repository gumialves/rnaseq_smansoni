#!/bin/bash
#SBATCH --job-name=fastqc_post
#SBATCH --output=logs/qc_post/fastqc_post_%A_%a.out
#SBATCH --error=logs/qc_post/fastqc_post_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-75%5

set -euo pipefail

# =====================
# Ativar conda
# =====================
source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate
conda activate rna-tools

# =====================
# Caminhos
# =====================
PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
SCRATCH_DIR="/scratch/Schisto-epigenetics/gustavo"
TRIM_DIR="${SCRATCH_DIR}/trimmed"
QC_POST_DIR="${SCRATCH_DIR}/fastqc_post"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$QC_POST_DIR" logs/qc_post


if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

line_num=$((SLURM_ARRAY_TASK_ID + 1))

# Extrai sample_name (coluna 3)
sample_name=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3}' "$META")

if [ -z "$sample_name" ]; then
    echo "[INFO] Nenhuma amostra encontrada para linha $line_num"
    exit 0
fi


replicates=( $(ls "${TRIM_DIR}/${sample_name}_R[0-9]_R1_trimmed.fastq.gz" 2>/dev/null || true) )

if [ ${#replicates[@]} -eq 0 ]; then
    echo "[ERRO] Nenhuma replicata encontrada para ${sample_name}"
    exit 1
fi


for r1_file in "${replicates[@]}"; do
    rep=$(basename "$r1_file" | sed -E "s/${sample_name}_R([0-9]+)_R1_trimmed.fastq.gz/\1/")
    r2_file="${TRIM_DIR}/${sample_name}_R${rep}_R2_trimmed.fastq.gz"

    if [[ -f "$r1_file" && -f "$r2_file" ]]; then
        echo "[INFO] Processando ${sample_name} (rep${rep})"
        fastqc "$r1_file" "$r2_file" --outdir "$QC_POST_DIR" --threads "${SLURM_CPUS_PER_TASK}"
    else
        echo "[ERRO] Arquivos trimmados não encontrados para ${sample_name} (rep${rep})"
        exit 1
    fi
done

echo "[OK] FastQC concluído para ${sample_name}"
