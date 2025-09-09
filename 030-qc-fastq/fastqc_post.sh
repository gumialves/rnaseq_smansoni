#!/bin/bash
#SBATCH --job-name=fastqc_post
#SBATCH --output=logs/qc_post/fastqc_post_%A_%a.out
#SBATCH --error=logs/qc_post/fastqc_post_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-180%5

set -euo pipefail

source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate
conda activate rna-tools

PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
SCRATCH_DIR="/scratch/Schisto-epigenetics/gustavo"
TRIM_DIR="${SCRATCH_DIR}/trimmed"
QC_POST_DIR="${SCRATCH_DIR}/fastqc_post"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$QC_POST_DIR"

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

line_num=$((SLURM_ARRAY_TASK_ID + 1))


sample_info=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3 "\t" $8 "\t" $13}' "$META")

if [ -z "$sample_info" ]; then
    echo "[INFO] Nenhuma amostra encontrada para SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# Divide a informação extraída
sample_name=$(echo "$sample_info" | cut -f1)
r1_acc=$(echo "$sample_info" | cut -f2)
r2_acc=$(echo "$sample_info" | cut -f3)

# Função para encontrar arquivos trimmed com padrões flexíveis
find_trimmed() {
    local base_acc="$1"
    local read_end="$2"
    local patterns=()

    # Padrões baseados nos arquivos reais encontrados
    patterns+=("${TRIM_DIR}/${base_acc}_${read_end}_val_${read_end}.fq.gz")
    patterns+=("${TRIM_DIR}/${base_acc}_${read_end}_trimmed.fq.gz")
    
    # Padrões alternativos
    patterns+=("${TRIM_DIR}/${base_acc}_val_${read_end}.fq.gz")
    patterns+=("${TRIM_DIR}/${base_acc}_R${read_end}_val_${read_end}.fq.gz")
    patterns+=("${TRIM_DIR}/${base_acc}_R${read_end}_trimmed.fq.gz")
    
    # Padrões com substituição de caracteres especiais
    local clean_acc=$(echo "$base_acc" | sed 's/#/_/g')
    patterns+=("${TRIM_DIR}/${clean_acc}_${read_end}_val_${read_end}.fq.gz")
    patterns+=("${TRIM_DIR}/${clean_acc}_${read_end}_trimmed.fq.gz")
    patterns+=("${TRIM_DIR}/${clean_acc}_val_${read_end}.fq.gz")
    patterns+=("${TRIM_DIR}/${clean_acc}_R${read_end}_val_${read_end}.fq.gz")
    patterns+=("${TRIM_DIR}/${clean_acc}_R${read_end}_trimmed.fq.gz")

    # Busca por qualquer arquivo que contenha o accession
    patterns+=("${TRIM_DIR}/*${base_acc}*val_${read_end}*.fq.gz")
    patterns+=("${TRIM_DIR}/*${base_acc}*trimmed*.fq.gz")
    patterns+=("${TRIM_DIR}/*${clean_acc}*val_${read_end}*.fq.gz")
    patterns+=("${TRIM_DIR}/*${clean_acc}*trimmed*.fq.gz")

    # Tenta encontrar o arquivo
    for pattern in "${patterns[@]}"; do
        # Usa nullglob para evitar expansão de padrões não encontrados
        shopt -s nullglob
        local files=($pattern)
        if [ ${#files[@]} -gt 0 ]; then
            echo "${files[0]}"
            return 0
        fi
    done

    echo ""
    return 1
}

# Encontra os arquivos trimmed R1 e R2
echo "[INFO] Procurando arquivos trimmed para ${sample_name} (R1: ${r1_acc}, R2: ${r2_acc})"
R1=$(find_trimmed "$r1_acc" "1")
R2=$(find_trimmed "$r2_acc" "2")

if [[ -n "$R1" && -n "$R2" && -f "$R1" && -f "$R2" ]]; then
    echo "[INFO] Arquivos encontrados:"
    echo "[INFO] R1: ${R1}"
    echo "[INFO] R2: ${R2}"
    echo "[INFO] Processando amostra ${SLURM_ARRAY_TASK_ID}: ${sample_name}"

    fastqc "$R1" "$R2" --outdir "$QC_POST_DIR" --threads "${SLURM_CPUS_PER_TASK}"

    echo "[OK] FastQC concluído para ${sample_name}."
else
    echo "[ERRO] Arquivos trimmed não encontrados para ${sample_name}"
    echo "[DEBUG] Tentando listar arquivos no diretório TRIM_DIR..."
    ls -la "${TRIM_DIR}" | grep -E "(${r1_acc}|${r2_acc})" | head -20 || true
    exit 1
fi
