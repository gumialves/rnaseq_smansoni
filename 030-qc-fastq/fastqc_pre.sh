#!/bin/bash
#SBATCH --job-name=fastqc_pre
#SBATCH --output=logs/fastqc_pre_%A_%a.out
#SBATCH --error=logs/fastqc_pre_%A_%a.err
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
RAW_DIR="${SCRATCH_DIR}/fastq_ftp"
QC_PRE_DIR="${SCRATCH_DIR}/fastqc_pre"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$QC_PRE_DIR"

# Verifica se é um job array válido
if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

# Pula o cabeçalho (linha 1) e pega a linha correspondente ao array task ID
line_num=$((SLURM_ARRAY_TASK_ID + 1))

# Extrai informações da amostra usando awk
sample_info=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3 "\t" $8 "\t" $13}' "$META")

if [ -z "$sample_info" ]; then
    echo "[INFO] Nenhuma amostra encontrada para SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# Divide a informação extraída
sample_name=$(echo "$sample_info" | cut -f1)
r1_acc=$(echo "$sample_info" | cut -f2)
r2_acc=$(echo "$sample_info" | cut -f3)

# Função para encontrar arquivos com padrões flexíveis
find_fastq() {
    local base_acc="$1"
    local read_end="$2"
    local patterns=()
    
    # Gera vários padrões possíveis
    patterns+=("${RAW_DIR}/${base_acc}_${read_end}.fastq.gz")
    patterns+=("${RAW_DIR}/${base_acc}_R${read_end}_001.fastq.gz")
    patterns+=("${RAW_DIR}/${base_acc}_${read_end}.fq.gz")
    patterns+=("${RAW_DIR}/${base_acc}_R${read_end}.fastq.gz")
    patterns+=("${RAW_DIR}/${base_acc}_R${read_end}.fq.gz")
    
    # Padrões com substituição de caracteres especiais
    local clean_acc=$(echo "$base_acc" | sed 's/#/_/g')
    patterns+=("${RAW_DIR}/${clean_acc}_${read_end}.fastq.gz")
    patterns+=("${RAW_DIR}/${clean_acc}_R${read_end}_001.fastq.gz")
    patterns+=("${RAW_DIR}/${clean_acc}_${read_end}.fq.gz")
    patterns+=("${RAW_DIR}/${clean_acc}_R${read_end}.fastq.gz")
    patterns+=("${RAW_DIR}/${clean_acc}_R${read_end}.fq.gz")
    
    # Busca por qualquer arquivo que contenha o accession
    patterns+=("${RAW_DIR}/*${base_acc}*_${read_end}*.fastq.gz")
    patterns+=("${RAW_DIR}/*${base_acc}*_R${read_end}_*.fastq.gz")
    patterns+=("${RAW_DIR}/*${clean_acc}*_${read_end}*.fastq.gz")
    patterns+=("${RAW_DIR}/*${clean_acc}*_R${read_end}_*.fastq.gz")
    
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

# Encontra os arquivos R1 e R2
echo "[INFO] Procurando arquivos para ${sample_name} (R1: ${r1_acc}, R2: ${r2_acc})"
R1=$(find_fastq "$r1_acc" "1")
R2=$(find_fastq "$r2_acc" "2")

if [[ -n "$R1" && -n "$R2" && -f "$R1" && -f "$R2" ]]; then
    echo "[INFO] Arquivos encontrados:"
    echo "[INFO] R1: ${R1}"
    echo "[INFO] R2: ${R2}"
    echo "[INFO] Processando amostra ${SLURM_ARRAY_TASK_ID}: ${sample_name}"
    
    fastqc "$R1" "$R2" --outdir "$QC_PRE_DIR" --threads "${SLURM_CPUS_PER_TASK}"
    
    echo "[OK] FastQC concluído para ${sample_name}."
else
    echo "[ERRO] Arquivos FASTQ não encontrados para ${sample_name}"
    echo "[DEBUG] Tentando listar arquivos no diretório RAW_DIR..."
    ls -la "${RAW_DIR}" | grep -E "(${r1_acc}|${r2_acc})" | head -20 || true
    exit 1
fi
