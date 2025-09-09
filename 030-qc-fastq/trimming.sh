#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=logs/trim/trimming_%A_%a.out
#SBATCH --error=logs/trim/trimming_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=18:00:00
#SBATCH --array=1-180%5

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
RAW_DIR="/scratch/Schisto-epigenetics/gustavo/fastq_ftp"
TRIM_DIR="/scratch/Schisto-epigenetics/gustavo/trimmed"
FASTQC_POST_DIR="/scratch/Schisto-epigenetics/gustavo/fastqc_post"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$TRIM_DIR" "$FASTQC_POST_DIR" "logs"

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

# Pula o cabeçalho (linha 1) e pega a linha correspondente ao array task ID
line_num=$((SLURM_ARRAY_TASK_ID + 1))

# Extrai informações da amostra - usando colunas 8 e 13 para os accessions ENA
sample_info=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3 "\t" $8 "\t" $13}' "$META")

if [ -z "$sample_info" ]; then
    echo "[INFO] Nenhuma amostra encontrada para SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# Divide a informação extraída
sample_name=$(echo "$sample_info" | cut -f1)
r1_acc=$(echo "$sample_info" | cut -f2)
r2_acc=$(echo "$sample_info" | cut -f3)

# Construir os caminhos dos arquivos
R1="${RAW_DIR}/${r1_acc}_1.fastq.gz"
R2="${RAW_DIR}/${r2_acc}_2.fastq.gz"

if [[ -f "$R1" && -f "$R2" ]]; then
    echo "[INFO] Processando amostra ${SLURM_ARRAY_TASK_ID}: ${sample_name}"
    echo "[INFO] Arquivo R1: ${R1}"
    echo "[INFO] Arquivo R2: ${R2}"
    
    # Executar Trim Galore
    trim_galore --paired \
        --quality 20 \
        --fastqc \
        --fastqc_args "--outdir ${FASTQC_POST_DIR} --threads ${SLURM_CPUS_PER_TASK}" \
        --cores ${SLURM_CPUS_PER_TASK} \
        --output_dir "$TRIM_DIR" \
        --length 20 \
        --stringency 1 \
        --clip_R1 15 \
        --clip_R2 15 \
        --three_prime_clip_R1 5 \
        --three_prime_clip_R2 5 \
        "$R1" "$R2"
    
    # Mover e renomear relatórios de trimming
    mv "${TRIM_DIR}/${r1_acc}_1_val_1.fq.gz" "${TRIM_DIR}/${sample_name}_R1_trimmed.fastq.gz"
    mv "${TRIM_DIR}/${r2_acc}_2_val_2.fq.gz" "${TRIM_DIR}/${sample_name}_R2_trimmed.fastq.gz"
    mv "${TRIM_DIR}/${r1_acc}_1_val_1.fastqc.html" "${TRIM_DIR}/${sample_name}_R1_trimming_report.html"
    mv "${TRIM_DIR}/${r2_acc}_2_val_2.fastqc.html" "${TRIM_DIR}/${sample_name}_R2_trimming_report.html"
    
    echo "[OK] Trimming concluído para ${sample_name}."
else
    echo "[ERRO] Arquivos FASTQ não encontrados para ${sample_name}"
    echo "[DEBUG] R1 esperado: ${R1}"
    echo "[DEBUG] R2 esperado: ${R2}"
    exit 1
fi
