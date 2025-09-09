#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --error=logs/multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00


source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate
conda activate rna-tools


# =======================================
# 2. Caminhos
# =======================================
SCRATCH_DIR="/scratch/Schisto-epigenetics/gustavo"

PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
FILE_DIR="${PROJECT_DIR}/030-qc-fastq"
QC_DIR="${SCRATCH_DIR}/fastqc_post"
OUT_DIR="${FILE_DIR}/multiqc_report_post"

mkdir -p "$OUT_DIR"

# =======================================
# 3. Executar MultiQC
# =======================================
echo "[INFO] Iniciando MultiQC..."
multiqc "${QC_DIR}" -o "$OUT_DIR" --interactive

echo "[INFO] Relatório HTML gerado em: ${OUT_DIR}/multiqc_report.html"
