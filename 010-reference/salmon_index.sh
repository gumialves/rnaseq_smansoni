#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=../010-reference/logs/salmon_index_%j.out
#SBATCH --error=../010-reference/logs/salmon_index_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail


source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate rna-tools


REF_DIR="../010-reference"
DATA_DIR="${REF_DIR}/data"
INDEX_DIR="${REF_DIR}/salmon_index"
LOG_DIR="${REF_DIR}/logs"

mkdir -p "$DATA_DIR" "$INDEX_DIR" "$LOG_DIR"


TRANSCRIPTS_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa.gz"
TRANSCRIPTS_FA="${DATA_DIR}/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa.gz"


if [ ! -f "$TRANSCRIPTS_FA" ]; then
    echo "[INFO] Baixando FASTA de transcritos..."
    wget -O "$TRANSCRIPTS_FA" "$TRANSCRIPTS_URL"
else
    echo "[INFO] FASTA de transcritos já existe em $TRANSCRIPTS_FA"
fi

echo "[INFO] Criando índice do Salmon..."
salmon index \
    -t "$TRANSCRIPTS_FA" \
    -i "$INDEX_DIR" \
    -p $SLURM_CPUS_PER_TASK \
    -k 31 #Rever o k-mer-size aplicado

echo "[OK] Índice do Salmon criado em $INDEX_DIR"
