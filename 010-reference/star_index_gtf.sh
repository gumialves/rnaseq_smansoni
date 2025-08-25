#!/bin/bash
#SBATCH --job-name=star_index_gtf
#SBATCH --output=logs/star_index_gtf_%j.out
#SBATCH --error=logs/star_index_gtf_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=08:00:00

set -euo pipefail

PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
REF_DIR="${PROJECT_DIR}/010-reference"
REF_DATA="${REF_DIR}/data"
STAR_INDEX_DIR="${REF_DIR}/star_index_gtf"

FA="${REF_DATA}/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"
GFF3="${REF_DATA}/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
GTF="${REF_DATA}/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gtf"


source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate

conda activate rna-tools


mkdir -p "${REF_DIR}/logs" "${STAR_INDEX_DIR}"

if [[ -f "$GFF3" && ! -f "$GTF" ]]; then
  echo "[INFO] Convertendo GFF3 para GTF..."
  gffread "$GFF3" -T -o "$GTF"
fi

echo "[INFO] Criando índice STAR com GTF..."
STAR --runMode genomeGenerate \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${FA}" \
  --sjdbGTFfile "${GTF}" \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM 170000000000

echo "[OK] Índice STAR criado em ${STAR_INDEX_DIR}"
