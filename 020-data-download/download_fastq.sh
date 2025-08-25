#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=logs/fastq_download_%j.log
#SBATCH --error=logs/fastq_download_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=6:00:00

# Diretório de destino
OUTDIR="/scratch/Schisto-epigenetics/gustavo/fastq_ftp"
LINKS_FILE="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni/020-data-download/ftp_links.txt"

# Criar diretório se não existir
mkdir -p $OUTDIR

# Ir para o diretório
cd $OUTDIR

echo "Iniciando download dos FASTQ para $OUTDIR"


cat $LINKS_FILE | xargs -n 1 -P 8 wget -c

echo "Download finalizado!"
