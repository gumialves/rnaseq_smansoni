#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --output=logs/deseq2_%j.out
#SBATCH --error=logs/deseq2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate r-analysis

Rscript deseq2_analysis.R
