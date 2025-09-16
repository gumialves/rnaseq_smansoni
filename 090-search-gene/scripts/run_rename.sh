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

./rename.py
