#!/bin/bash

#SBATCH --job-name=100_integrated_analysis
#SBATCH --output=logs/100_integrated_analysis_%A_%a.out
#SBATCH --error=logs/100_integrated_analysis_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail

# =====================
# Configuração - Isso é só um template, favor lembrar de atualizar isso pro main, inclusive os jobs do sbatch
# =====================
source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate r-analysis

echo "Running 01 - Prepare Data"
Rscript 01-prepare-data/prepare.R

echo ""
echo ""
echo "Running 02 - Life Cycle"
Rscript 02-life-cycle/lifecycle_LTR.R

echo ""
echo ""
echo "Running 03 - Sex Dimorphism"
Rscript 03-sex-dimorphism/dimorphism_stage_interact.R

echo ""
echo ""
echo "Running 04 - Visualization"
Rscript 04-funtional-groups/vizualizegroups.R

echo ""
echo ""
echo "Running 05"
Rscript 05-candidate-scoring/integrative.R

echo "Pipeline finished."
