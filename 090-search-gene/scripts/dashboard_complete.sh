#!/bin/bash
# run_complete_analysis.sh
# Script mestre que executa o array job E depois o dashboard

set -euo pipefail

PROJECT_DIR="/home/ra236875@bio.ib.unicamp.br/rnaseq_smansoni"
SCRIPTS_DIR="$PROJECT_DIR/090-search-gene/scripts"
cd "$SCRIPTS_DIR"

echo "=== EXECUÇÃO COMPLETA DA ANÁLISE PREM ==="
echo "Data de início: $(date)"

# 1. Executar o array job
echo "1. Submetendo array job..."
ARRAY_JOB_ID=$(sbatch --parsable run_integrated_pipeline.sh)

echo "   Array job submetido com ID: $ARRAY_JOB_ID"

# 2. Executar o dashboard DEPOIS do array job terminar
echo "2. Agendando dashboard para execução após array job..."
DASHBOARD_JOB_ID=$(sbatch --parsable --dependency=afterany:$ARRAY_JOB_ID generate_final_dashboard.sh)

echo "   Dashboard agendado com ID: $DASHBOARD_JOB_ID"

# 3. Mostrar informações
echo ""
echo "=== RESUMO DA SUBMISSÃO ==="
echo "Array Job ID: $ARRAY_JOB_ID"
echo "Dashboard Job ID: $DASHBOARD_JOB_ID"
echo ""
echo "Para monitorar:"
echo "  squeue -u $USER"
echo ""
echo "Para ver os resultados do dashboard:"
echo "  ls -la $PROJECT_DIR/090-search-gene/results/dashboard.*"
echo ""
echo "Data: $(date)"
