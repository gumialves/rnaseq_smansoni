#!/bin/bash

# pipeline_config.sh - Configurações centrais para o pipeline de RNA-seq de S. mansoni

# Diretório base do projeto
export PROJECT_DIR="$(dirname "$(dirname "$(readlink -f "$0")")")"

# Diretórios principais
export REF_DIR="${PROJECT_DIR}/010-reference"
export DATA_DIR="${REF_DIR}/data"
export LOG_DIR="${REF_DIR}/logs"
export ENVS_DIR="${PROJECT_DIR}/envs"

# URLs para o genoma e transcritos de referência
export GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa.gz"
export TRANSCRIPTS_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa.gz"
export GFF3_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gz"

# Nomes de arquivos de referência
export GENOME_FA_GZ="$(basename "$GENOME_URL")"
export GENOME_FA="${GENOME_FA_GZ%.gz}"
export TRANSCRIPTS_FA_GZ="$(basename "$TRANSCRIPTS_URL")"
export TRANSCRIPTS_FA="${TRANSCRIPTS_FA_GZ%.gz}"
export GFF3_GZ="$(basename "$GFF3_URL")"
export GFF3="${GFF3_GZ%.gz}"
export GTF="${GFF3%.gff3}.gtf"

# Caminhos completos para os arquivos de referência
export REF_GENOME_FA="${DATA_DIR}/${GENOME_FA}"
export REF_TRANSCRIPTS_FA="${DATA_DIR}/${TRANSCRIPTS_FA}"
export REF_GFF3="${DATA_DIR}/${GFF3}"
export REF_GTF="${DATA_DIR}/${GTF}"

# Diretórios de índice
export SALMON_INDEX_DIR="${REF_DIR}/salmon_index"
export STAR_INDEX_DIR="${REF_DIR}/star_index"
export STAR_INDEX_GTF_DIR="${REF_DIR}/star_index_gtf"

# Configurações de recursos (pode ser sobrescrito por SLURM)
export THREADS=${SLURM_CPUS_PER_TASK:-8}
export MEMORY=${SLURM_MEM:-32G}

# Ambiente Conda
export CONDA_BASE="$(conda info --base)"
export RNA_TOOLS_ENV="rna-tools"
export R_ANALYSIS_ENV="r-analysis"

# Função para ativar ambiente Conda
activate_conda_env() {
    local env_name="$1"
    if [ -z "$CONDA_BASE" ]; then
        echo "[ERRO] Conda não encontrado. Certifique-se de que o Conda está instalado e no PATH."
        exit 1
    fi
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "$env_name"
    if [ $? -ne 0 ]; then
        echo "[ERRO] Falha ao ativar o ambiente Conda: $env_name"
        exit 1
    fi
}

# Função para verificar a existência de um comando
check_command() {
    local cmd="$1"
    if ! command -v "$cmd" &> /dev/null; then
        echo "[ERRO] Comando '$cmd' não encontrado. Certifique-se de que o ambiente Conda correto está ativado e a ferramenta instalada."
        exit 1
    fi
}

