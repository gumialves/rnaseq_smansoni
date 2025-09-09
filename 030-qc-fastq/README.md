# Etapa 030 - Qualidade das Leituras (QC) e Trimming

Esta etapa avalia a qualidade das leituras **antes e depois do trimming**, e gera um relatório consolidado usando **MultiQC**.  

## Estrutura da Pasta
030-qc-fastq/
│── run_fastqc_pre.slurm # QC antes do trimming
│── run_trimming.slurm # Trimming com Trim Galore
│── run_fastqc_post.slurm # QC depois do trimming
│── run_multiqc.slurm # Consolidação com MultiQC
│── logs/ # Logs das submissões SLURM
│── fastqc_pre/ # Resultados FastQC (pré-trimming)
│── fastqc_post/ # Resultados FastQC (pós-trimming)
│── multiqc_report/ # Relatório consolidado MultiQC

markdown


## Entrada Necessária
- Arquivos **FASTQ** 
- Arquivo `RNAseq_metadata.tsv` com informações das amostras:
Internal Sample ID Library Name Sample Name ENA Sample Accession i7 Tag Sequence i5 Tag Sequence Lane 1 filename Lane 1 ENA Run Accession ...


## Passos do Pipeline
### 1. Rodar QC inicial (antes do trimming)

```bash
sbatch 030-qc-fastq/run_fastqc_pre.sh
Saída: 030-qc-fastq/fastqc_pre
```
###2. Rodar trimming com Trim Galore
```bash

sbatch 030-qc-fastq/run_trimming.slurm
Remove adaptadores

Realiza QC pós-trimming automático

Saída: /scratch/Schisto-epigenetics/gustavo/trimmed
```
###3. Rodar QC pós-trimming explicitamente
```bash

sbatch 030-qc-fastq/run_fastqc_post.slurm
Saída: 030-qc-fastq/fastqc_post
```
4. Consolidar com MultiQC
```bash
sbatch 030-qc-fastq/run_multiqc.slurm
Saída: 030-qc-fastq/multiqc_report/multiqc_report.html
```
# Como interpretar os resultados
## MultiQC compara pré vs pós trimming:

Verifique gráficos de per-base quality: deve melhorar após trimming
Adapter Content deve reduzir após trimming
GC Content deve permanecer estável
Use este relatório para decidir se trimming foi efetivo e necessário.
