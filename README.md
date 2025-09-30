# Pipeline de RNA-seq para Schistosoma mansoni

Este repositório contém um pipeline para análise de dados de RNA-seq do parasita *Schistosoma mansoni*, com foco na investigação da expressão gênica em diferentes estágios do ciclo de vida e entre sexos, bem como na atuação de reguladores epigenéticos.

## Objetivo do Projeto

O principal objetivo deste projeto é:
1. Investigar a **expressão gênica** (via mRNA) nos diferentes estágios do ciclo de vida e entre sexos do parasita *Schistosoma mansoni*.
2. Analisar a atuação de **reguladores epigenéticos** ao longo desses estágios e sexos, com base em seus perfis de expressão gênica.

## Estrutura do Pipeline

O pipeline é organizado em módulos sequenciais, cada um responsável por uma etapa específica da análise:

- `main.sh`: Script principal que orquestra a execução de todas as etapas do pipeline.
- `config/`: Diretório para arquivos de configuração e metadados do projeto.
- `envs/`: Contém definições de ambientes Conda para gerenciamento de dependências.
- `010-reference/`: Download e indexação do genoma de referência.
- `020-data-download/`: Download dos arquivos FASTQ do European Nucleotide Archive (ENA).
- `030-qc-fastq/`: Controle de qualidade (QC) e trimming dos dados FASTQ.
- `040-alignment/`: Etapa de alinhamento das leituras ao genoma de referência.
- `050-quantification/`: Quantificação da expressão gênica.
- `060-deg-analysis/`: Análise de Expressão Gênica Diferencial (DEG) utilizando DESeq2.
- `070-dtu-analysis/`: (Opcional) Análise de Uso Diferencial de Transcritos (DTU) e isoformas.
- `080-splicing/`: (Opcional) Análise de splicing alternativo.
- `090-search-gene/`: Scripts para busca e análise de genes específicos.

## Fluxo de Trabalho

1.  **Preparar referência**: Baixar FASTA e GTF do WormBase, criar índices para Salmon e STAR.
2.  **Baixar dados**: Download direto do ENA, verificar integridade via md5.
3.  **QC inicial**: FastQC individual + MultiQC global, trimming (fastp) se necessário.
4.  **Alinhamento e Quantificação**: Salmon (modo quasi-mapping) para contagens de transcritos, tximport para sumarizar ao nível de gene.
5.  **Análise DEG**: DESeq2 para comparação entre estágios/sexos, PCA, heatmaps, volcano plots.
6.  **DTU / isoformas (opcional)**: DRIMSeq, IsoformSwitchAnalyzeR para switches funcionais.
7.  **Splicing alternativo (opcional)**: rMATS, SUPPA2 ou MAJIQ para eventos canônicos/complexos.

## Referências

