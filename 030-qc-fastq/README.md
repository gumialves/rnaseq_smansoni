

---

Este documento descreve como é feito o **controle de qualidade (FastQC)** e o **trimming (Trim Galore)** nos dados de RNA-seq utilizados ao longo do pipleine, além de instruções para análise dos relatórios com **MultiQC**.

---

## Estrutura dos arquivos de entrada

* Os arquivos FASTQ brutos seguem o padrão:

  ```
  ERR11178266_1.fastq.gz
  ERR11178266_2.fastq.gz
  ```

  Onde:

  * `ERR11178266` = accession de run/lane no ENA.
  * `_1` e `_2` = leituras paired-end.

* Cada **amostra biológica** (ex.: `X_Eggs_R1`) pode ter **uma ou mais lanes (runs)**.
  Essas relações estão descritas no metadata:

  ```
  RNAseq_metadata.tsv
  ```

  informações importantes:

  * **Sample Name** (ex.: `X_Eggs_R1`, `F_Cercariae_R3`)
  * **Run accessions** (ex.: `ERR11178266, ERR11178356`) → runs que pertencem à mesma amostra.

---

## Fluxo de análise

### 1. **FastQC pré-trimming**

Script: `fastqc_pre.sh`

* Roda FastQC diretamente nos arquivos brutos (`ERRxxx_1/2.fastq.gz`).
* Cada run/lane é avaliado individualmente.
* Resultados salvos em:

---

### 2. **Trimming**

Script: `trimming.sh`

#### Etapas

1. **Concatenar lanes por amostra**

   * Todos os `_1.fastq.gz` da mesma amostra são unidos em um único arquivo.
   * Todos os `_2.fastq.gz` idem.
   * Exemplo:

     ```bash
     zcat ERR11178266_1.fastq.gz ERR11178356_1.fastq.gz | gzip > X_Eggs_R1_R1.fastq.gz
     zcat ERR11178266_2.fastq.gz ERR11178356_2.fastq.gz | gzip > X_Eggs_R1_R2.fastq.gz
     ```

2. **Rodar Trim Galore**

   ```bash
   trim_galore --paired \
       --quality 20 \
       --length 20 \
       --cores 8 \
       --output_dir trimmed/ \
       X_Eggs_R1_R1.fastq.gz X_Eggs_R1_R2.fastq.gz
   ```

3. **Renomear saída**

   * Arquivos resultantes:

     ```
     X_Eggs_R1_R1_trimmed.fastq.gz
     X_Eggs_R1_R2_trimmed.fastq.gz
     X_Eggs_R2_R1_trimmed.fastq.gz
     X_Eggs_R2_R2_trimmed.fastq.gz
     ```
* → ()_R1 - Representa a replicata 
* → ()_R1_R2 - Representa a run. ex: 'primera replica segunda lane'

#### Parâmetros utilizados

* `--paired` → indica que os dados são paired-end.
* `--quality 20` → corta bases de qualidade < Q20 no final da leitura.
* `--length 20` → descarta leituras menores que 20 nt após o trimming.
* `--cores 8` → usa 8 threads para paralelização.
* `--output_dir` → onde salvar os arquivos trimmados.

---

### 3. **FastQC pós-trimming**

Script: `fastqc_post.sh`

* Roda FastQC nos arquivos **já trimmados e renomeados**.
* Exemplo de entrada:

  ```
  X_Eggs_R1_R1_trimmed.fastq.gz
  X_Eggs_R1_R2_trimmed.fastq.gz
  ```
* Resultados salvos em:

  ```
  /scratch/Schisto-epigenetics/gustavo/fastqc_post/
  ```

---

### 4. **MultiQC**

Após rodar os passos acima, use o MultiQC para consolidar todos os relatórios:

```bash
# Pré-trimming
multiqc /scratch/Schisto-epigenetics/gustavo/fastqc_pre -o multiqc_pre

# Pós-trimming
multiqc /scratch/Schisto-epigenetics/gustavo/fastqc_post -o multiqc_post
```

---

## Resultados do MultiQC

### Pré-trimming (multiqc\_pre/)

* **Per base quality**: verifica se os arquivos têm queda de qualidade no final.
* **Adapter content**: confirma se havia contaminação por adaptadores.
* **Sequence duplication**: avalia duplicatas (pode estar alto em RNA-seq devido a genes muito expressos).
* **Overrepresented sequences**: identifica possíveis rRNA/artefatos.

### Pós-trimming (multiqc\_post/)

* **Comparar com pré-trimming**:

  * A qualidade por base deve melhorar (menos queda no final).
  * Adaptadores devem desaparecer ou reduzir fortemente.
  * Comprimento médio das leituras deve ser menor (cortes aplicados).
* **Duplicação**: não deve mudar muito, pois é mais característica da biblioteca biológica do que do trimming.

---
