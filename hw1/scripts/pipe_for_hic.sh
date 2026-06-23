#!/bin/bash
set -euo pipefail

THREADS=4
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

declare -A SAMPLES=(
  [MoPh7]=85
  [MoPh11]=86
  [MoPh14]=87
  [MoPh15]=88
)

mkdir -p data/raw
mkdir -p data/trimmed
mkdir -p data/juicer
mkdir -p results/fastqc
mkdir -p results/cutadapt
mkdir -p results/hic

echo "=== Hi-C pipeline start ==="

for sample in "${!SAMPLES[@]}"
do
    S=${SAMPLES[$sample]}

    echo ">>> Processing ${sample} (S${S})"

    # Если финальный hic уже существует — пропускаем образец
    if [ -f "results/hic/${sample}.inter_30.hic" ]; then
        echo ">>> SKIP ${sample} (already completed)"
        continue
    fi


    # 1. DOWNLOAD

    if [ ! -f "data/raw/${sample}_R1.fastq.gz" ]; then
        wget --no-check-certificate \
          -O "data/raw/${sample}_R1.fastq.gz" \
          "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20${sample}_S${S}_L001_R1_001.fastq.gz"
    fi

    if [ ! -f "data/raw/${sample}_R2.fastq.gz" ]; then
        wget --no-check-certificate \
          -O "data/raw/${sample}_R2.fastq.gz" \
          "https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics/Copy%20of%20${sample}_S${S}_L001_R2_001.fastq.gz"
    fi

    # 2. FASTQC

    fastqc \
      "data/raw/${sample}_R1.fastq.gz" \
      "data/raw/${sample}_R2.fastq.gz" \
      -o results/fastqc

    # 3. CUTADAPT

    if [ ! -f "data/trimmed/${sample}_R1.trimmed.fastq.gz" ]; then

        cutadapt \
          -q 20 \
          -m 70 \
          -a "${ADAPTER}" \
          -A "${ADAPTER}" \
          -o "data/trimmed/${sample}_R1.trimmed.fastq.gz" \
          -p "data/trimmed/${sample}_R2.trimmed.fastq.gz" \
          "data/raw/${sample}_R1.fastq.gz" \
          "data/raw/${sample}_R2.fastq.gz" \
          > "results/cutadapt/${sample}.cutadapt.log" 2>&1
    fi

    # 4. JUICER INPUT STRUCTURE

    mkdir -p "data/juicer/${sample}/fastq"

    ln -sf \
      "$(pwd)/data/trimmed/${sample}_R1.trimmed.fastq.gz" \
      "data/juicer/${sample}/fastq/${sample}_R1.fastq.gz"

    ln -sf \
      "$(pwd)/data/trimmed/${sample}_R2.trimmed.fastq.gz" \
      "data/juicer/${sample}/fastq/${sample}_R2.fastq.gz"

    # 5. JUICER RUN

    bash tools/juicer/scripts/juicer.sh \
      -D "$(pwd)/tools/juicer" \
      -d "$(pwd)/data/juicer/${sample}" \
      -g T2T_human \
      -z "$(pwd)/data/reference/T2T_human.fna" \
      -p "$(pwd)/data/reference/chrom.sizes" \
      -y "$(pwd)/data/reference/restriction_sites_DpnII.txt" \
      -s DpnII \
      -t "${THREADS}"

    # 6. FINAL HIC FILE

    if [ -f "data/juicer/${sample}/aligned/inter_30.hic" ]; then

        cp \
          "data/juicer/${sample}/aligned/inter_30.hic" \
          "results/hic/${sample}.inter_30.hic"

        echo ">>> DONE ${sample}"

    else
        echo ">>> ERROR: inter_30.hic not found for ${sample}"
        exit 1
    fi

done

echo "=== ALL SAMPLES FINISHED ==="
