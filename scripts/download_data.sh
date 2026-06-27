#!/bin/bash
# Download RNA-seq FASTQ files from SRA
# Dataset: GSE317136 (Ewing sarcoma cell lines vs hMSC)

set -e  # exit on error

cd "$(dirname "$0")/.."  # move to project root (ewing_sarcoma/)

# Ewing sarcoma samples (A-673, control/no doxycycline)
EWING_SAMPLES=("SRR36910393" "SRR36910395" "SRR36910397")

# Normal hMSC samples
HMSC_SAMPLES=("SRR36910429" "SRR36910431" "SRR36910433")

ALL_SAMPLES=("${EWING_SAMPLES[@]}" "${HMSC_SAMPLES[@]}")

mkdir -p data/raw_fastq

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    echo ">>> Downloading $SAMPLE..."
    prefetch "$SAMPLE" -O data/raw_fastq/

    echo ">>> Converting $SAMPLE to FASTQ..."
    fasterq-dump "data/raw_fastq/$SAMPLE/$SAMPLE.sra" \
        --outdir "data/raw_fastq/$SAMPLE/"
done

echo "Download complete for all samples."
