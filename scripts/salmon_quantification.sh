#!/bin/bash
# Quantify transcript expression with Salmon

set -e

cd "$(dirname "$0")/.."

SAMPLES=("SRR36910393" "SRR36910395" "SRR36910397" "SRR36910429" "SRR36910431" "SRR36910433")

mkdir -p quantification

for SAMPLE in "${SAMPLES[@]}"; do
    echo ">>> Quantifying $SAMPLE with Salmon..."
    salmon quant \
        -i data/reference/salmon_index \
        -l A \
        -1 "data/raw_fastq/$SAMPLE/${SAMPLE}_1.fastq" \
        -2 "data/raw_fastq/$SAMPLE/${SAMPLE}_2.fastq" \
        --validateMappings \
        -o "quantification/$SAMPLE"
done

echo "Quantification complete for all samples."
