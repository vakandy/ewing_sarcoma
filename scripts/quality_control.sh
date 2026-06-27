#!/bin/bash
# Quality control with FastQC + MultiQC

set -e

cd "$(dirname "$0")/.."

SAMPLES=("SRR36910393" "SRR36910395" "SRR36910397" "SRR36910429" "SRR36910431" "SRR36910433")

mkdir -p qc/fastqc_reports
mkdir -p qc/multiqc_report

for SAMPLE in "${SAMPLES[@]}"; do
    echo ">>> Running FastQC on $SAMPLE..."
    fastqc "data/raw_fastq/$SAMPLE/${SAMPLE}_1.fastq" \
           "data/raw_fastq/$SAMPLE/${SAMPLE}_2.fastq" \
           --outdir qc/fastqc_reports/
done

echo ">>> Running MultiQC to summarize all reports..."
multiqc qc/fastqc_reports/ --outdir qc/multiqc_report/

echo "Quality control complete."
