#!/bin/bash
# Build reference transcriptome + Salmon index

set -e

cd "$(dirname "$0")/.."

mkdir -p data/reference
cd data/reference

echo ">>> Downloading reference genome (GRCh38)..."
wget -nc https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo ">>> Downloading gene annotation (GTF)..."
wget -nc https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

echo ">>> Decompressing reference files..."
gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -k Homo_sapiens.GRCh38.109.gtf.gz

echo ">>> Extracting transcriptome with gffread..."
gffread Homo_sapiens.GRCh38.109.gtf \
    -g Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -w transcriptome.fa

echo ">>> Building Salmon index..."
salmon index \
    -t transcriptome.fa \
    -i salmon_index \
    -p 4

echo "Reference preparation complete."
