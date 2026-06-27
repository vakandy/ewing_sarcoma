# Ewing Sarcoma GO Enrichment Analysis
#
# This script regenerates res and sig directly (same code as in
# ewing_sarcoma_deseq2_analysis.R) instead of reading saved CSV files.


# Load libraries
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tidyr)


# Set up file paths
base_directory <- "/home/vivek/Documents/ewing_sarcoma"
samples <- c("SRR36910393", "SRR36910395", "SRR36910397",  # Ewing sarcoma
             "SRR36910429", "SRR36910431", "SRR36910433")  # hMSC normal
files <- file.path(base_directory, "quantification", samples, "quant.sf")
names(files) <- samples
print(file.exists(files))


# tx2gene mapping 
txdb <- makeTxDbFromGFF(
  file.path(base_directory, "data/reference/Homo_sapiens.GRCh38.109.gtf"),
  format = "gtf"
)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
head(tx2gene)


# Import Salmon counts with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
cat("Dimensions of count matrix:", dim(txi$counts), "\n")


# Sample metadata
coldata <- data.frame(
  sample = samples,
  condition = factor(c("Ewing", "Ewing", "Ewing", "hMSC", "hMSC", "hMSC")),
  row.names = samples
)
print(coldata)


# DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Ewing", "hMSC"))
res <- res[order(res$padj), ]
summary(res)


# Filtering significant DEGs
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1.5)
cat("Number of significant DEGs:", nrow(sig), "\n")
dir.create(file.path(base_directory, "results"), showWarnings = FALSE)


# GO enrichment analysis
sig_genes <- rownames(sig) # Significant gene list (Ensembl IDs)
universe_genes <- rownames(res)[!is.na(res$padj)] # all genes tested (non-NA padj)
go_results <- enrichGO(
  gene          = sig_genes,
  universe      = universe_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",          # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE           # Convert IDs to gene symbols in output
)

go_df <- as.data.frame(go_results)
write.csv(go_df, file.path(base_directory, "results/GO_enrichment_full.csv"), row.names = FALSE) #Create CSV
cat("Number of significantly enriched GO terms:", nrow(go_df), "\n")


# GO enrichment summary table (top 20 terms)
go_df_top20 <- go_df %>% slice_head(n = 20)
go_table_top20 <- go_df_top20 %>%
  select(ID, Count, Description, GeneRatio, BgRatio, RichFactor,
         FoldEnrichment, zScore, pvalue, p.adjust, qvalue)
go_table_top20
