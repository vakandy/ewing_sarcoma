# Load Libraries
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)


# Set up file paths
base_dir <- "/home/vivek/Documents/ewing_sarcoma"

samples <- c("SRR36910393", "SRR36910395", "SRR36910397",  # Ewing sarcoma
             "SRR36910429", "SRR36910431", "SRR36910433")  # hMSC normal

files <- file.path(base_dir, "quantification", samples, "quant.sf")
names(files) <- samples
print(file.exists(files))


# tx2gene mapping from GTF
txdb <- makeTxDbFromGFF(
  file.path(base_dir, "data/reference/Homo_sapiens.GRCh38.109.gtf"),
  format = "gtf"
)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
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


# Filter significant DEGs
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1.5)
cat("Number of significant DEGs:", nrow(sig), "\n")

dir.create(file.path(base_dir, "results"), showWarnings = FALSE)

write.csv(as.data.frame(res), 
          file.path(base_dir, "results/DESeq2_all_results.csv"))
write.csv(as.data.frame(sig), 
          file.path(base_dir, "results/DESeq2_significant_DEGs.csv"))


# Heatmap
top_n <- 50
top_genes <- head(rownames(sig), top_n)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = top_genes,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first") # Convert Ensembl IDs to gene symbols

gene_labels <- ifelse(is.na(gene_symbols), top_genes, gene_symbols) # Fall back to Ensembl ID if no symbol is found

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat))) # Build expression matrix

rownames(mat_scaled) <- gene_labels # Replace row names with gene symbols

sample_labels <- c("SRR36910393" = "SRR36910393 (Ewing_1)",
                   "SRR36910395" = "SRR36910395 (Ewing_2)",
                   "SRR36910397" = "SRR36910397 (Ewing_3)",
                   "SRR36910429" = "SRR36910429 (hMSC_1)",
                   "SRR36910431" = "SRR36910431 (hMSC_2)",
                   "SRR36910433" = "SRR36910433 (hMSC_3)") # Sample labels include SRR number + sample label for columns

colnames(mat_scaled) <- sample_labels[colnames(mat_scaled)]

annotation_col <- data.frame(
  Condition = coldata$condition,
  row.names = sample_labels[rownames(coldata)]
)

ann_colors <- list(
  Condition = c(Ewing = "#D62728", hMSC = "#1F77B4")
)

heatmap_colors <- colorRampPalette(c("#053061", "#2166AC", "#92C5DE", 
                                     "#FFFFFF", 
                                     "#F4A582", "#B2182B", "#67001F"))(100) # Color palette

pheatmap(mat_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         color = heatmap_colors,
         breaks = seq(-2, 2, length.out = 101),
         main = paste("Top", top_n, "DEGs: Ewing Sarcoma vs hMSC"),
         filename = file.path(base_dir, paste0("results/heatmap_top", top_n, "_DEGs_final.pdf")),
         width = 10,
         height = 10
) # Heatmap save as PDF

pheatmap(mat_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         color = heatmap_colors,
         breaks = seq(-2, 2, length.out = 101),
         main = paste("Top", top_n, "DEGs: Ewing Sarcoma vs hMSC"),
         filename = file.path(base_dir, paste0("results/heatmap_top", top_n, "_DEGs_final.png")),
         width = 10,
         height = 10,
         res = 300
) # Heatmap saved as PNG


