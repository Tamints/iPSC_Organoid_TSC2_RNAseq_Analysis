library(stringr)

sample_names <- colnames(expr)

meta <- data.frame(
  sample_id = sample_names,
  condition = ifelse(str_detect(sample_names, "CTRL"), "CTRL", "TSC"),
  timepoint = str_extract(sample_names, "Day[0-9]+"),
  row.names = sample_names,
  stringsAsFactors = FALSE
)

# Convert condition to factor and set CTRL as reference
meta$condition <- factor(meta$condition)
meta$condition <- relevel(meta$condition, ref = "CTRL")

# Sanity check
head(meta)

## PCA

# Load required libraries
library(DESeq2)
library(ggplot2)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = expr,
  colData = meta,
  design = ~ condition
)

# Pre-filter very low-expression genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Variance stabilizing transformation (for PCA)
vsd <- vst(dds, blind = TRUE)

# Generate PCA data
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA Plot
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of GSE247367 Bulk RNA-seq")
