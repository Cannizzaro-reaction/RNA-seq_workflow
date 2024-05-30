# Install and load necessary R packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pheatmap", force = TRUE)
# set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com"))
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)

# Input of raw data
counts <- read.table("result/expression/raw_counts.txt", header = TRUE, row.names = 1)
colnames(counts) <- c("untreated_1", "untreated_2", "untreated_3", "gentamicin_1", "gentamicin_2", "gentamicin_3")

# Set up data frame for sample information
sample_info <- data.frame(
    row.names = colnames(counts),
    condition = c("untreated", "untreated", "untreated", "gentamicin", "gentamicin", "gentamicin")
)

# Set up DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)

# Data preprocessing (low quality data trimming & normalization)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

# Extract the results and arrange the results in order
res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Set up p value and log2 fold change value
pCutoff <- 10e-6
log2FCCutoff <- 2.0

########## Volcano Map ##########
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'Fold Change'),
                ylab = bquote(~-Log[10]~ 'P-value'),
                pCutoff = pCutoff,
                FCcutoff = log2FCCutoff,
                pointSize = 3.0,
                labSize = 3.0,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 0.5,
                legendLabels = c('NS', 'Log2 FC', 'Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)

# Screen out the genes with significant expression difference for heat map
sigGenes <- resOrdered[which(resOrdered$padj < pCutoff & abs(resOrdered$log2FoldChange) > log2FCCutoff), ]
sigGeneNames <- rownames(sigGenes)
sigCounts <- counts[sigGeneNames, ]
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
normalized_counts <- assay(vsd)
sigNormalizedCounts <- normalized_counts[sigGeneNames, ]

# set proper width and height
num_genes <- nrow(sigNormalizedCounts)
num_samples <- ncol(sigNormalizedCounts)
max_genes_per_page <- 50
max_samples_per_page <- 10
cellwidth <- 35
cellheight <- 300 / min(max_genes_per_page, num_genes) / 10

########## Heat Map ##########
pheatmap(sigNormalizedCounts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = num_genes <= max_genes_per_page, 
         show_colnames = num_samples <= max_samples_per_page, 
         annotation_col = sample_info,
         cellwidth = cellwidth,
         cellheight = cellheight)

########## MA Plot ##########
plotMA(res, ylim = c(-5, 5), main = "MA Plot", colSig = "darkorange", colNonSig = "dodgerblue", 
       xlab = "Mean of Normalized Counts", ylab = "Fold Change")

# Output
legend("topright", legend = c("Significant", "Non-Significant"), col = c("darkorange", "dodgerblue"), pch = 16)