setwd("D:/R/DESEq2_Differential_Expression")  # Set working directory

# Clear variables ------------------------------------------------------------
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)
library(DESeq2)

# Load the expression matrix -------------------------------------------------
load("count_and_TPM.Rdata")

expr <- mRNA_count
expr <- expr[rowMeans(expr) > 1, ]

# Create grouping information ------------------------------------------------
table(substr(colnames(expr), 14, 16))  # Check the number of tumor and normal samples
Tumor <- grep('01A', colnames(expr))  # Find the positions of tumor samples
Tumor
group <- factor(c(rep("Normal", times = 11), rep("Tumor", times = 152)))  # Create a factor variable for groups

Data <- data.frame(row.names = colnames(expr),  # Create a grouping data frame
                   group = group)
View(Data)
# The gene expression matrix (expr) and sample grouping data frame (Data) are ready

# Start differential expression analysis -------------------------------------
# Step 1: Build the DESeq2 object (dds)
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = Data,
                              design = ~ group)
# Step 2: Perform differential expression analysis
dds2 <- DESeq(dds)
res <- results(dds2, contrast = c("group", "Tumor", "Normal"))  # Tumor first, control second
## Or: res = results(dds)
res <- res[order(res$pvalue), ]
summary(res)
my_result <- as.data.frame(res)  # Convert to a data frame for easier viewing
my_result <- na.omit(my_result)  # Remove rows with missing values

# Step 3: Save the differential expression results
library(dplyr)
my_result$Gene_symbol <- rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol', colnames(my_result)[1:(dim(my_result)[2] - 1)], everything())
rownames(my_result) <- NULL
write.csv(my_result, file = "my_result_DESeq2.csv")

# DEG filtering --------------------------------------------------------------
my_result$regulate <- ifelse(my_result$padj > 0.05, "unchanged",
                       ifelse(my_result$log2FoldChange > 2, "up-regulated",
                              ifelse(my_result$log2FoldChange < -2, "down-regulated", "unchanged")))
table(my_result$regulate)

# Extract upregulated and downregulated genes
DEG_deseq2 <- subset(my_result, padj < 0.05 & abs(log2FoldChange) > 2) 
upgene <- DEG_deseq2[DEG_deseq2$regulate == 'up-regulated', ]
downgene <- DEG_deseq2[DEG_deseq2$regulate == 'down-regulated', ]
write.csv(DEG_deseq2, file = "DEG_deseq2.csv")

# Visualize the results ------------------------------------------------------
plot(my_result$log2FoldChange, -log2(my_result$padj))  # Basic volcano plot
library(ggplot2)
library(ggrepel)

p <- ggplot(data = my_result, aes(x = log2FoldChange, y = -log10(padj), color = regulate)) + 
  geom_point(shape = 16, size = 2) + 
  theme_set(theme_set(theme_bw(base_size = 20))) + 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(plot.title = element_text(size = 15, hjust = 2.5)) + 
  theme_classic() +
  scale_colour_manual(values = c('#86D47D', '#DFE0DF', '#C34B99')) +
  geom_vline(xintercept = c(-2, 2), lty = 4, col = "gray", lwd = 0.8) +  # logFC threshold
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "gray", lwd = 0.6) +  # adj.p.val threshold
  labs(title = 'TCGA-ESCA') +  # Title at the top left corner
  annotate("text", x = upgene$log2FoldChange[1:3], y = -log10(upgene$padj[1:3]), 
           label = upgene$Gene_symbol[1:3], size = 5.0) +
  annotate("text", x = downgene$log2FoldChange[1:3], y = -log10(downgene$padj[1:3]), 
           label = downgene$Gene_symbol[1:3], size = 5.0)
plot(p)

# Alternative volcano plot ---------------------------------------------------
library(ggpubr)
my_result$v <- -log10(my_result$padj)  # Create a new column 'v'

# Select genes to display
my_select <- subset(my_result, padj < 1e-20 & abs(log2FoldChange) > 4)
my_select <- my_select$Gene_symbol

# Plotting
ggscatter(my_result,
          x = "log2FoldChange",
          y = "v",
          ylab = "-log10(adjusted p-value)",
          size = 2,
          color = "regulate",
          palette = c('#86D47D', '#DFE0DF', '#C34B99'),
          label = "Gene_symbol",
          label.select = c("GPR155"))
