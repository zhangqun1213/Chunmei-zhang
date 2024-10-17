setwd("D:/R/Differential_Expression_Analysis")  # Set working directory

# Load required R packages ---------------------------------------------------
library(limma)
library(ggplot2)  # For plotting volcano plots
library(pheatmap)  # For plotting heatmaps

# Input expression matrix and group file -------------------------------------
expr_data <- read.table("Matrix_of_expression.txt", header = TRUE,
                        row.names = 1, sep = "\t")
group <- read.csv("group.csv", header = TRUE, row.names = 1, sep = ",")

# Construct the design matrix ------------------------------------------------
design <- model.matrix(~0 + factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)

# Construct the contrast matrix ----------------------------------------------
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

# Build linear fit models ----------------------------------------------------
fit <- lmFit(expr_data, design)  # Non-linear least squares fitting
fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)  # Use empirical Bayes to adjust variance in t-tests
DEG <- topTable(fit2, coef = 1, n = Inf)
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                   ifelse(DEG$logFC > 1, "up-regulated",
                          ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
table(DEG$regulate)
write.table(data.frame(gene_symbol = rownames(DEG), DEG), 
            file = "DEG_result.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)

# Separate upregulated and downregulated genes -------------------------------
DE_1_0.05 <- DEG[DEG$P.Value < 0.05 & abs(DEG$logFC) > 1, ]
upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated", ]
downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated", ]
write.csv(upGene_1_0.05, "upGene_1_0.05.csv")
write.csv(downGene_1_0.05, "down-regulated.csv")

# Plot the volcano plot ------------------------------------------------------
pdf("volcano.pdf")
ggplot(DEG, aes(x = logFC, y = -log10(P.Value))) +  # x-axis: logFC, y-axis: -log10(P.Value)
  geom_point(alpha = 0.6, size = 3.5, aes(color = regulate)) +  # Set point transparency and size
  ylab("-log10(P.Value)") +  # y-axis label
  scale_color_manual(values = c("blue", "grey", "red")) +  # Define point colors
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +  # logFC threshold lines
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +  # p-value threshold line
  theme_bw()  # Apply black and white theme
dev.off()

# Plot another volcano plot using ggvolcano ----------------------------------
library(ggVolcano)
Genes <- rownames(DEG)
DEG$Genes <- Genes  # Add row names as a new column
ggvolcano(data = DEG, x = "logFC", y = "P.Value", label = "Genes",
          label_number = 10, output = FALSE,
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC", y_lab = "-Log10P.Value")

# Plot the heatmap -----------------------------------------------------------
DEG_genes <- DEG[DEG$P.Value < 0.05 & abs(DEG$logFC) > 1, ]
DEG_gene_expr <- expr_data[rownames(DEG_genes), ]
pdf("pheatmap.pdf")
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Define color palette
         scale = "row",  # Row normalization
         border_color = NA,  # No border lines
         fontsize = 10,  # Set font size
         show_rownames = FALSE) 
dev.off()
