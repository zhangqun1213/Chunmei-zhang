setwd("D:/bilibiliR/11_ssGSEA_Immune_Infiltration")  # Set working directory
rm(list = ls(all = TRUE))  # Clear all variables

# Install and load required R packages ---------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")

library(GSVA)  # Algorithm for ssGSEA
library(tidyverse)  # Data manipulation
library(ggpubr)  # Plotting
library(ggplot2)  # Plotting
library(pheatmap)  # For heatmaps

# ssGSEA is an algorithm that performs analysis based on reference gene sets.

# Prepare data: Transcriptome + Immune Cell Types + Group Data ---------------
DEG_expr <- read.csv("DEG_expr.csv", row.names = 1)
group <- read.csv("group.csv", row.names = 1)
markergenes <- read.csv("markergenes.csv")
table(markergenes$Cell.type)  # View immune cell types in the marker genes

# Data inspection
boxplot(DEG_expr, outline = FALSE, notch = FALSE, las = 2)
qx <- as.numeric(quantile(DEG_expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || 
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
  DEG_expr[which(DEG_expr <= 0)] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}  

# Data processing
# 'lapply' applies a function to each element of a list/vector.
geneset <- split(markergenes, markergenes$Cell.type)

im_geneset <- lapply(geneset, function(x) {
  gene <- x$Metagene
  unique(gene)
})
lapply(im_geneset[1:3], head)

save(im_geneset, file = "im_geneset.Rdata")

DEG_expr <- as.matrix(DEG_expr) 

# Perform ssGSEA analysis ----------------------------------------------------
result <- gsva(DEG_expr, im_geneset, method = "ssgsea")
result1 <- as.data.frame(t(result))

write.csv(result1, "ssGSEA_result.csv")

# Visualize results ----------------------------------------------------------
# Heatmap --------------------------------------------------------------------
pdf(file = "heatmap.pdf", width = 10, height = 8)
pheatmap(result1,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border = "black", main = "Heatmap",
         show_rownames = TRUE, show_colnames = TRUE,
         cexCol = 1, scale = 'row', cluster_col = TRUE, cluster_row = FALSE,
         angle_col = "45", legend = FALSE, legend_breaks = c(-3, 0, 3),
         fontsize_row = 10, fontsize_col = 10)
dev.off()  

# Grouped box plot -----------------------------------------------------------
# Data processing
data <- cbind(result1, group)
data <- data[, c(28, 29, 1:27)]
data <- pivot_longer(data, cols = 3:29, names_to = "celltype", values_to = "proportion")

# Plotting
pdf(file = "Grouped_Boxplot_with_Significance.pdf", width = 10, height = 8)
ggboxplot(data = data, 
          x = "celltype", y = "proportion", 
          combine = TRUE, merge = FALSE, color = "black", 
          fill = "group", palette = c("#1C3EDF", "#DF1C26"),
          xlab = "ssGSEA", ylab = "Expression", 
          bxp.errorbar = FALSE, bxp.errorbar.width = 0.2, 
          linetype = "solid", width = 0.8, notch = FALSE, 
          outlier.shape = 20, repel = TRUE, label.rectangle = TRUE, 
          ggtheme = theme_pubr()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  stat_compare_means(label = "p.signif", method = "t.test", 
                     ref.group = ".all.", hide.ns = FALSE, 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")))
dev.off()
