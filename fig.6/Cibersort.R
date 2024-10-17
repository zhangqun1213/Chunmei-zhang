setwd("D:/R/Immune_Infiltration")  # Set working directory

rm(list = ls(all = TRUE))  # Clear all variables

# Install and load the CIBERSORT package -------------------------------------
# install.packages('devtools')
# library(devtools)
# if(!require(CIBERSORT)) devtools::install_github("Moonerss/CIBERSORT")
# devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
library(ggplot2)  # For plotting
library(pheatmap)  # For heatmaps
library(ggpubr)  # For stacked bar plots
library(reshape2)  # For data transformation
library(tidyverse)  # For data manipulation

# Prepare data: Transcriptome + Immune Cell Types + Group Data ----------------
DEG_expr <- read.csv("DEG_expr.csv", row.names = 1)
group <- read.csv("group.csv", row.names = 1)

# Data inspection
boxplot(DEG_expr, outline = FALSE, notch = FALSE, las = 2)
qx <- as.numeric(quantile(DEG_expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || 
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  DEG_expr[DEG_expr <= 0] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}   

# Normalize data if means are inconsistent -----------------------------------
library(limma)
DEG_expr <- normalizeBetweenArrays(DEG_expr)
boxplot(DEG_expr, outline = FALSE, notch = FALSE, las = 2)

# Immune cell type reference data --------------------------------------------
LM22_local <- read.table("LM22.txt", header = TRUE, row.names = 1, sep = "\t")
data(LM22)  # Load LM22 data from CIBERSORT
all(LM22 == LM22_local)  # Check if the downloaded file matches the built-in LM22 data

# Perform CIBERSORT analysis -------------------------------------------------
result <- cibersort(sig_matrix = LM22, mixture_file = DEG_expr, perm = 1, QN = TRUE)
result <- as.data.frame(result)
write.csv(result, "cibersort_result.csv")

# Visualize results ----------------------------------------------------------
result1 <- result[, 1:ncol(LM22)]
result1 <- result1[, apply(result1, 2, function(x) { sum(x) > 0 })]

# Heatmap
pdf(file = "Heatmap.pdf", width = 10, height = 8)
pheatmap(result1,
         color = colorRampPalette(c(rep("skyblue", 3.5), "#FEFCFB", rep("#ED5467", 3.5)))(50),
         border = "skyblue", main = "Heatmap",
         show_rownames = TRUE, show_colnames = TRUE,
         cexCol = 1, scale = 'row', cluster_col = TRUE, cluster_row = FALSE,
         angle_col = "45", legend = FALSE, legend_breaks = c(-3, 0, 3),
         fontsize_row = 10, fontsize_col = 10)
dev.off()

# Stacked bar plot -----------------------------------------------------------
identical(rownames(result1), group$Samples)
data <- cbind(rownames(result1), result1)
colnames(data)[1] <- "Samples"
data <- melt(data, id.vars = "Samples")
colnames(data) <- c('Samples', 'celltype', 'proportion')

mycolors <- c('#D4E2A7', '#88D7A4', '#A136A1', '#BAE8BC', '#C757AF', 
              '#DF9FCE', '#D5E1F1', '#305691', '#B6C2E7', '#E8EFF7',
              '#9FDFDF', '#EEE0F5', '#267336', '#98CEDD', '#CDE2EE',
              '#DAD490', '#372E8A', '#4C862D', '#81D5B0', '#BAE8C9',
              '#A7DCE2', '#AFDE9C')

pdf(file = "stacked_bar_chart.pdf", width = 10, height = 8)
ggplot(data, aes(Samples, proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = mycolors) +
  ggtitle("Proportion of Immune Cells") + theme_gray() +
  theme(axis.ticks.length = unit(3, 'mm'), axis.title.x = element_text(size = 11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  guides(fill = guide_legend(title = "Types of Immune Cells"))
dev.off()

# Box plot grouped by samples ------------------------------------------------
data1 <- cbind(result1, group)
data1 <- pivot_longer(data1, cols = 3:24, names_to = "celltype", values_to = "proportion")

pdf(file = "Grouped_Boxplot.pdf", width = 10, height = 8)
ggboxplot(data = data1, x = "celltype", y = "proportion", fill = "group",
          palette = NULL, title = "TME Cell Composition", ylab = "Cell Composition",
          ggtheme = theme_pubr()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

# Box plot for a single group ------------------------------------------------
data2 <- data1[data1$group == 'Tumor', ]

pdf(file = "Single_Group_Boxplot.pdf", width = 10, height = 8)
ggboxplot(data = data2, x = "celltype", y = "proportion", fill = "celltype",
          title = "TME Cell Composition", ggtheme = theme_pubr()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

# Enhanced grouped box plot with significance --------------------------------
data3 <- pivot_longer(data1, cols = 3:24, names_to = "celltype", values_to = "proportion")

pdf(file = "Enhanced_Boxplot_with_Significance.pdf", width = 10, height = 8)
ggboxplot(data = data3, x = "celltype", y = "proportion", fill = "group", 
          palette = c("#81D5B0", "#ED5462"), title = "TME Cell Composition",
          ylab = "Cell Composition", ggtheme = theme_pubr()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  stat_compare_means(label = "p.signif", method = "t.test", aes(group = group), hide.ns = TRUE)
dev.off()
