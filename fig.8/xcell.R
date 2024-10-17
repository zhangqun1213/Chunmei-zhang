# Set working directory ------------------------------------------------------
setwd("D:/bilibiliR/16_xCell_Immune_Infiltration")

# Load required R packages ---------------------------------------------------
library(devtools)
# install_github('dviraran/xCell')
library(xCell)
library(tidyverse)  # Data processing
library(readxl)  # Read Excel files
library(ggstatsplot)  # Plotting
library(ggpubr)  # Plotting
library(pheatmap)  # Heatmap plotting
library(survival)  # Survival analysis
library(survminer)  # Survival curves

# Load expression matrix -----------------------------------------------------
load("TPM.Rdata")
# Remove duplicates based on sample names
expr_use <- mRNA_TPM_new
expr_use = expr_use[, sort(colnames(expr_use))]
k = !duplicated(str_sub(colnames(expr_use), 1, 12))
table(k)
expr_use = expr_use[, k]
table(substr(colnames(expr_use), 14, 16))  # View tumor and normal sample counts
colnames(expr_use) <- substr(colnames(expr_use), 1, 12)
expr_data <- as.matrix(expr_use)  # Convert to matrix

# Retrieve xCell gene list (10,800 genes)
genes <- xCell.data$genes

# Ensure the dataset contains at least 5,000 of the xCell genes for accuracy
xCell.data$signatures  # View gene sets
geneset <- read_excel("489_Cell_Type_Gene_Features.xlsx")

# Begin analysis -------------------------------------------------------------
scores <- rawEnrichmentAnalysis(expr_data, 
                                signatures = xCell.data$signatures, 
                                genes = xCell.data$genes, 
                                parallel.sz = 4, 
                                parallel.type = "SOCK")
tscores <- transformScores(scores, 
                           fit.vals = xCell.data$spill$fv, 
                           scale = TRUE)
result <- spillOver(tscores, 
                    K = xCell.data$spill$K, 
                    alpha = 0.5)

# One-step analysis
result <- xCellAnalysis(expr_data)

# Process the data -----------------------------------------------------------
result1 <- as.data.frame(t(result))
result1 <- result1[, c(65:67, 1:64)]
summary(result1[, 1:3])
result1 <- result1[, 1:3]

# Data exploration -----------------------------------------------------------
View(expr_use)  # Expression matrix
View(result1)  # Immune, stroma, and microenvironment scores
surv_data <- read.csv("survival_data.csv", row.names = 1)  # Survival data

# Group data by immune score and survival status -----------------------------
expr_use <- t(expr_use)
combined_data1 <- cbind(result1, expr_use)
rownames(surv_data) <- surv_data$submitter_id
surv_data <- surv_data[, -1]
combined_data2 <- merge(surv_data, combined_data1, by = 0, all = FALSE)

data_use <- combined_data2[, 1:36]
rownames(data_use) <- data_use$Row.names
data_use <- data_use[, -1]
table(data_use$status)

# Plot boxplots --------------------------------------------------------------
ggplot(data_use, aes(x = status, y = ImmuneScore, color = status)) + 
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.1)) +
  theme(text = element_text(size = 12), axis.title = element_text(size = 18)) + 
  xlab("Status") + 
  ylab("Immune Score") +
  scale_x_discrete(labels = c("Alive", "Dead")) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("#60B2F2", "#E12939")) + 
  stat_compare_means(method = "t.test")

ggboxplot(data_use, "status", "StromaScore", 
          color = "status", palette = c("#E12939", "#2F74E5")) + 
  stat_compare_means(method = "t.test")

ggbetweenstats(data = data_use, x = status, y = MicroenvironmentScore, plot.type = "box")

# Plot heatmap ---------------------------------------------------------------
Immune_median_score <- median(data_use$ImmuneScore)
data_use$Immune_group <- cut(data_use$ImmuneScore, 
                             breaks = c(-Inf, Immune_median_score, Inf), 
                             labels = c("low", "high"))
data_use_1 <- data_use[, c(36, 6:35)]
row_index <- order(data_use_1$Immune_group)
data_use_1 <- data_use_1[row_index, ]
table(data_use_1$Immune_group)
heat_data <- t(data_use_1[, 2:31])
scaled_heat_data <- scale(heat_data)
group <- data_use_1[, 1, drop = FALSE]

pheatmap(scaled_heat_data,
         color = colorRampPalette(c("#0084E4", "white", "#E12939"))(50),
         annotation_col = group,
         cluster_row = TRUE,
         cluster_col = FALSE,
         border = "black",
         main = "Heatmap",
         show_rownames = TRUE,
         show_colnames = FALSE,
         cexCol = 1,
         scale = 'row',
         angle_col = 90,
         legend = TRUE,
         legend_breaks = c(-5, 0, 5),
         fontsize_row = 10,
         fontsize_col = 10)

# Plot survival curves -------------------------------------------------------
data_use_surv <- data_use[, c(36, 1, 2)]
time_year <- data_use_surv$time / 365
data_use_surv$time_year <- time_year
data_use_surv$status <- ifelse(data_use_surv$status == "Dead", 1, 0)
fit <- survfit(Surv(time_year, status) ~ Immune_group, data = data_use_surv)

ggsurvplot(fit, 
           data = data_use_surv,
           surv.median.line = "hv",
           size = 1,
           cex.lab = 2,
           break.time.by = 1,
           xlim = c(0, 6),
           axis.title.x = element_text(size = 5), 
           axis.title.y = element_text(size = 5),
           palette = c("#F04972", "#7BB0E0"),
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE,
           xlab = "Follow-up years",  
           ylab = "Survival probability",
           risk.table.col = "strata",
           legend.labs = c("Immune score low", "Immune score high"),
           risk.table.height = 0.3,
           ggtheme = theme_bw())
