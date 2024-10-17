# Set the working directory --------------------------------------------------
setwd("D:/R/Heatmap_Plot")

# Install required packages --------------------------------------------------
# Check if the 'psych' package is installed
if (!require("psych", quietly = TRUE)) {
  # If not installed, use install.packages to install
  install.packages("psych")
}

# Check if the 'pheatmap' package is installed
if (!require("pheatmap", quietly = TRUE)) {
  # If not installed, use install.packages to install
  install.packages("pheatmap")
}

# Load required libraries ----------------------------------------------------
library(psych)
library(pheatmap)

# Load example data ----------------------------------------------------------
load("D:/bilibiliR/41_Significant_Heatmap/TCGA-ESCA_mrna_expr_tpm.rdata")
cli <- read.csv("clinical_data.csv", header = TRUE, row.names = 1, sep = ",")

# Extract part of the data for visualization ---------------------------------
mydata <- mrna_expr_tpm[1:30,]  # Select genes of interest for plotting
cli_use <- cli[,1:5]  # Select clinical traits for correlation analysis

# Transpose the expression matrix
mydata <- as.data.frame(t(mydata)) 

# Convert clinical variables to numeric
cli_use$status <- ifelse(cli_use$status == "Dead", 1, 0)  # Dead = 1, Alive = 0
cli_use$time <- as.numeric(cli_use$time)  # Convert to numeric
cli_use$gender <- ifelse(cli_use$gender == "male", 1, 0)  # Male = 1, Female = 0
cli_use$age <- as.numeric(cli_use$age)  # Convert to numeric

cli_use$stage_copy <- ifelse(cli_use$stage_copy == "Stage I", 1, 
                        ifelse(cli_use$stage_copy == "Stage II", 2,
                               ifelse(cli_use$stage_copy == "Stage III", 3, 
                                      ifelse(cli_use$stage_copy == "Stage IV", 4,
                                             ifelse(cli_use$stage_copy == "stage x", 4, NA)))))

# View the processed data
mydata[1:5, 1:5]  # Sample names as row names, gene names as column names
cli_use[1:5, 1:5]  # Sample names as row names, clinical traits as column names

# Compute correlations -------------------------------------------------------
result <- corr.test(x = mydata,  # Expression matrix for correlation
                    y = cli_use,  # Clinical data (must be numeric)
                    method = "pearson"  # Use Pearson correlation (default)
                    )

# Extract R-values and P-values
result_R <- result$r  # Correlation coefficients
result_P <- round(result$p, 3)  # Round P-values to 3 decimal places

# Prepare data for visualization ---------------------------------------------
# Transpose the data
data <- t(result_R)
P <- t(result_P)

# Replace P-values with significance symbols
P[P >= 0 & P < 0.001] <- "***"
P[P >= 0.001 & P < 0.01] <- "**"
P[P >= 0.01 & P < 0.05] <- "*"
P[P >= 0.05 & P <= 1] <- ""

# Heatmap visualization ------------------------------------------------------
pheatmap(data,
         display_numbers = P,  # Display significance symbols on the heatmap
         color = colorRampPalette(c("#128BAE", "white", "#E8870D"))(100),  # Color gradient
         border = "white",  # Border color
         main = "Heatmap",  # Title of the heatmap
         show_rownames = TRUE,  # Show row names
         show_colnames = TRUE,  # Show column names
         cexCol = 1,  # Scaling factor for column labels
         scale = 'none',  # No scaling (options: 'row', 'column', 'none')
         angle_col = "45",  # Angle of column labels
         legend = TRUE,  # Show legend
         legend_breaks = c(-0.2, 0, 0.2),  # Data range in the legend
         fontsize_row = 10,  # Font size for row labels
         fontsize_col = 10)  # Font size for column labels
