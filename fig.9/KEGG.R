# Set working directory ------------------------------------------------------
setwd("D:/R/KEGG_enrichment_analysis")

# Load required R packages --------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Check if the clusterProfiler package is installed
if (!require("clusterProfiler", quietly = TRUE)) {
  # If not installed, use BiocManager to install
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
# Check if the org.Hs.eg.db package is installed
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  # If not installed, use BiocManager to install
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
# Check if the ggplot2 package is installed
if (!require("ggplot2", quietly = TRUE)) {
  # If not installed, use install.packages to install
  install.packages("ggplot2")
}
# Check if the dplyr package is installed
if (!require("dplyr", quietly = TRUE)) {
  # If not installed, use install.packages to install
  install.packages("dplyr")
}
# Check if the ggsci package is installed
if (!require("ggsci", quietly = TRUE)) {
  # If not installed, use install.packages to install
  install.packages("ggsci")
}
# Load required R packages --------------------------------------------------
library(clusterProfiler) # Main package for enrichment analysis
library(org.Hs.eg.db)    # For species annotation information
library(ggplot2)         # For plotting
library(dplyr)           # For data manipulation

# Read data: results of differential analysis --------------------------------
upGene <- read.csv("upGene_1_0.05.csv", row.names = 1)
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", # Input data type
              toType = c("ENTREZID"), # Target data type
              OrgDb = org.Hs.eg.db) # Species

# KEGG enrichment analysis ---------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", # KEGG database
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

# Bar plot of KEGG enrichment results
barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        title = "KEGG_enrichment") # Title

# Dot plot of KEGG enrichment results
dotplot(KEGG)

# Issue: The KEGG database API has been updated, and Y's clusterProfiler package needs to be updated ----------------------
# Solution: Uninstall the clusterProfiler package used for enrichment analysis, then install the latest version
.libPaths() 
remove.packages("clusterProfiler", lib = file.path("your_R_package_installation_path"))
remove.packages("BiocManager", lib = file.path("your_R_package_installation_path"))
# Reinstall the necessary R packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("clusterProfiler")

# If it still doesn't work, your R version may be too old
install.packages("installr")
library(installr)
updateR()
# During the update, remember to select to transfer the old version of R packages to the new version

# Beautified version -----------------------------------------------------------
df <- as.data.frame(KEGG)

# Bubble plot ---------------------------------------------------------------
# Sort and select the top 10 most significant entries
df_sorted <- df[order(df$p.adjust), ][1:10, ]
colnames(df_sorted)
df_sorted$Description
df_sorted$Description <- gsub(" - Homo sapiens \\(human\\)", "", df_sorted$Description)

# Generate scatter plot
ggplot(df_sorted, aes(x = GeneRatio, y = Description, size = -log10(p.adjust), color = p.adjust)) +
  geom_point(alpha = 0.7) + # Add scatter plot layer with transparency
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "The most enrichment KEGG pathway", x = "GeneRatio", y = "KEGG pathway", size = "-log10(P-value)", color = "P-value") + 
  theme_minimal() + # Set the theme to minimal
  theme(axis.text.y = element_text(face = 'bold', size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold")) # Adjust text sizes
