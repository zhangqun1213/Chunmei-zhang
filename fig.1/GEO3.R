setwd("D:/R/GEO_Annotations")  # Set working directory

# Install required R packages ------------------------------------------------
# Check if the 'tidyverse' package is installed
if (!require("tidyverse", quietly = TRUE)) {
  # If not installed, use the install.packages function to install it
  install.packages("tidyverse")
}

# Check if the 'GEOquery' package is installed
if (!require("GEOquery", quietly = TRUE)) {
  # If not installed, install it using BiocManager
  BiocManager::install("GEOquery")
}

# Check if the 'ggVennDiagram' package is installed
if (!require("ggVennDiagram", quietly = TRUE)) {
  # If not installed, use the install.packages function to install it
  install.packages("ggVennDiagram")
}

# Load required R packages ---------------------------------------------------
library(GEOquery)  # For data download
library(tidyverse)  # For data cleaning
library(ggVennDiagram)  # For drawing Venn diagrams and visualizing intersections

# Download sample data -------------------------------------------------------
# 1. Expression matrix
# Download the GSE138043 gene expression dataset from GEO
eset <- getGEO('GSE138043', destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)  

eset <- eset[[1]]  # Select the first sample dataset
probes_expr <- exprs(eset)  # Get the gene expression matrix 'probes_expr'
dim(probes_expr)  # Check the dimensions of the gene expression matrix
data_use <- as.data.frame(probes_expr)  # Convert the gene expression matrix to a data frame 'data_use'

# 2. Annotation information
# Download platform GPL5175 annotation information from GEO
ann_info <- getGEO(GEO = 'GPL5175', destdir = ".")

anno_geoquery <- Table(ann_info)  # Use the GEOquery package to convert 'ann_info' into the dataset 'anno_geoquery'
colnames(anno_geoquery)  # View column names
test <- anno_geoquery[c(1, 11, 15), c(1, 10)]  # Extract partial information for testing

# Split annotation information (test) ----------------------------------------
# Split the 'gene_assignment' column (keep all results)
test_split_1 <- test %>% 
  separate_rows(gene_assignment, sep = " /// ")  %>%  # Split the 'gene_assignment' column by "///" into new rows, and store the results in a new data frame 'split_df'
  separate(gene_assignment, into = c("Ensembl_ID", "gene_symbol", 
                                     "gene_description", "location", "ENTREZ_ID"), 
           sep = " // ")  # Split the 'gene_assignment' column by "//" into five columns: 'Ensembl_ID', 'gene_symbol', 'gene_description', 'location', and 'ENTREZ_ID'

# Split the 'gene_assignment' column (keep only the first probe's corresponding gene symbol)
test_split_2 <- test %>% 
  separate_rows(gene_assignment, sep = " /// ")  %>%  
  separate(gene_assignment, into = c("Ensembl_ID", "gene_symbol", 
                                     "gene_description", "location", "ENTREZ_ID"), 
           sep = " // ")  %>%  
  distinct(ID, .keep_all = TRUE)  # Keep only the first occurrence of each unique ID

# Split annotation information (final version) -------------------------------
colnames(anno_geoquery)  # View column names
anno_info <- anno_geoquery[, c("ID", "gene_assignment")]  # Select only the 'ID' and 'gene_assignment' columns
anno_info <- anno_info[anno_info$gene_assignment != "---", ]  # Remove rows with '---' in 'gene_assignment'

mydata <- anno_info %>%
  separate_rows(gene_assignment, sep = " /// ") %>%  
  separate(gene_assignment, into = c("Ensembl_ID", "gene_symbol", 
                                     "gene_description", "location", "ENTREZ_ID"), 
           sep = " // ") %>%  
  distinct(ID, .keep_all = TRUE)  # Keep only the first occurrence of each unique ID

# Find the intersection of vectors -------------------------------------------
intersect_array <- intersect(rownames(data_use), mydata$ID)

# Simple visualization of a Venn diagram -------------------------------------
list_ID <- list('Expression Matrix Probes' = rownames(data_use),
                'Annotation Probes' = as.character(mydata$ID))
ggVennDiagram(list_ID, label_alpha = 0) + 
  scale_fill_gradient(low = "gray", high = "skyblue")  
# 'ggVennDiagram' generates a radial Venn diagram. It takes a list as input.
# 'label_alpha = 0': Controls the transparency of the labels, setting it to 0 hides them.
# 'scale_fill_gradient' adjusts the fill color of the Venn diagram, setting the low values to gray and high values to sky blue.
