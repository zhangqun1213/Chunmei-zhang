# Set the working directory --------------------------------------------------
setwd("D:/bilibiliR/38_IOBR_Immune_Infiltration_Part1")

# Code to install the IOBR package -------------------------------------------
# options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

depens <- c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2", "devtools",
            'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
            'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", 
            "timeROC", "pracma")

for (i in 1:length(depens)) {
  depen <- depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(version = "3.17")(depen, update = FALSE)
}

if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")

# Load the IOBR package ------------------------------------------------------
library(IOBR)

# Load other required R packages ---------------------------------------------
library(tidyverse)

# Load example data ----------------------------------------------------------
load("TCGA-ESCA_mrna_expr_tpm.rdata")  # Load standardized RNA-seq data

# Backup the data
mydata <- mrna_expr_tpm

# Explore Tumor Microenvironment (TME) Deconvolution Algorithms --------------
tme_deconvolution_methods

# View available options for signature score calculation ---------------------
signature_score_calculation_methods
# Methods include ssGSEA, PCA, and Z-score.

# Explore TME-related gene sets ----------------------------------------------
names(signature_tme)

# Explore metabolism-related gene sets ---------------------------------------
names(signature_metabolism)

# Explore gene sets related to biomedical research (e.g., m6A, exosomes) -----
names(signature_tumor)

# Explore immune cell signature gene sets ------------------------------------
signature_collection

# View citations for the reference gene sets ---------------------------------
signature_collection_citation

# Part I: Deconvolution Algorithms

# CIBERSORT immune cell infiltration analysis -------------------------------
cibersort <- deconvo_tme(eset = mydata,  # Gene expression matrix with genes as rows and samples as columns
                         method = "cibersort",  # Use 'CIBERSORT' method
                         arrays = FALSE,  # Data is not from microarrays
                         perm = 100)  # Permutations for statistical analysis (recommended ≥100)

cell_bar_plot(input = cibersort[1:20,],  # Plot results for the first 20 samples
              title = "Cell Fraction", 
              legend.position = "bottom",  # Set legend position
              palette = 3,  # Select color palette
              coord_filp = TRUE,  # Flip the coordinate axes
              show_col = FALSE)  # Hide cell type colors

# Save the plot as a PDF
ggsave("Cibersort.pdf", width = 10, height = 8)

# xCell immune cell infiltration analysis ------------------------------------
xcell <- deconvo_tme(eset = mydata,  # Gene expression matrix
                     method = "xcell",  # Use 'xcell' method
                     arrays = FALSE)  # Data is not from microarrays
head(xcell)

# EPIC immune cell infiltration analysis -------------------------------------
epic <- deconvo_tme(eset = mydata,  # Gene expression matrix
                    method = "epic",  # Use 'epic' method
                    arrays = FALSE)  # Data is not from microarrays
head(epic)

cell_bar_plot(input = epic[1:20,],  # Plot results for the first 20 samples
              title = "Cell Fraction", 
              legend.position = "bottom", 
              palette = 3, 
              coord_filp = TRUE, 
              show_col = FALSE)

# Save the plot as a PDF
ggsave("epic.pdf", width = 10, height = 8)

# Other methods available for exploration:
# MCPcounter               EPIC              xCell 
# "mcpcounter"             "epic"            "xcell" 
# CIBERSORT        CIBERSORT Absolute          IPS 
# "cibersort"         "cibersort_abs"         "ips" 
# ESTIMATE                 SVR               lsei 
# "estimate"              "svr"             "lsei" 
# TIMER                quanTIseq 
# "timer"             "quantiseq" 

# Divider --------------------------------------------------------------
# Function usage guide (no need to run this section) --------------------------
deconvo_tme(eset,  # Gene expression matrix with genes as rows and samples as columns
            project = NULL,  # Project name to distinguish datasets (optional)
            method = tme_deconvolution_methods,  # Choose deconvolution method
            arrays = FALSE,  # Set to TRUE for microarray data
            tumor = TRUE,  # Use tumor-specific matrices/processes
            perm = 1000,  # Permutations for statistical analysis (recommended ≥100)
            reference,  # Immune cell reference matrix (e.g., lm22, lm6)
            scale_reference,  # Scale reference data if TRUE
            plot = FALSE,  # Whether to plot the results
            scale_mrna,  # Adjust mRNA content for different cell types
            group_list = NULL,  # List of sample tumor types
            platform = "affymetrix",  # Platform type, e.g., "affymetrix"
            absolute.mode = FALSE,  # Run in absolute mode for CIBERSORT or SVR
            abs.method = "sig.score")  # Choose 'no.sumto1' or 'sig.score' for absolute mode
