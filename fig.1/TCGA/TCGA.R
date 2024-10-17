setwd("D:/R/TCGA_Data_Download")  # Set working directory

# Load required R packages ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

# Download data from GDC TCGA ------------------------------------------------
getGDCprojects()$project_id  # View available tumor datasets
project <- "TCGA-ESCA"
TCGA_ESCA <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query = TCGA_ESCA, method = "api")

expr <- GDCprepare(query = TCGA_ESCA)

names(expr@assays)
# For differential analysis, we need "unstranded" data, which is the count data.
# For count-based differential expression analysis, DESeq2 and edgeR are recommended. 
# You can use one or both to find common genes.

# Extract specific data based on needs ---------------------------------------
counts <- as.data.frame(assay(expr))  # Extract counts data
TPM <- as.data.frame(assay(expr, i = "tpm_unstrand"))  # Extract TPM data
# fpkm <- as.data.frame(assay(expr, i = "fpkm_unstrand"))  # Extract FPKM data

# Get additional information
data <- as.data.frame(rowRanges(expr))
colnames(data)
mydata <- data[, c("gene_type", "gene_name")]
table(mydata$gene_type)
expr_count <- cbind(gene_type = data$gene_type, gene_name = data$gene_name, counts)
expr_TPM <- cbind(gene_type = data$gene_type, gene_name = data$gene_name, TPM)

# Download clinical data -----------------------------------------------------
project <- "TCGA-ESCA"
query <- GDCquery(project = project,
                  data.category = "Clinical",
                  data.format = "bcr xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

# Save data for future use ---------------------------------------------------
save(expr_count, expr_TPM, clinical, file = "TCGA_ESCA.Rdata")

# Load saved data ------------------------------------------------------------
load("TCGA_ESCA.Rdata")

# Load additional packages ---------------------------------------------------
library(dplyr)  # For data manipulation

# Filter for protein-coding mRNAs --------------------------------------------
table(expr_count$gene_type)
mRNA <- expr_count[expr_count$gene_type == "protein_coding", ]

# Clean count data -----------------------------------------------------------
mRNA <- mRNA[, -1]
dim(mRNA)  # Check the number of rows
mRNA$mean <- rowMeans(mRNA[, 2:175])  # Calculate the mean expression value for each row
mRNA <- arrange(mRNA, desc(mean))  # Sort by mean value in descending order
mRNA <- mRNA %>% distinct(gene_name, .keep_all = TRUE)  # Remove duplicate gene names
mRNA <- select(mRNA, -c(mean))  # Remove the 'mean' column

# Separate tumor and normal samples ------------------------------------------
rownames(mRNA) <- mRNA$gene_name
mRNA <- mRNA[, -1]
table(substr(colnames(mRNA), 14, 16))

Tumor <- grep('01A', colnames(mRNA))
Tumor  # Position of tumor samples
Tumor_mRNA <- mRNA[, Tumor]

Normal <- grep('11A', colnames(mRNA))
Normal  # Position of normal samples
Normal_mRNA <- mRNA[, Normal]

mRNA_count <- cbind(Normal_mRNA, Tumor_mRNA)
# write.csv(mRNA_count, file = "mRNA_count.csv")

# Clean TPM data -------------------------------------------------------------
mRNA_TPM <- expr_TPM[expr_TPM$gene_type == "protein_coding", ]
mRNA_TPM <- mRNA_TPM[, -1]
dim(mRNA_TPM)  # Check the number of rows
mRNA_TPM$mean <- rowMeans(mRNA_TPM[, 2:175])  # Calculate the mean expression value for each row
mRNA_TPM <- arrange(mRNA_TPM, desc(mean))  # Sort by mean value in descending order
mRNA_TPM <- mRNA_TPM %>% distinct(gene_name, .keep_all = TRUE)  # Remove duplicate gene names
mRNA_TPM <- select(mRNA_TPM, -c(mean))  # Remove the 'mean' column

# Separate tumor and normal samples ------------------------------------------
rownames(mRNA_TPM) <- mRNA_TPM$gene_name
mRNA_TPM <- mRNA_TPM[, -1]
table(substr(colnames(mRNA_TPM), 14, 16))

Tumor <- grep('01A', colnames(mRNA_TPM))
Tumor  # Position of tumor samples
Tumor_mRNA_TPM <- mRNA_TPM[, Tumor]

Normal <- grep('11A', colnames(mRNA))
Normal  # Position of normal samples
Normal_mRNA_TPM <- mRNA_TPM[, Normal]

mRNA_TPM_new <- cbind(Normal_mRNA_TPM, Tumor_mRNA_TPM)

# Save cleaned data for future analysis --------------------------------------
save(mRNA_count, mRNA_TPM_new, file = "count_and_TPM.Rdata")
