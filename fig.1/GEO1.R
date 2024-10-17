setwd("D:/R/GEO_Data_Merging")  # Set working directory

# Problem to solve -----------------------------------------------------------
# Problem 1: Merge multiple GEO datasets (eliminate batch effect)
# Problem 2: Read downloaded data and perform ID conversion
# Problem 3: Handle network issues where annotation data download fails 
#            (use AnnoProbe for annotation of genomic data)

# Load required R packages ---------------------------------------------------
library(devtools)  # Install R packages from GitHub
library(GEOmirror)  # GEO mirror
library(GEOquery)  # Download GEO data
#---------------
# install_github("jmzeng1314/AnnoProbe")
library(AnnoProbe)  # Annotate genomic data
library(limma)  # Normalization
library(sva)  # Remove batch effects

# The 'sva' package explores and removes batch effects in gene expression data.
# Batch effects refer to variables unrelated to the primary research interest but 
# still present in the data, possibly due to differences in sample source or collection time.

# Define the color palette
mycolors <- c('#54D9ED','#324B4F','#95B0B5','#D1BCFE','#FCFCD4','#F78ABA',
              '#00A9CB','#8EF8B6','#ED546C','#1BBEB9','#DE90DE','#00A44F')
              
# Download Dataset 1 ---------------------------------------------------------
# Read annotation files downloaded manually
expr_data <- read.table("GSE31821_series_matrix.txt", sep = "\t", 
                        header = TRUE, fill = TRUE, row.names = 1, 
                        comment.char = '!')  
# Lines starting with '!' are treated as comments
Anno_Info <- read.table("GPL570-55999.txt", sep = "\t", quote = "", 
                        header = TRUE, fill = TRUE)

# Simple data arrangement
colnames(Anno_Info)
Anno_Info <- Anno_Info[, c("ID", "Gene.Symbol", "ENTREZ_GENE_ID")]  # Select required columns

eset_1 <- expr_data
boxplot(eset_1, outline = FALSE, notch = FALSE, las = 2, col = '#54D9ED')  # Simple boxplot
dim(eset_1)  # Check number of samples and genes

# If data is not uniform
# eset_1 <- normalizeBetweenArrays(eset_1)

# Data cleaning and ID conversion
same1 <- match(rownames(eset_1), Anno_Info$ID)  # Find matching elements
eset_1$Gene_Symbol <- Anno_Info[same1, "Gene.Symbol"]  # Add Gene_Symbol to expression matrix
eset_1 <- eset_1[!duplicated(eset_1$Gene_Symbol), ]  # Keep the first occurrence and remove duplicates

# For GEO data, we take the mean of duplicate gene symbols;
# For TCGA data, we take the maximum expression value;
# Here, we directly remove duplicate genes.
eset_1 <- eset_1[!is.na(eset_1$Gene_Symbol) & !grepl('//', eset_1$Gene_Symbol), ]
eset_1 <- eset_1[eset_1$Gene_Symbol != "", ]  # Remove empty values
rownames(eset_1) <- eset_1$Gene_Symbol  # Modify row names
eset_1 <- eset_1[, -7]  # Remove unnecessary columns (replace with the last column in your matrix)

save(eset_1, file = "GSE31821.Rdata")

# Download Dataset 2 ---------------------------------------------------------
eset_2_raw <- geoChina('GSE66724')
Anno_Info_2 <- idmap('GPL570', type = 'soft')

# Extract expression matrix
eset_2_raw <- eset_2_raw[[1]]
eset_2 <- exprs(eset_2_raw)  # Extract expression matrix
eset_2 <- as.data.frame(eset_2)
dim(eset_2)

boxplot(eset_2, outline = FALSE, notch = FALSE, las = 2, col = '#E93639')  # Simple boxplot

# If data is not uniform
# eset_2 <- normalizeBetweenArrays(eset_2)

# Data cleaning and ID conversion
same2 <- match(rownames(eset_2), Anno_Info_2$ID)  # Find matching elements
eset_2$Gene_Symbol <- Anno_Info_2[same2, "symbol"]  # Add Gene_Symbol to expression matrix
eset_2 <- eset_2[!duplicated(eset_2$Gene_Symbol), ]  # Keep the first occurrence and remove duplicates

eset_2 <- eset_2[!is.na(eset_2$Gene_Symbol) & !grepl('//', eset_2$Gene_Symbol), ]
eset_2 <- eset_2[eset_2$Gene_Symbol != "", ]  # Remove empty values
rownames(eset_2) <- eset_2$Gene_Symbol  # Modify row names
eset_2 <- eset_2[, -17]  # Remove unnecessary columns

# Data merging ---------------------------------------------------------------
identical(eset_1, eset_2)  # Check if the two datasets are identical (they are not)
eset_2 <- eset_2 + 2  # Adjust for demonstration purposes

same_genes <- intersect(rownames(eset_1), rownames(eset_2))  # Find common genes
eset_merge <- cbind(eset_1[same_genes, ], eset_2[same_genes, ])  # Merge by common genes

col <- c(rep("#54D9ED", 6), rep("#E93639", 16))  # Blue for first 6 samples, red for the next 16
par(mar = c(8, 3, 1, 2) + 0.1)  # Adjust margins to display X-axis labels
boxplot(eset_merge, col = col, las = 2)  # Overall boxplot

# Read group data ------------------------------------------------------------
group_1 <- read.csv("group.csv", header = TRUE)  # Manually created group data

phenoDat <- pData(eset_2_raw)  # Extract clinical data from GSE66724
colnames(phenoDat)  # View available clinical information
group_2 <- phenoDat[, c("geo_accession", "treatment_protocol_ch1")]  # Extract relevant columns
colnames(group_2) <- c("Sample", "group")  # Rename columns
group_2$group <- "atrial fibrillation"  # Assign group label

# Merge group data
group <- rbind(group_1, group_2)

# Start batch effect removal -------------------------------------------------
GSE <- c(rep('GSE31821', 6), rep('GSE66724', 16))
group_list <- group$group

table(group_list, GSE)

data <- eset_merge
batch <- c(rep('GSE31821', 6), rep('GSE66724', 16))
design <- model.matrix(~group_list)

# Remove batch effect using the sva package
expr_limma <- removeBatchEffect(data, batch = batch, design = design)

boxplot(expr_limma, col = col, las = 2)  # Plot after batch effect removal

# Compare before and after batch effect removal
par(mfrow = c(1, 2))
boxplot(eset_merge, col = col, las = 2, main = "Before")
boxplot(expr_limma, col = col, las = 2, main = "After")
