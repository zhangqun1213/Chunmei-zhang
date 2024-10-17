setwd("D:/R/GEO_Data_Download")  # Set working directory

# Load required R packages ---------------------------------------------------
# install.packages("devtools")
library(devtools)
# install.packages("GEOquery")
# install_github("jmzeng1314/GEOmirror")
# install_github("jmzeng1314/idmap1")

library(GEOmirror)
library(idmap1)
library(GEOquery)
options('download.file.method.GEOquery' = 'libcurl')
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(dplyr)
library(limma)
library(tidyr)

# Alternative download method ------------------------------------------------
install_github("jmzeng1314/AnnoProbe")
library(AnnoProbe)

# Download GEO data and annotation information -------------------------------
eset <- geoChina('GSE53819')  # Method 1: Use domestic mirror (requires GEOmirror package)
eset <- getGEO('GSE53819', destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)  # Method 2: Recommended (requires GEOquery package)
ann_info <- getGEO(GEO = 'GPL6480', destdir = ".")

# Get the expression matrix --------------------------------------------------
eset <- eset[[1]]
probes_expr <- exprs(eset)
dim(probes_expr)
data_use <- as.data.frame(probes_expr)

# Get clinical information ---------------------------------------------------
phenoDat <- pData(eset)
dim(phenoDat)
Samples <- row.names(phenoDat)
phenoDat <- cbind(Samples, phenoDat)
colnames(phenoDat)
group <- phenoDat[, c("Samples", "source_name_ch1")] 
rownames(group) <- NULL
group$source_name_ch1 <- c(rep("Tumor", 18), rep("Normal", 18))
colnames(group) <- c("Samples", "group")
group
write.csv(phenoDat, "clinicaldata.csv")

# Check if log transformation is needed --------------------------------------
boxplot(data_use, las = 2)  # Display row names completely with 'las'
qx <- as.numeric(quantile(data_use, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))  
# Get data distribution and sample quantiles
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)  # Determine if log transformation is needed
LogC
if (LogC) { 
  data_use[which(data_use <= 0)] <- NaN
  data_use <- log2(data_use) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}  # Perform log2 transformation if needed

# Data inspection: Boxplot, PCA plot, and Hierarchical Clustering -------------

# 1. Boxplot
boxplot(data_use, las = 2)  # Display row names completely

# Standardize if the means are inconsistent ----------------------------------
data_use <- normalizeBetweenArrays(data_use)
boxplot(data_use, outline = FALSE, notch = TRUE, las = 2)

# 2. PCA plot
PCA_result <- PCA(t(data_use), graph = FALSE)
fviz_pca_ind(PCA_result,
             geom.ind = c("point", "text"),  # Show points and sample names
             mean.point = FALSE,  # Do not display the center point
             repel = TRUE,  # Adjust sample names to avoid overlap
             col.ind = group$group,  # Assign different colors by group
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE,  # Add confidence ellipses
             legend.title = "Groups")

# 3. Hierarchical Clustering Plot
Sample_clust <- dist(t(data_use))  # Calculate distance between variables
hc <- hclust(Sample_clust)  # Perform hierarchical clustering
plot(hc, hang = -1, cex = 0.8)

# Get annotation information -------------------------------------------------
Meta(ann_info)$title
anno_geoquery <- Table(ann_info)

# Convert probe IDs to gene symbols ------------------------------------------
colnames(anno_geoquery)
id_symbol <- anno_geoquery %>%
  select(ID, GENE_SYMBOL) %>%
  separate(GENE_SYMBOL, c("GENE_SYMBOL", NA), sep = "///")  # Split gene names
id_symbol[, 2] <- trimws(id_symbol[, 2])  # Remove whitespace from the 'gene_symbol' column

# Add gene names to the expression matrix
data_use$probe_id <- rownames(data_use)
data_with_name <- merge(data_use, id_symbol, by.x = "probe_id", by.y = "ID")  
# Use merge function to add gene names
dim(data_with_name)
data_with_name <- data_with_name[data_with_name$GENE_SYMBOL != "", ]  # Remove rows with empty gene symbols
dim(data_with_name)
table(duplicated(data_with_name[, ncol(data_with_name)]))

# Take the average value for duplicate gene symbols
data_with_name <- avereps(data_with_name[, -c(1, ncol(data_with_name))], 
                          ID = data_with_name$GENE_SYMBOL)  
# Check if there are still duplicate gene symbols
table(duplicated(rownames(data_with_name)))

# Check the expression levels for internal controls
data_with_name['GAPDH', ]  
data_with_name['ACTB', ]

# Save the expression matrix and group data ----------------------------------
write.table(data.frame(gene_symbol = rownames(data_with_name), data_with_name),
            file = "Matrix_of_expression.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.csv(group, "group.csv")
