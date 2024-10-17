# Set working directory ------------------------------------------------------
setwd("D:/R/GO_Enrichment_Circle")

# Load required R packages -------------------------------------------------
# Check if the GOplot package is installed
if (!require("GOplot", quietly = TRUE)) {
  # If not installed, use install.packages to install
  install.packages("GOplot")
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Check if the clusterProfiler package is installed
if (!require("clusterProfiler", quietly = TRUE)) {
  # If not installed, use BiocManager to install
  BiocManager::install("clusterProfiler")
}

# Check if the org.Hs.eg.db package is installed
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  # If not installed, use BiocManager to install
  BiocManager::install("org.Hs.eg.db")
}

# Load R packages -----------------------------------------------------------
library(clusterProfiler) # Main package for enrichment analysis
library(org.Hs.eg.db) # For species annotation information
library(GOplot) # For visualization

# Import differential analysis results ---------------------------------------
DEG <- read.table("DEG_result.txt", header = TRUE)

Threshold <- 0.0005 # Set threshold

DEG_use <- DEG[DEG$adj.P.Val < Threshold, ] # Filtering step

colnames(DEG_use) # View column names, basis for selecting columns below
colnames(DEG_use)[1] <- "SYMBOL" # Rename column for merging
DEG_use <- DEG_use[, c("SYMBOL", "logFC")] # Keep only gene names and log fold change columns

# Read data: results of gene differential analysis ---------------------------
Genes <- bitr(DEG_use$SYMBOL,  
              fromType = "SYMBOL", # Input data type
              toType = c("ENTREZID"), # Type of data to convert
              OrgDb = org.Hs.eg.db) # Species database

data_use <- merge(DEG_use, Genes, by = "SYMBOL") # Merge data by SYMBOL column, remove genes without matching ENTREZID

# GO enrichment analysis and organization -----------------------------------
GO <- enrichGO(gene = data_use$ENTREZID, # Input genes' ENTREZID
               OrgDb = org.Hs.eg.db, # Annotation information from org.Hs.eg.db database
               keyType = "ENTREZID", # Specify the gene identifier type to use
               ont = "all", # Type of ontology for enrichment analysis: "BP" (Biological Process), "CC" (Cell Component), "MF" (Molecular Function) or "all" (all three)
               pAdjustMethod = "BH", # p-value adjustment method for multiple tests (Benjamini-Hochberg method)
               pvalueCutoff = 1, # Threshold for p-value (genes below this threshold will be considered significant)
               qvalueCutoff = 1, # Threshold for q-value (q-value below this threshold will be considered significant)
               minGSSize = 5, # Minimum size of gene sets (number of genes)
               maxGSSize = 5000, # Maximum size of gene sets
               readable = TRUE) # Whether to convert ENTREZID to gene symbols for better readability

# Convert GO results to data frame format
GO_result <- as.data.frame(GO) 
# Select significantly enriched pathways based on given threshold (Threshold)
GO_result <- GO_result[(GO_result$pvalue < 0.05 & GO_result$p.adjust < 0.05), ] 
# Create a data frame containing the GO enrichment analysis results
go_result <- data.frame(Category = GO_result$ONTOLOGY, 
                        ID = GO_result$ID, # Pathway ID in the enrichment analysis results
                        Term = GO_result$Description, # Pathway description in the enrichment analysis results
                        Genes = gsub("/", ", ", GO_result$geneID), # Genes in the enrichment results, multiple IDs separated by commas
                        adj_pval = GO_result$p.adjust) # Adjusted p-value in the enrichment analysis results
# Create a data frame containing gene IDs and differential expression values
genelist <- data.frame(ID = data_use$SYMBOL, logFC = data_use$logFC) 
# Set the row names of the genelist data frame to ID
row.names(genelist) = genelist[, 1] 
# Generate a visualization chart for pathway enrichment analysis results
circ <- circle_dat(go_result, genelist) 
head(circ)

# Visualization preparation and plotting --------------------------------------
nrow(go_result) # Check how many pathways are there
term_num = 3 # Set the number of GO entries to display in the chart  

nrow(genelist) # Check how many genes there are
gene_num = nrow(genelist)  # Set the maximum number of genes to display in the chart 

# Generate a chord diagram for the pathway enrichment analysis results
chord <- chord_dat(circ, # A data frame containing the pathway enrichment analysis results
                   genelist[1:gene_num, ], # A data frame containing gene IDs and differential expression values, only the top gene_num genes
                   go_result$Term[1:term_num] # A vector containing pathway descriptions, only the top term_num pathways
)

# Start plotting (generate PDF for higher clarity)

# 1. GO enrichment chord diagram
pdf(file = "GO_Enrichment_Circle.pdf", width = 12, height = 12)
GOChord(chord, 
        title = "GOcircos",
        # Title
        space = 0.01,
        # Set spacing between genes to 0.01
        gene.order = 'logFC',    
        # Order genes by logFC value
        gene.size = 2,
        # Size of gene labels
        nlfc = 1,
        # Define the number of logFC columns (default = 1)
        lfc.col = c('#F5053C', '#EDF5F7', 'cyan1'),
        # Fill colors for logFC: c (low value color, midpoint color, high value color)
        lfc.min = -3,
        # Set minimum value for logFC ratio (default = -3)
        lfc.max = 3,
        # Set maximum value for logFC ratio (default = 3)
        gene.space = 0.2,       
        # Distance of genes from the circle
        border.size = 0.2,  
        # Size of the functional area border
        process.label = 12,
        # Size of the legend
        limit = c(1, 5)
        # Control the number of GO entries and genes in the chord diagram, limit excessive or insufficient relationships, ensuring clear and readable visualization.
        # The first value is 1, meaning each gene must be assigned at least one GO entry
        # The second value is 5, meaning each displayed GO entry must be associated with at least 5 genes.
)       
dev.off()

# 2. GO circular dendrogram
pdf(file = "GO_Circular_Dendrogram.pdf", width = 10, height = 8)
GOCluster(circ, 
          metric = "euclidean",
          # Set distance metric to Euclidean
          clust.by = 'logFC',
          # Cluster based on logFC values
          lfc.col = c('#F5053C', '#EDF5F7', 'cyan1'),
          # Fill colors for logFC: c (low value color, midpoint color, high value color)
          as.character(go_result[1:term_num, 3])
          # Select the top term_num GO terms from go_result
)
dev.off()

# 3. Tai Chi plot
GOCircle(circ, nsub = 8)

go_result$ID[1:20]
Select_GO_Term <- c('GO:0042113',
                     'GO:0042100',
                     'GO:0048731',
                     'GO:0030183',
                     'GO:0003341')
GOCircle(circ, nsub = Select_GO_Term)
# - The outermost shows GO entries, the outer circle represents genes (upregulated/downregulated)

# 4. Bubble plot
GOBubble(circ,
         labels = 9, # Threshold for labeling (GO:...) in the plot; terms with -log(p.adj) > 6 will be annotated, larger values show fewer displays
         display = 'multiple', # Show bubbles for multiple variables on the same plot
         colour = c('#F5053C', '#F9F871', 'cyan1'), # Set colors
         bg.col = TRUE # Whether to add background color
)
GOBubble(circ)
# The Z-score is calculated by normalizing the expression level of the gene in the pathway compared to the expression levels of all genes in the pathway.
