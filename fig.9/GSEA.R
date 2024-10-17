# Set working directory ------------------------------------------------------
setwd("D:/R/GSEA_Analysis")

# Load required R packages --------------------------------------------------
# Uncomment the following lines if packages are not installed
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("topGO")
#BiocManager::install("DOSE")
#BiocManager::install("pathview")

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(pathview)
library(topGO)

# Create genelist -------------------------------------------------------------
deg <- read.csv("deg.csv", header = TRUE, sep = ",")

names(deg)[1] <- "SYMBOL"
deg1 <- bitr(geneID = deg$SYMBOL, 
             fromType = "SYMBOL", 
             toType = "ENTREZID", 
             OrgDb = org.Hs.eg.db)
DEG <- merge(deg, deg1, by = "SYMBOL")

DEG1 <- DEG[, c(1, 2, 9)]
colnames(DEG1)
dega <- DEG1[, c(3, 2)]
geneList <- dega[, 2]
names(geneList) <- as.character(dega[, 1])
head(geneList)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)
de <- names(geneList)[abs(geneList) > 2] # Filter by absolute value

# Advanced GO enrichment analysis --------------------------------------------
Genes <- bitr(deg$SYMBOL,  
              fromType = "SYMBOL", 
              toType = c("ENTREZID"), 
              OrgDb = org.Hs.eg.db)

GO <- enrichGO(gene = Genes$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               keyType = "ENTREZID", 
               ont = "MF", # Options: BP/CC/MF 
               pAdjustMethod = "BH", 
               pvalueCutoff = 1, 
               qvalueCutoff = 1, 
               minGSSize = 5, 
               maxGSSize = 5000, 
               readable = TRUE)

# Generate GO DAG plot -------------------------------------------------------
View(as.data.frame(GO))
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont = "MF", readable = TRUE)
goplot(ego)

GO_MF <- simplify(GO, cutoff = 0.05, by = "pvalue", select_fun = min)
plotGOgraph(GO_MF, 
            firstSigNodes = 10, 
            useInfo = "all", 
            sigForAll = TRUE, 
            useFullNames = TRUE)

# Use upset plot for dataset visualization
upsetplot(GO)

# Gene-concept network visualization -----------------------------------------
pdf(file = "Gene_Concept_Network_1.pdf", width = 10, height = 10)
cnetplot(GO, categorySize = "pvalue", foldChange = geneList)
dev.off()

pdf(file = "Gene_Concept_Network_2.pdf", width = 13, height = 8)
cnetplot(GO, foldChange = geneList, circular = TRUE, colorEdge = TRUE)
dev.off()

pdf(file = "Gene_Concept_Network_3.pdf", width = 13, height = 8)
p1 <- cnetplot(GO, node_label = "category")
p2 <- cnetplot(GO, node_label = "gene")
p3 <- cnetplot(GO, node_label = "all")
p4 <- cnetplot(GO, node_label = "none")
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = LETTERS[1:4])
dev.off()

# Heatmap visualization ------------------------------------------------------
pdf(file = "Heatmap_1.pdf", width = 13, height = 8)
heatplot(GO, showCategory = 30, foldChange = geneList)
dev.off()

# Advanced KEGG pathway visualization -----------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID, 
                   organism = "hsa", 
                   keyType = "kegg", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1)

View(as.data.frame(KEGG))

hsa04115 <- pathview(gene.data = geneList, 
                     pathway.id = "hsa04115", 
                     species = "hsa", 
                     limit = list(gene = max(abs(geneList)), cpd = 1))

# GSEA analysis --------------------------------------------------------------
kk2 <- gseKEGG(geneList = geneList, 
               organism = 'hsa', 
               keyType = "kegg", 
               exponent = 1, 
               minGSSize = 10, 
               maxGSSize = 500, 
               pvalueCutoff = 0.05, 
               pAdjustMethod = "none", 
               by = "fgsea")

class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af <- as.data.frame(kk2@result)
write.table(af, file = "all_GSEA.xls", sep = "\t", quote = FALSE, col.names = TRUE)

# Visualization of top GSEA results -------------------------------------------
View(as.data.frame(kk2@result))
num <- 5
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore), num)])
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore), num)])
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore), num), 
                                                  tail(order(kk2@result$enrichmentScore), num))])

# Visualize a specific pathway ------------------------------------------------
gseaplot2(kk2, 
          title = "DNA replication", 
          "hsa03030", 
          color = "red", 
          base_size = 20, 
          subplots = 1:3, 
          pvalue_table = TRUE)

# Ridge plot of GSEA results --------------------------------------------------
ridgeplot(kk2, 
          showCategory = 20, 
          fill = "p.adjust", 
          core_enrichment = TRUE, 
          label_format = 32)
gseaplot(kk2, geneSetID = 1, by = "all", title = kk2$Description[1])

# gseGO analysis --------------------------------------------------------------
gsea_go <- gseGO(geneList = geneList, 
                 OrgDb = org.Hs.eg.db, 
                 ont = "ALL", 
                 minGSSize = 10, 
                 maxGSSize = 500, 
                 pvalueCutoff = 0.05, 
                 verbose = FALSE)

View(as.data.frame(gsea_go))

# Plot gseGO results ---------------------------------------------------------
gseaplot2(gsea_go, 
          geneSetID = "GO:0000070", 
          title = "GSE_GO", 
          color = "blue", 
          base_size = 11, 
          rel_heights = c(1.5, 0.5, 1), 
          subplots = 1:3, 
          pvalue_table = TRUE, 
          ES_geom = "line")
