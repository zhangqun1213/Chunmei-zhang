setwd("D:/R/GSEA分析")


# 加载所需要的R包 ----------------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

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

# 制作genelist --------------------------------------------------------------
deg<-read.csv("deg.csv",header = T,sep = ",")

names(deg)[1] <- "SYMBOL"
deg1 <- bitr(geneID = deg$SYMBOL,fromType = "SYMBOL",
             toType = "ENTREZID",OrgDb = org.Hs.eg.db)
DEG <- merge(deg,deg1,by = "SYMBOL")

DEG1<-DEG[,c(1,2,9)]
colnames(DEG1)
dega<-DEG1[,c(3,2)]
geneList = dega[,2]
names(geneList) = as.character(dega[,1])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
de <- names(geneList)[abs(geneList) > 2]#abs计算绝对值

# GO富集分析进阶版 -----------------------------------------------------------------
Genes <- bitr(deg$SYMBOL,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "MF",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 1,   #pvalue的阈值
               qvalueCutoff = 1, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol

# ###诱导 GO DAG 图 ----------------------------------------------------------
View(as.data.frame(GO))#以数据框的形式展示GO的结果
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="MF", readable=TRUE)#大家也可以试试BP、CC等
goplot(ego)

GO_MF <- simplify(GO, cutoff=0.05, by="pvalue", select_fun=min)
plotGOgraph(GO_MF,
            firstSigNodes = 10,
            useInfo = "all",
            sigForAll = TRUE,
            useFullNames = TRUE)
#注释信息：基础版GO富集分析只展示富集最显著的条目，这里我们挖掘GO term上下层级
           #关系以及富集程度的DAG图
           #→ 代表GO term上下层级关系
           #⚪  代表富集程度未在前10的条目
           #方框 代表富集程度在前10的条目 颜色代表富集程度

#要展示数据集以及数据集中包含的元素，数据集多时可以用upset plot
upsetplot(GO)

#基因概念网络--我认为最好玩的！！！！
pdf(file="基因概念网络1.pdf",width = 10,height = 10)
cnetplot(GO, categorySize="pvalue", foldChange=geneList)
dev.off()

pdf(file="基因概念网络2.pdf",width = 13,height = 8)
cnetplot(GO, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
dev.off()

pdf(file="基因概念网络3.pdf",width = 13,height = 8)
p1 <- cnetplot(GO, node_label="category") 
p2 <- cnetplot(GO, node_label="gene") 
p3 <- cnetplot(GO, node_label="all") 
p4 <- cnetplot(GO, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
dev.off()
#热图
pdf(file="热图1.pdf",width = 13,height = 8)
heatplot(GO, showCategory = 30,foldChange=geneList)
dev.off()

# KEGG通路展示进阶版 -------------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

View(as.data.frame(KEGG))#以数据框的形式展示KEGG的结果

hsa04115 <- pathview(gene.data = geneList,
                     pathway.id = "hsa04115", #上述结果中的hsa04750通路
                     species = "hsa",
                     limit = list(gene=max(abs(geneList)), cpd=1))

# GSEA分析 ------------------------------------------------------------------
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
af=as.data.frame(kk2@result)
write.table(af,file=paste0("all_GSEA.xls"),sep="\t",quote=F,col.names=T)

# 结果可视化--分别取GSEA结果的前5个后5个展示 -----------------------------------------------
View(as.data.frame(kk2@result))
num=5
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])

gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
dev.off()


# 单独展示某一个条目 ---------------------------------------------------------------
gseaplot2(kk2,
          title = "DNA replication",  #设置标题
          "hsa03030", #绘制hsa04658通路的结果，通路名称与编号对应
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, 
          pvalue_table = T) # 显示p值

# GSEA 结果表达分布的脊线图 ---------------------------------------------------------
#跑分和预排序是可视化 GSEA 结果的传统方法,可视化基因集的分布和富集分数
ridgeplot(kk2,
          showCategory = 20,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 32)
gseaplot(kk2, geneSetID = 1, by = "all", title = kk2$Description[1])

# gsego分析 -----------------------------------------------------------------
gsea_go <- gseGO(geneList     = geneList,#根据LogFC排序后的基因列表
               OrgDb        = org.Hs.eg.db,
               ont          = "ALL",#GO分析的模块
               minGSSize    = 10,#最小基因集的基因数
               maxGSSize    = 500,#最大基因集的基因数
               pvalueCutoff = 0.05,#p值的阈值
               verbose      = FALSE)#是否输出提示信息
View(as.data.frame(gsea_go))
#绘图
gseaplot2(
  gsea_go,
  geneSetID="GO:0000070",
  title = "GSE_GO",
  color = "blue",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = TRUE,
  ES_geom = "line")
