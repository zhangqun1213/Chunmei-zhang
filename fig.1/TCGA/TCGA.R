setwd("D:/R/TCGA数据下载")

# 加载所需要的R包 ----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

# 从GDC TCGA下载数据 -----------------------------------------------------------
getGDCprojects()$project_id #查看有哪些肿瘤数据
project="TCGA-ESCA"
TCGA_ESCA <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query = TCGA_ESCA,method="api")

expr <- GDCprepare(query=TCGA_ESCA)

names(expr@assays)
#差异分析，我们要"unstranded"的数据，就是count的数据。
#基于count数据的差异分析,推荐是DESeq2和edgeR,取其中一种，或者取交集

# 按需求提取不同的数据 --------------------------------------------------------------
counts <- as.data.frame(assay(expr)) #提取counts数据
TPM <- as.data.frame(assay(expr,i = "tpm_unstrand")) #提取TPM数据
#fpkm <- as.data.frame(assay(expr,i = "fpkm_unstrand")) #提取FPKM数据

#获取其它信息数据
data=as.data.frame(rowRanges(expr))
colnames(data)
mydata <- data[,c("gene_type","gene_name")]
table(mydata$gene_type)
expr_count <-  cbind(gene_type=data$gene_type,gene_name=data$gene_name,counts)
expr_TPM <- cbind(gene_type=data$gene_type,gene_name=data$gene_name,TPM)

# 临床数据的下载 -----------------------------------------------------------------
project = "TCGA-ESCA"
query <- GDCquery(project = project,
                  data.category = "Clinical",
                  data.format = "bcr xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

# 数据存储 --------------------------------------------------------------------
save(expr_count,expr_TPM,clinical,file = "TCGA_ESCA.Rdata")


# 数据读取 --------------------------------------------------------------------
load("TCGA_ESCA.Rdata")

# 加载所需要的包 -----------------------------------------------------------------
library(dplyr) #数据整理
#只取protein coding，编码蛋白的mRNA进行分析
table(expr_count$gene_type)
mRNA <- expr_count[expr_count$gene_type=="protein_coding",]

# count数据清洗 ---------------------------------------------------------------
mRNA <- mRNA[,-1]
dim(mRNA)#查看一共多少行
mRNA$mean <- rowMeans(mRNA[,2:175])#计算每一行表达量的平均值
mRNA <- arrange(mRNA,desc(mean))#按mean值大小降序排列
mRNA <- mRNA %>% distinct(gene_name, .keep_all = T)
mRNA <- select(mRNA,-c(mean))#删除之前生成的mean一行

#将肿瘤和非肿瘤样本区分开排序
rownames(mRNA) <- mRNA$gene_name
mRNA <- mRNA[,-1]
table(substr(colnames(mRNA),14,16))

Tumor <- grep('01A',colnames(mRNA))
Tumor  #肿瘤样本所处位置
Tumor_mRNA <- mRNA[,Tumor]

Normal <- grep('11A',colnames(mRNA))
Normal #正常样本所处位置
Normal_mRNA <- mRNA[,Normal]

mRNA_count <- cbind(Normal_mRNA,Tumor_mRNA)
#write.csv(mRNA_count,file = "mRNA_count.csv")

# TPM数据清洗 -----------------------------------------------------------------
mRNA_TPM <- expr_TPM[expr_TPM$gene_type=="protein_coding",]
mRNA_TPM <- mRNA_TPM[,-1]
dim(mRNA_TPM)#查看一共多少行
mRNA_TPM$mean <- rowMeans(mRNA_TPM[,2:175])#计算每一行表达量的平均值
mRNA_TPM <- arrange(mRNA_TPM,desc(mean))#按mean值大小降序排列
mRNA_TPM <- mRNA_TPM %>% distinct(gene_name, .keep_all = T)
mRNA_TPM <- select(mRNA_TPM,-c(mean))#删除之前生成的mean一行
#将肿瘤样本和非肿瘤样本分开排序
rownames(mRNA_TPM) <- mRNA_TPM$gene_name
mRNA_TPM <- mRNA_TPM[,-1]
table(substr(colnames(mRNA_TPM),14,16))

Tumor <- grep('01A',colnames(mRNA_TPM))
Tumor #肿瘤样本所处位置
Tumor_mRNA_TPM <- mRNA_TPM[,Tumor]

Normal <- grep('11A',colnames(mRNA))
Normal #正常样本所处位置
Normal_mRNA_TPM <- mRNA_TPM[,Normal]

mRNA_TPM_new <- cbind(Normal_mRNA_TPM,Tumor_mRNA_TPM)

# 保存文件用于下次的分析 -------------------------------------------------------------
save(mRNA_count,mRNA_TPM_new,file = "count_and_TPM.Rdata")