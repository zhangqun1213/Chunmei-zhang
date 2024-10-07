setwd("D:/R/DESEq2差异表达")

# 清空变量 --------------------------------------------------------------------
rm(list=ls(all=TRUE))
options(stringsAsFactors = F)
library(DESeq2)
# 读取表达矩阵 ------------------------------------------------------------------
load("count_and_TPM.Rdata")

expr <- mRNA_count
expr <- expr[rowMeans(expr) > 1,]

# 创建分组信息 ------------------------------------------------------------------
table(substr(colnames(expr),14,16))#查看肿瘤样本和正常样本数量
Tumor <- grep('01A',colnames(expr))#查看肿瘤样本所处位置
Tumor
group <- factor(c(rep("Normal",times=11),rep("Tumor",times=152)))#创建分组因子变量

Data <- data.frame(row.names = colnames(expr), #创建分组数据框
                   group = group)
View(Data)
#基因表达矩阵(expr)、样本分组数据框(Data)都已经准备完毕


# 开始进行差异表达分析 --------------------------------------------------------------
#第一步：构建DEseq2对象(dds)
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = Data,
                              design = ~ group)
#第二步：开始差异分析
dds2 <- DESeq(dds)
res <- results(dds2, contrast=c("group", "Tumor", "Normal"))#肿瘤在前，对照在后
##或者res= results(dds)
res <- res[order(res$pvalue),]
summary(res)
my_result <- as.data.frame(res)#转成容易查看的数据框
my_result <- na.omit(my_result)#删除倍数为0的值
#第三步：保存差异分析的结果
library(dplyr)
my_result$Gene_symbol<-rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
rownames(my_result) <- NULL
write.csv(my_result,file="my_result_DESeq2.csv")


# DEG的筛选 ------------------------------------------------------------------
my_result$regulate <- ifelse(my_result$padj > 0.05, "unchanged",
                       ifelse(my_result$log2FoldChange > 2, "up-regulated",
                              ifelse(my_result$log2FoldChange < -2, "down-regulated", "unchanged")))
table(my_result$regulate)
#可以把上调基因和下调基因取出放在一块
DEG_deseq2 <-subset(my_result, padj < 0.05 & abs(log2FoldChange) > 2) 
upgene <- DEG_deseq2[DEG_deseq2$regulate=='up-regulated',]
downgene <- DEG_deseq2[DEG_deseq2$regulate=='down-regulated',]
write.csv(DEG_deseq2,file= "DEG_deseq2.csv")

# 结果可视化 -------------------------------------------------------------------
plot(my_result$log2FoldChange,-log2(my_result$padj))#基础函数查看火山图
library(ggplot2)
library(ggrepel)
p <- ggplot(data=my_result, aes(x=log2FoldChange, y=-log10(padj),color=regulate)) + 
  geom_point(shape = 16, size=2) + 
  theme_set(theme_set(theme_bw(base_size=20))) + 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(plot.title = element_text(size=15,hjust = 2.5)) + 
  theme_classic()+
  scale_colour_manual(values = c('#86D47D','#DFE0DF','#C34B99'))+
  geom_vline(xintercept = c(-2,2),lty=4,col ="gray",lwd=0.8)+ #logFC分界线
  geom_hline(yintercept=-log10(0.05),lty=2,col = "gray",lwd=0.6)+#adj.p.val分界线
  #annotate("text",label="886 genes are up-regulated",x = 5.5,y = 30)+    #可以去掉
  #annotate("text",label="623 genes are down-regulated",x = 5.7,y = 29)+   #可以去掉
  labs(title='TCGA-ESCA')+#图片左上角的名称
  annotate("text",x=upgene$log2FoldChange[1:3],y=(-log10(upgene$padj[1:3])),label=upgene$Gene_symbol[1:3], size=5.0)+
  annotate("text",x=downgene$log2FoldChange[1:3],y=(-log10(downgene$padj[1:3])),label=downgene$Gene_symbol[1:3], size=5.0)
plot(p)



#换一种方式绘制火山图
library(ggpubr)
my_result$v <- -log10(my_result$padj)#生成新的一列v
#挑选要展示的基因
my_select <- subset(my_result, padj < 0.00000000000000000001 & abs(log2FoldChange) > 4)
my_select <- my_select$Gene_symbol
#绘图
ggscatter(my_result,
          x = "log2FoldChange",
          y = "v",
          ylab = "-log10(adjust p-value)",
          size = 2,
          color = "regulate",
          palette = c('#86D47D','#DFE0DF','#C34B99'),
          label = "Gene_symbol",
          label.select = c("GPR155"))
