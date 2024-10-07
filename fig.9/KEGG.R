setwd("D:/R/KEGG富集分析")

# 下载所需要的R包 ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 检查 clusterProfiler 包是否已经安装
if (!require("clusterProfiler", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
# 检查 org.Hs.eg.db 包是否已经安装
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
# 检查 ggplot2 包是否已经安装
if (!require("ggplot2", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggplot2")
}
# 检查 dplyr 包是否已经安装
if (!require("dplyr", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("dplyr")
}
# 检查 ggsci 包是否已经安装
if (!require("ggsci", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggsci")
}
# 加载所需要的R包 ----------------------------------------------------------------
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#绘图所需
library(dplyr)#整理数据

# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene <- read.csv("upGene_1_0.05.csv",row.names = 1)
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# KEGG富集分析 ----------------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        title = "KEGG_enrichment")#标题
dotplot(KEGG)


# 问题所在：KEGG数据库API更新了 Y叔的clusterProfiler包也需要更新-----------------------------------------------
#解决办法就是卸载掉富集分析所用的clusterProfiler包，再安装最新版本的
.libPaths() 
remove.packages("clusterProfiler",lib = file.path("你的R包安装路径"))
remove.packages("BiocManager",lib = file.path("你的R包安装路径"))
#重新安装相应的R包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("clusterProfiler")

#假如依然不行，可能是你的R版本太旧了
install.packages("installr")
library(installr)
updateR()
#更新过程中，记得选择把旧版本的R包转移到新版本

# 美化版本 --------------------------------------------------------------------
df <- as.data.frame(KEGG)

# 气泡图 ---------------------------------------------------------------------
# 排序并选择最显著的 10 个条目
df_sorted <- df[order(df$p.adjust), ][1:10, ]
colnames(df_sorted)
df_sorted$Description
df_sorted$Description <- gsub(" - Homo sapiens \\(human\\)", "", df_sorted$Description)
###
ggplot(df_sorted, aes(x=GeneRatio, y=Description, size=-log10(p.adjust), color=p.adjust)) + ## 生成散点图，以GeneRatio为x轴，Description为y轴，-log10(p.adjust)为点的大小，p.adjust为点的颜色
  geom_point(alpha=0.7) + # 绘制散点图层，设置透明度为0.7
  scale_color_gradient(low="blue", high="red") +
  labs(title="The most enrichment KEGG pathway", x="GeneRatio", y="KEGG pathway", size="-log10(P-value)", color="P-value") + #设置颜色映射为渐变色，低值为蓝色，高值为红色
  theme_minimal() + #设置图表的标题、x轴和y轴标签，以及点的大小和颜色的图例标签
  theme(axis.text.y = element_text(face = 'bold',size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold")) #设置主题为最小化主题，调整y轴的文本大小和字体粗细，以及标题、轴标签的字体大小
