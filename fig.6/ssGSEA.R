setwd("D:/bilibiliR/11、ssGSEA免疫浸润")
rm(list=ls(all=TRUE))

# 加载所需要的R包 ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")

library(GSVA)#算法
library(tidyverse)#数据处理
library(ggpubr)#绘图
library(ggplot2)#绘图
library(pheatmap)#绘制热图

#ssGSEA是一种算法，根据我们输入的参考基因集进行相应的分析。

# 准备数据:转录组+免疫细胞类型列表+分组数据 -------------------------------------
DEG_expr <- read.csv("DEG_expr.csv",row.names = 1)
group <- read.csv("group.csv",row.names = 1)
markergenes <- read.csv("markergenes.csv")
table(markergenes$Cell.type)
#数据检验
boxplot(DEG_expr,outline=F, notch=F , las=2)
qx <- as.numeric(quantile(DEG_expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)#判断是否进行log的标准
LogC
if (LogC) { 
  DEG_expr[which(DEG_expr <= 0)] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}  

#数据处理
#盛夏的果实函数科普：
#lapply 函数是 R 语言的一个内置函数，用于对列表或向量中的每一个元素执行函数。
geneset <- split(markergenes,markergenes$Cell.type)

im_geneset <- lapply(geneset, function(x){
  gene = x$Metagene
  unique(gene)
})
lapply(im_geneset[1:3], head)

save(im_geneset,file = "im_geneset.Rdata")

DEG_expr <- as.matrix(DEG_expr) 
# 开始进行ssGSEA分析 ------------------------------------------------------------
result <- gsva(DEG_expr,im_geneset,method = "ssgsea")
result1 <- as.data.frame(t(result))

write.csv(result1,"ssGSEA_result.csv")

# 结果可视化 -------------------------------------------------------------------
# 热图 ----------------------------------------------------------------------
pdf(file="heatmap.pdf",width = 10,height = 8)
pheatmap(result1,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border="black",#边框颜色
         main = "Heatmap",#指定图表的标题
         show_rownames = T,#是否展示行名
         show_colnames = T,#是否展示列名
         cexCol = 1,#指定列标签的缩放比例。
         scale = 'row',#指定是否应按行方向或列方向居中和缩放，或不居中和缩放。对应的值为row, column和none。
         cluster_col=T,#分别指定是否按列和行聚类。
         cluster_row=F,
         angle_col = "45",#指定列标签的角度。
         legend = F,#指定是否显示图例。
         legend_breaks=c(-3,0,3),#指定图例中显示的数据范围为-3到3。
         fontsize_row = 10,#分别指定行标签和列标签的字体大小。
         fontsize_col = 10)
dev.off()  

# 分组箱式图 -------------------------------------------------------------------
#数据处理
data <- cbind(result1,group)
colnames(data)
data <- data[,c(28,29,1:27)]
data <- pivot_longer(data = data,
                      cols = 3:29,
                      names_to = "celltype",
                      values_to = "proportion")
#开始绘图
pdf(file="分组箱式图+显著性P值.pdf",width = 10,height = 8)
ggboxplot(data = data,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",#绘制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图。
          merge = FALSE,#是否将相同值的分组合并。
          color = "black",#箱形图边框的颜色。
          fill = "group",#箱形图填充色。
          palette = c("#1C3EDF","#DF1C26"),#颜色调色板。
          title = NULL,#图形标题。
          xlab = "ssGSEA",#x 轴标签。
          ylab = "Expression",#y 轴标签
          bxp.errorbar = FALSE,#是否在箱形图中绘制误差条。
          bxp.errorbar.width = 0.2,#误差条宽度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的标签
          short.panel.labs = TRUE,#是否将分面标签缩短
          linetype = "solid",#线条类型
          size = NULL,#图形大小。
          width = 0.8,#箱形图的宽度。
          notch = FALSE,#是否在箱形图中绘制刻度。
          outlier.shape = 20,#异常值标记的形状。
          select = NULL,#要绘制的变量
          remove = NULL,#不要绘制的变量。
          order = NULL,#箱形图的排序方式。
          error.plot = "pointrange",#如何绘制误差，可以是 "pointrange" 或 "errorbar"。
          label = NULL,#要添加的标签
          font.label = list(size = 12, color = "black"),#标签的字体属性
          label.select = NULL,#要添加标签的数据点
          repel = TRUE,#是否使用 repel 库的功能使标签互不重叠
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
          stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = F,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")))
dev.off()  
