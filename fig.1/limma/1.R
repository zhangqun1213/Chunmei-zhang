setwd("D:/R/差异表达分析")

# 加载所需要的R包 ----------------------------------------------------------------
library(limma)
library(ggplot2) #用于绘制火山图
library(pheatmap) #用于绘制热图

# 输入表达矩阵和分组文件 -------------------------------------------------------------
expr_data<-read.table("Matrix_of_expression.txt",header = T,
                           row.names = 1,sep = "\t")
group<-read.csv("group.csv",header = T,row.names = 1,sep = ",")

# #构建分组矩阵--design ---------------------------------------------------------
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)

# #构建比较矩阵——contrast -------------------------------------------------------
contrast.matrix <- makeContrasts(Tumor-Normal,levels = design)

# #线性拟合模型构建 ---------------------------------------------------------------
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf)
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                   ifelse(DEG$logFC > 1, "up-regulated",
                          ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
table(DEG$regulate)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = "DEG_result.txt",
            sep = "\t",quote = F,row.names = F,col.names = T)

# 区分上下调基因 -----------------------------------------------------------------
DE_1_0.05 <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>1,]
upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated",]
downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated",]
write.csv(upGene_1_0.05,"upGene_1_0.05.csv")
write.csv(downGene_1_0.05,"down-regulated.csv")

# 火山图的绘制 ------------------------------------------------------------------
pdf("volcano.pdf")
ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ #x轴logFC,y轴adj.p.value
  geom_point(alpha=0.6,size=3.5,aes(color=regulate))+ #点的透明度，大小
  ylab("-log10(P.Value)")+ #y轴的说明
  scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ #logFC分界线
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ #adj.p.val分界线
  theme_bw()  #火山图绘制
dev.off()

# ggvolcano绘制另一种火山图 ----------------------------------------------------------
library(ggVolcano)
Genes <- rownames(DEG)
DEG$Genes <- Genes#将行名添加为一列
ggvolcano(data = DEG,x = "logFC",y = "P.Value",label = "Genes",
          label_number = 10,output = FALSE,
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value")

# 热图的绘制 -------------------------------------------------------------------
DEG_genes <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>1,]
DEG_gene_expr <- expr_data[rownames(DEG_genes),]
pdf("pheatmap.pdf")
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue","white","red"))(100), #颜色
         scale = "row", #归一化的方式
         border_color = NA, #线的颜色
         fontsize = 10, #文字大小
         show_rownames = F) 
dev.off()

