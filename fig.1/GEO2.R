setwd("D:/R/GEO数据下载")

# 加载所需要的R包 ----------------------------------------------------------------
#install.packages("devtools")
library(devtools)
#install.packages("GEOquery")
#install_github("jmzeng1314/GEOmirror")
#install_github("jmzeng1314/idmap1")

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

# 另外的下载方式 -----------------------------------------------------------------
install_github("jmzeng1314/AnnoProbe")
library(AnnoProbe)

# 下载GEO数据和注释信息 ---------------------------------------------------------
eset <- geoChina('GSE53819')#方法一，调用国内镜像，需要GEOmirror包
eset = getGEO('GSE53819', destdir=".", AnnotGPL = F, getGPL = F)#方法二，更推荐这种，需要GEOquery包
ann_info <- getGEO(GEO = 'GPL6480',destdir = ".")

# 获取表达矩阵 ------------------------------------------------------------------
eset
eset <- eset[[1]]
probes_expr <- exprs(eset)
dim(probes_expr)
data_use <- as.data.frame(probes_expr)

# 获取临床信息 ------------------------------------------------------------------
phenoDat <- pData(eset)
dim(phenoDat)
Samples <- row.names(phenoDat)
phenoDat <- cbind(Samples,phenoDat)
colnames(phenoDat)
group <- phenoDat[,c("Samples","source_name_ch1")] 
rownames(group) <- NULL
group$source_name_ch1 <- c(rep("Tumor",18),rep("Normal",18))
colnames(group) <- c("Samples","group")
group
write.csv(phenoDat,"clinicaldata.csv",)
# 判断表达矩阵是否需要对数处理 --------------------------------------------------
boxplot(data_use,las=2)#las完整展示行名
qx <- as.numeric(quantile(data_use, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) #判断是否进行log的标准
LogC
if (LogC) { 
  data_use[which(data_use <= 0)] <- NaN
  data_use <- log2(data_use) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}   #进行log2取对数

# 数据检验：箱线图，PCA图、层次聚类图 --------------------------------------------
#①箱线图绘制
boxplot(data_use,las=2)#las完整展示行名
# 均数不一致需标准化的情况 ------------------------------------------------------------
data_use=normalizeBetweenArrays(data_use)
boxplot(data_use,outline=FALSE, notch=T , las=2)
#②PCA绘制
PCA_result <- PCA(t(data_use), graph = F)
fviz_pca_ind(PCA_result,
             geom.ind = c("point","text"), #展示点和样本名
             mean.point = F, #不展示中心点
             repel = T, #将样本名根据位置进行展示，保证不出绘图边界
             col.ind = group$group, # 按分组给不同的颜色
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = T, # 绘制置信圈
             legend.title = "Groups")
#③层次聚类图
Sample_clust <- dist(t(data_use)) #计算变量间距离
hc <- hclust(Sample_clust) #hclust进行聚类
plot(hc,
     hang = -1,
     cex = 0.8)

# 获取注释信息 ------------------------------------------------------------------
Meta(ann_info)$title
anno_geoquery <- Table(ann_info)

# 探针ID转换为gene_symbol -------------------------------------------------------
colnames(anno_geoquery)
id_symbol <- anno_geoquery %>%
  select(ID,GENE_SYMBOL) %>%
  separate(GENE_SYMBOL,c("GENE_SYMBOL", NA),sep ="///") #对基因名进行操作
id_symbol[,2] <- trimws(id_symbol[,2]) #去除gene_symbol列的空格
#####
data_use$probe_id <- rownames(data_use)#添加一列用于添加基因名
data_with_name <- merge(data_use,
                        id_symbol,
                        by.x = "probe_id",
                        by.y = "ID")  #merge函数添加基因名
dim(data_with_name)
data_with_name <- data_with_name[data_with_name$GENE_SYMBOL != "",] #删除gene_symbol为空的数据
dim(data_with_name)
table(duplicated(data_with_name[,ncol(data_with_name)]))
data_with_name <- avereps(data_with_name[,-c(1,ncol(data_with_name))],
                          ID = data_with_name$GENE_SYMBOL) #对重复gene_symbol取平均值
table(duplicated(rownames(data_with_name))) #再判断gene_symbol是否有重复值
data_with_name['GAPDH',] #检测内参的表达量
data_with_name['ACTB',]
write.table(data.frame(gene_symbol=rownames(data_with_name),data_with_name),
            file = "Matrix_of_expression.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
write.csv(group,"group.csv")
