setwd("D:/R/GEO数据集合并")

# 这段代码要解决的问题 --------------------------------------------------------------
#问题一：多个GEO数据集合并（消除批次效应）
#问题二：从官网下载的数据，读入+ID转换
#问题三：网络问题，某些数据集下载注释信息失败（AnnoProbe对基因组数据进行注释)

# 加载所需要的R包 ----------------------------------------------------------------
library(devtools)#安装github上的R包
library(GEOmirror)#GEO镜像
library(GEOquery)#GEO数据下载
#---------------
#install_github("jmzeng1314/AnnoProbe")
library(AnnoProbe)#对基因组数据进行注释
library(limma)#标准化
library(sva)#去批次

#sva包是一种用于探索和消除在基因表达数据中批次效应的方法。
#批次效应是指那些不属于研究中的主要兴趣，但是在数据中仍然存在的变量。
#这些变量可能是由于实验过程中的差异而导致的，如样本来源、采集时间等。

#定义用得到的色板
mycolors <- c('#54D9ED','#324B4F','#95B0B5','#D1BCFE','#FCFCD4','#F78ABA',
              '#00A9CB','#8EF8B6','#ED546C','#1BBEB9','#DE90DE','#00A44F')
              
# 下载数据集一 ------------------------------------------------------------------
#读取自己从官网下载的注释文件
expr_data <- read.table("GSE31821_series_matrix.txt", sep = "\t",header = T, fill = T,row.names = 1,comment.char = '!')#comment.char = '!'表示文件中以'!'开头的行都被视为注释
Anno_Info <- read.table("GPL570-55999.txt",sep = "\t",quote = "",header = T,fill = T)
#简单整理
colnames(Anno_Info)
Anno_Info <- Anno_Info[,c("ID","Gene.Symbol","ENTREZ_GENE_ID")]#提取自己需要的列

eset_1 <- expr_data
boxplot(eset_1,outline=F, notch=F , las=2,col='#54D9ED')#绘制简易箱式图
dim(eset_1)#查看样本数和基因数
#假如数据不均一
#eset_1 <- normalizeBetweenArrays(eset_1)

#数据清洗与ID转换
same1 <- match(rownames(eset_1),Anno_Info$ID)#找出相同的元素的位置
eset_1$Gene_Symbol <- Anno_Info[same1,c("Gene.Symbol")]#添加Gene_Symbol进入表达矩阵
eset_1 <- eset_1[!duplicated(eset_1$Gene_Symbol),]# 会保留 "eset_1" 中 "Gene_Symbol" 列中第一次出现的行，并将其他重复出现的行删除。

#在前面GEO数据ID转换时，我们采取重复gene_symbol取平均值的方法；
#TCGA数据ID转换时，我们取基因表达量最大值的方法；
#在这里，我们采取直接删除重复基因的方法。
eset_1 <- eset_1[!is.na(eset_1$Gene_Symbol)&!grepl('//',eset_1$Gene_Symbol),]
eset_1 <- eset_1[eset_1$Gene_Symbol!="",]#删除空值
rownames(eset_1) <- eset_1$Gene_Symbol#修改行名
colnames(eset_1)
eset_1 <- eset_1[,-7]#删除不需要的列（这里换成你自己表达矩阵的最后一列）

save(eset_1,file = "GSE31821.Rdata")

# 下载数据集二 ------------------------------------------------------------------
eset_2_raw <- geoChina('GSE66724')
Anno_Info_2 <- idmap('GPL570',type = 'soft')
#获取表达矩阵
eset_2_raw
eset_2_raw <- eset_2_raw[[1]]
eset_2 <- exprs(eset_2_raw)#提取表达矩阵
eset_2 <- as.data.frame(eset_2)
dim(eset_2)

boxplot(eset_2,outline=F, notch=F , las=2,col='#E93639')#绘制简易箱式图
#假如数据不均一
#eset_2 <- normalizeBetweenArrays(eset_2)

#数据清洗与ID转换
same2 <- match(rownames(eset_2),Anno_Info_2$ID)#找出相同的元素的位置
eset_2$Gene_Symbol <- Anno_Info_2[same2,c("symbol")]#添加Gene_Symbol进入表达矩阵
eset_2 <- eset_2[!duplicated(eset_2$Gene_Symbol),]# 会保留 "eset_2" 中 "Gene_Symbol" 列中第一次出现的行，并将其他重复出现的行删除。

eset_2 <- eset_2[!is.na(eset_2$Gene_Symbol)&!grepl('//',eset_2$Gene_Symbol),]
eset_2 <- eset_2[eset_2$Gene_Symbol!="",]#删除空值
rownames(eset_2) <- eset_2$Gene_Symbol#修改行名
colnames(eset_2)
eset_2 <- eset_2[,-17]#删除不需要的列（这里换成你自己表达矩阵的最后一列）


# 数据集合并 -------------------------------------------------------------------
identical(eset_1,eset_2)#查看两个数据集是否一致，答案是否定的
eset_2 <- eset_2+2 #演示需要，人为调整
same_genes <- intersect(rownames(eset_1),rownames(eset_2))#查看两个表达矩阵共同的基因
eset_merge <- cbind(eset_1[same_genes,],eset_2[same_genes,])#提取共同基因所在的行并进行合并
col <- c(rep("#54D9ED", 6), rep("#E93639", 16))#前六列样本显示蓝色，后面 16 列样本显示红色
par(mar = c(8, 3, 1, 2) + 0.1)#调整图中四周留白区域，用于完整显示X轴标签
boxplot(eset_merge, col=col,las=2)#整体箱式图绘制


# 读取分组数据 ------------------------------------------------------------------
group_1 <- read.csv("group.csv",header = T)#手动制作的分组数据

phenoDat <- pData(eset_2_raw)#提取GSE66724的临床信息（分组数据）
colnames(phenoDat)#查看有什么临床信息
group_2 <- phenoDat[,c("geo_accession","treatment_protocol_ch1")]#提取有用的列（此处修改为你要提取的列名）
colnames(group_2) <- c("Sample","group")#修改列名
rownames(group_2) <- NULL#去除行名
group_2$group <- "atrial fibrillation"
#合并分组数据
group <- rbind(group_1,group_2)


# 正式开始去批次 -----------------------------------------------------------------
GSE <- c(rep('GSE31821',6),rep('GSE66724',16))
GSE
group_list <- group$group
group_list
table(group_list,GSE)

data <- eset_merge
batch <- c(rep('GSE31821',6),rep('GSE66724',16))
design <- model.matrix(~group_list)
#用sva包的removeBatchEffect去除批次效应
expr_limma <- removeBatchEffect(data,batch = batch,design = design)
boxplot(expr_limma, col=col,las=2)#整体箱式图绘制
#对比看去除批次效应前后效果
par(mfrow=c(1, 2))
boxplot(eset_merge, col=col,las=2,main = "before")
boxplot(expr_limma, col=col,las=2,main = "after")
