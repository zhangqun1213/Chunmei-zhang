setwd("D:/R/GEO注释")

# 下载所需要的R包 ----------------------------------------------------------------
# 检查 tidyverse 包是否已经安装
if (!require("tidyverse", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("tidyverse")
}

# 检查 GEOquery 包是否已经安装
if (!require("GEOquery", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  BiocManager::install("GEOquery")
}

# 检查 ggVennDiagram 包是否已经安装
if (!require("ggVennDiagram", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggVennDiagram")
}

# 加载所需要的R包 ----------------------------------------------------------------
library(GEOquery) # 数据下载
library(tidyverse) # 数据清洗
library(ggVennDiagram) # 绘制venn图，交集可视化

# 下载示例数据 ------------------------------------------------------------------
# ①表达矩阵
# 从GEO数据库下载GSE138043基因表达数据集 
eset = getGEO('GSE138043', destdir=".", AnnotGPL = F, getGPL = F)  

eset <- eset[[1]] # 选择第一个样本数据集
probes_expr <- exprs(eset) # 获得该数据集的基因表达矩阵 probes_expr  
dim(probes_expr)  # 查看基因表达矩阵的维度 
data_use <- as.data.frame(probes_expr) # 将基因表达矩阵转换为数据框data_use 

# ②注释信息
# 从GEO数据库下载平台GPL5175的注释信息ann_info
ann_info <- getGEO(GEO = 'GPL5175',destdir = ".")

anno_geoquery <- Table(ann_info) # 使用GEOquery包将ann_info转换为数据集anno_geoquery
colnames(anno_geoquery)# 查看列名
test <- anno_geoquery[c(1,11,15),c(1,10)]# 提取其中部分信息测试


# 拆分注释信息（测试） ------------------------------------------------------------------
# 拆分 gene_assignment 列（保留所有结果）
test_split_1 = test %>% 
  separate_rows(gene_assignment, sep = " /// ")  %>%  # 将`test`数据框中的`gene_assignment`列按照"///"分隔符进行拆分，生成新的行，并将拆分后的结果存储为新的数据框`split_df`
  separate(gene_assignment, into = c("Ensembl_ID", "gene_symbol", "gene_description", "location", "ENTREZ_ID"), sep = " // ") #将`gene_assignment`列按照"//"分隔符进行拆分，将拆分后的结果分别存储为`Ensembl_ID`、`gene_symbol`、`gene_description`、`location`、`ENTREZ_ID`五个列

# 拆分 gene_assignment 列（只保留第一个探针对应的Gene symbol）
test_split_2 = test %>% 
  separate_rows(gene_assignment, sep = " /// ")  %>%  #将`test`数据框中的`gene_assignment`列按照"///"分隔符进行拆分，生成新的行，并将拆分后的结果存储为新的数据框`split_df`
  separate(gene_assignment, into = c("Ensembl_ID", "gene_symbol", "gene_description", "location", "ENTREZ_ID"), sep = " // ")  %>%  #将`gene_assignment`列按照"//"分隔符进行拆分，将拆分后的结果分别存储为`Ensembl_ID`、`gene_symbol`、`gene_description`、`location`、`ENTREZ_ID`五个列
  distinct(ID, .keep_all = TRUE) 


# 拆分注释信息（正式） --------------------------------------------------------------
colnames(anno_geoquery)# 查看列名
anno_info <- anno_geoquery[,c("ID","gene_assignment")]# 只选取芯片列、注释信息列
anno_info <- anno_info[anno_info$gene_assignment != "---",]# 删除注释信息中的---

mydata <- anno_info %>%
  separate_rows(gene_assignment, sep = " /// ") %>% #将`anno_info`数据框中的`gene_assignment`列按照"///"分隔符进行拆分，生成新的行，并将拆分后的结果存储为新的数据框`split_df`
  separate(gene_assignment, into = c("Ensembl_ID", "gene_symbol", "gene_description", "location", "ENTREZ_ID"), sep = " // ") %>% #将`gene_assignment`列按照"//"分隔符进行拆分，将拆分后的结果分别存储为`Ensembl_ID`、`gene_symbol`、`gene_description`、`location`、`ENTREZ_ID`五个列
  distinct(ID, .keep_all = TRUE) 

# 向量的交集
intersect_array <- intersect(rownames(data_use), 
                             mydata$ID)
# 维恩图简单的可视化
list_ID <-list('表达矩阵芯片' = rownames(data_use),
               '注释信息芯片' = as.character(mydata$ID))
ggVennDiagram(list_ID, label_alpha=0) + scale_fill_gradient(low="gray", high = "skyblue")
# ggVennDiagram用于生成 Venn 径向图（radial Venn diagram）的函数，需要传入一个列表作为参数，
# label_alpha = 0：是控制标签透明度的参数，将其设为 0 表示隐藏标签。
# scale_fill_gradient调整 Venn 图填充颜色的函数，将低值的填充颜色设定为灰色，高值的填充颜色设定为天蓝色。