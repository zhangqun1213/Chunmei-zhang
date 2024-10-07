# 加载工作路径 ------------------------------------------------------------------
setwd("D:/bilibiliR/38、IOBR免疫浸润part1")

# IOBR包下载所需代码 ----------------------------------------------------------------
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

depens <- c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(version = "3.17")(depen,update = FALSE)
}

if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")

# 加载IOBR包 -----------------------------------------------------------------
library(IOBR)

# 脚本中需要的其它R包 --------------------------------------------------------------
library(tidyverse)

# 载入示例数据 ------------------------------------------------------------------
load("TCGA-ESCA_mrna_expr_tpm.rdata")# 输入标准化的RNA-seq数据

#习惯性备份
mydata <- mrna_expr_tpm

# 查看肿瘤微环境_反卷曲_算法 ------------------------------------------------------------
tme_deconvolution_methods

# 返回特征估计的可用参数选项--------------------------------------------------------------
signature_score_calculation_methods
#在特征评分评估过程中采用了三种方法，包括单样本基因集富集分析（ssGSEA）、主成分分析（PCA）和Z-评分。

# 肿瘤微环境相关特征集 --------------------------------------------------------------
names(signature_tme)

# 代谢相关基因集 -----------------------------------------------------------------
names(signature_metabolism)

# 与生物医学基础研究相关的基因集：例如m6A和外泌体 ---------------------------------
names(signature_tumor)

# 所有免疫细胞的特征基因集 ------------------------------------------------------------
signature_collection

# 查看参考基因集文献出处 -------------------------------------------------------------
signature_collection_citation



# Part I 反卷曲_算法

# Cibersort免疫细胞浸润 ---------------------------------------------------------
cibersort <- deconvo_tme(eset = mydata, # 基因表达矩阵，可以是一个数值型矩阵或数据框，行名为基因，列名为样本名
                         method = "cibersort", # 指定'CIBERSORT'方法进行去卷积
                         arrays = F, # 输入的不是芯片数据
                         perm = 100 # 统计分析的置换次数（建议≥100次）影响'CIBERSORT'方法。
)

cell_bar_plot(input = cibersort[1:20,], # 输入的去卷积结果，选择前20个样本进行绘图
              title = "Cell Fraction",
              legend.position = "bottom", # 图例位置
              palette = 3, # 调色板颜色方案
              coord_filp = T, # 是否翻转坐标轴
              show_col = F # 是否显示细胞类型的颜色
              )
# 将堆积柱状图保存为PDF文件
ggsave("Cibersort.pdf", width = 10, height = 8)


# xcell免疫细胞浸润分析 -----------------------------------------------------------
xcell <- deconvo_tme(eset = mydata, # 基因表达矩阵，可以是一个数值型矩阵或数据框，行名为基因，列名为样本名
                     method = "xcell",# 指定'xcell'方法进行去卷积
                     arrays = FALSE # 输入的不是芯片数据
                     )
head(xcell)

# epic免疫细胞浸润分析 ------------------------------------------------------------
epic <- deconvo_tme(eset = mydata, # 基因表达矩阵，可以是一个数值型矩阵或数据框，行名为基因，列名为样本名
                    method = "epic", # 指定'epic'方法进行去卷积
                    arrays = FALSE # 输入的不是芯片数据
                    )
head(epic)

cell_bar_plot(input = epic[1:20,], # 输入的去卷积结果，选择前20个样本进行绘图
              title = "Cell Fraction",
              legend.position = "bottom", # 图例位置
              palette = 3, # 调色板颜色方案
              coord_filp = T, # 是否翻转坐标轴
              show_col = F # 是否显示细胞类型的颜色
)
# 将堆积柱状图保存为PDF文件
ggsave("epic.pdf", width = 10, height = 8)

# 还有很多其它方法，大家可以尽情尝试
# MCPcounter               EPIC              xCell 
# "mcpcounter"             "epic"            "xcell" 
# CIBERSORT        CIBERSORT Absolute          IPS 
# "cibersort"         "cibersort_abs"         "ips" 
# ESTIMATE                 SVR               lsei 
# "estimate"              "svr"             "lsei" 
# TIMER                quanTIseq 
# "timer"             "quantiseq" 


# 分割线 ---------------------------------------------------------------------
# 函数用法介绍（此处不需要运行！！！） -------------------------------------------------------------
deconvo_tme(eset, # 基因表达矩阵，可以是一个数值型矩阵或数据框，行名基因名，列名为样本名
            project = NULL, # 项目名称，用于区分不同的数据集，默认为空
            method = tme_deconvolution_methods, # 字符串，指定去卷积方法，支持多种方法
            arrays = FALSE, # 适用于微阵列数据的模式，默认为FALSE。目前影响'CIBERSORT'， 'svr'和'xCell'方法。
            tumor = TRUE, # 适用于肿瘤样本的签名矩阵/过程，默认为TRUE。目前影响'EPIC'方法。
            perm = 1000, # 统计分析的置换次数（建议≥100次）。目前影响'CIBERSORT'和'svr_ref'方法。
            reference, # 免疫细胞基因矩阵，例如lm22、lm6或使用generateRef/generateRef_rnaseq生成的矩阵
            scale_reference, # 逻辑值，指示是否对参考文件进行缩放。如果为TRUE，则参考文件中的值将在行方向上进行居中和缩放。目前影响'svr'和'lsei'方法。
            plot = FALSE, # 是否绘制图形。目前影响'IPS'方法。
            scale_mrna, # 逻辑值。如果为FALSE，则禁用对不同细胞类型的mRNA含量进行校正。这由计算绝对分数的方法支持（EPIC和quanTIseq）。
            group_list = NULL, # 肿瘤类型列表的样本。
            platform = "affymetrix", # 字符串，指示平台类型。默认为"affymetrix"。目前影响'ESTIMATE'方法。
            absolute.mode = FALSE, # 在绝对模式下运行CIBERSORT或svr（默认值为FALSE）。
            abs.method = "sig.score", # 如果设置了absolute为TRUE，则选择方法：'no.sumto1'或'sig.score'。
)