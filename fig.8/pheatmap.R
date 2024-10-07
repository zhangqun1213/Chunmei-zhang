# 加载工作路径 ------------------------------------------------------------------
setwd("D:/R/热图绘制")

# 下载用到的R包 -----------------------------------------------------------------
# 检查 psych包是否已经安装
if (!require("psych", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("psych")
}

# 检查 pheatmap 包是否已经安装
if (!require("pheatmap", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("pheatmap")
}

# 加载需要的R包 -----------------------------------------------------------------
library(psych)
library(pheatmap)

# 载入示例数据 ------------------------------------------------------------------
load("D:/bilibiliR/41、显著性热图绘制/TCGA-ESCA_mrna_expr_tpm.rdata")
cli <- read.csv("临床数据.csv",header = T,row.names = 1,sep = ",")

# 提取部分数据用于绘图
mydata <- mrna_expr_tpm[1:30,] # 此步骤意在提取感兴趣的基因用于绘图，可以自行挑选感兴趣的基因
cli_use <- cli[,1:5] # 提取要进行相关性分析的临床性状

# 对表达矩阵数据进行转置
mydata <- as.data.frame(t(mydata)) 
# 将临床数据中的变量转为数值型
cli_use$status <- ifelse(cli_use$status == "Dead",1,0) #死亡为1，存活为0
cli_use$time <- as.numeric(cli_use$time)# 转为数值型
cli_use$gender <- ifelse(cli_use$gender == "male",1,0) #男性为1，女性为0
cli_use$age <- as.numeric(cli_use$age)# 转为数值型
cli_use$stage_copy <- ifelse(cli_use$stage_copy == "Stage I", 1, 
                        ifelse(cli_use$stage_copy == "Stage II", 2,
                               ifelse(cli_use$stage_copy == "Stage III", 3, 
                                      ifelse(cli_use$stage_copy == "Stage IV", 4,
                                             ifelse(cli_use$stage_copy == "stage x", 4, NA)))))
# 查看整理好的数据
mydata[1:5,1:5] # 行名样本名，列名基因名
cli_use[1:5,1:5] # 行名样本名，列名临床信息


# 计算相关性 -------------------------------------------------------------------
result <- corr.test(x = mydata,# 输入要计算相关性的表达矩阵
                    y = cli_use,# 输入临床信息文件（需要是数值型数据）
                    method = "pearson" #method=“pearson”是默认值。备选项是"spearman"和“kendall”
                    )
# 相关性R值
result_R <- result$r
# 相关性P值
result_P <- round(result$p,3) #对P值保留三位有效数字

# 可视化之前的准备(对数据进行转置)
data <- t(result_R)
P <- t(result_P)

# P值的替换（设定一定的阈值进行替换）
P [P >= 0 & P < 0.001] <- "***"
P [P >= 0.001 & P < 0.01] <- "**"
P [P >= 0.01 & P < 0.05] <- "*"
P [P >=0.05 & P <= 1] <- ""


# 热图可视化 -------------------------------------------------------------------
pheatmap(data,
         display_numbers = P,
         color = colorRampPalette(c("#128BAE", "white", "#E8870D"))(100),
         border="white",#边框颜色
         main = "Heatmap",#指定图表的标题
         show_rownames = TRUE,#是否展示行名
         show_colnames = TRUE,#是否展示列名
         cexCol = 1,#指定列标签的缩放比例。
         scale = 'none',#指定是否应按行方向或列方向居中和缩放，或不居中和缩放。对应的值为row, column和none。
         angle_col = "45",#指定列标签的角度。
         legend = TRUE,#指定是否显示图例。
         legend_breaks=c(-0.2,0,0.2),#指定图例中显示的数据范围为-0.2到0.2。
         fontsize_row = 10,#分别指定行标签和列标签的字体大小。
         fontsize_col = 10)
