setwd("D:/bilibiliR/16、xcell免疫浸润")

# 加载所需要的R包 ----------------------------------------------------------------
library(devtools)
#install_github('dviraran/xCell')
library(xCell)
library(tidyverse)#数据处理
library(readxl)#数据读取
library(ggstatsplot)#绘图
library(ggpubr)#绘图
library(pheatmap)#绘制热图
library(survival)#生存分析
library(survminer)#生存曲线

# 读取表达矩阵文件 ----------------------------------------------------------------
load("TPM.Rdata")
#表达矩阵以样本为核心去重
expr_use <- mRNA_TPM_new
expr_use = expr_use[,sort(colnames(expr_use))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
table(substr(colnames(expr_use),14,16))#查看肿瘤样本和正常样本数量
colnames(expr_use) <- substr(colnames(expr_use),1,12)
expr_data <- as.matrix(expr_use)#转为矩阵

#xCell考虑了10800个基因，使用如下代码可以获取这个基因列表。
genes <- xCell.data$genes
#xCell要求你的数据集至少包含这个10800个基因中的5000个，如果你的数据集包含的这些基因太少了，将会导致结果不准确
#看特征基因集
xCell.data$signatures
#xCell.data$signatures 表示一个基因集集合，其中包含了 489 个基因集，
#每个基因集都有一个名称和一些唯一标识符。这个集合还包含了一些其他信息，
#如基因标识符的类型和基因集集合的类型等
geneset <- read_excel("489细胞型基因特征.xlsx")

# 开始分析 --------------------------------------------------------------------
#rawEnrichmentAnalysis#原始富集分析
#transformScores#转换分数:将原始富集分数转换为类似于百分比的线性比例，xCell使用预先计算的校准参数进行转换。
#spillOver#溢出:应用于由xCell产生的溢出补偿，以减少相关性很强的的细胞类型之间的依赖性。 
scores <-  rawEnrichmentAnalysis(expr_data, signatures = xCell.data$signatures,
                                 genes = xCell.data$genes,parallel.sz = 4,
                                 parallel.type = "SOCK")
tscores <-  transformScores(scores,fit.vals = xCell.data$spill$fv,scale = T)
result <- spillOver(tscores,K = xCell.data$spill$K, alpha = 0.5)

#一步搞定
result <- xCellAnalysis(expr_data)


#数据简单处理
result1 <- as.data.frame(t(result))
result1 <- result1[,c(65:67,1:64)]
summary(result1[,1:3])
result1 <- result1[,1:3]

# 我们能用这个数据做什么呢 ------------------------------------------------------------
View(expr_use)#文件一：表达矩阵文件
View(result1)#文件二：样本的免疫评分+基质评分+微环境得分
surv_data <- read.csv("生存信息.csv",row.names = 1)#文件三：患者的生存信息


# 根据免疫评分/生存状态分组 -----------------------------------------------------------
#数据简单处理
expr_use <- t(expr_use)
combined_data1 <- cbind(result1,expr_use)
rownames(surv_data) <- surv_data$submitter_id
surv_data <- surv_data[,-1]
combined_data2 <- merge(surv_data, combined_data1, by=0, all=FALSE)

data_use <- combined_data2[,1:36]
rownames(data_use) <- data_use$Row.names
data_use <- data_use[,-1]
table(data_use$status)

# 绘制箱式图 -------------------------------------------------------------------
colnames(data_use)
ggplot(data_use, aes(x=status, y=ImmuneScore,color=status)) + 
  geom_boxplot() + geom_point(position=position_jitter(width=0.1)) +#添加数据点
  theme(text = element_text(size=12),
        axis.title = element_text(size=18)) + 
  xlab("status") + #为X轴和Y轴添加标签
  ylab("ImmuneScore") +
  scale_x_discrete(labels=c("Alive", "Dead")) + #更改X轴的标签
  theme_minimal(base_size=12, base_family="sans") + scale_color_manual(values = c("#60B2F2","#E12939")) + stat_compare_means(method = "t.test")#调整图形的主题+更改颜色图例的颜色+在图上添加组间比较

ggboxplot(data_use, "status", "StromaScore", color = "status",palette = c("#E12939", "#2F74E5")) + stat_compare_means(method = "t.test")
#ggstatsplot绘图
ggbetweenstats(data = data_use,x = status,y = MicroenvironmentScore,plot.type = "box")


# 绘制热图 --------------------------------------------------------------------
Immune_median_score <- median(data_use$ImmuneScore)
Immune_median_score
data_use$Immune_group <- cut(data_use$ImmuneScore, breaks = c(-Inf, Immune_median_score, Inf), labels = c("low", "high"))
data_use_1 <- data_use[,c(36,6:35)]
#-----------------------
row_index <- order(data_use_1$Immune_group)#按免疫评分分组排序
data_use_1 <- data_use_1[row_index, ]#重新排序
table(data_use_1$Immune_group)#查看高低分组的数目
heat_data <- t(data_use_1[,2:31])
heat_data <- as.data.frame(heat_data)
scaled_heat_data <- scale(heat_data)
#分组文件
group <- data_use_1[,1,drop=FALSE]

pheatmap(scaled_heat_data,
         color = colorRampPalette(c("#0084E4", "white", "#E12939"))(50),
         annotation_col = group,
         cluster_row = T,#分别指定是否按列和行聚类。
         cluster_col = FALSE,
         border="black",#边框颜色
         main = "Heatmap",#指定图表的标题
         show_rownames = TRUE,#是否展示行名
         show_colnames = F,#是否展示列名
         cexCol = 1,#指定列标签的缩放比例。
         scale = 'row',#指定是否应按行方向或列方向居中和缩放，或不居中和缩放。对应的值为row, column和none。
         angle_col = "90",#指定列标签的角度。
         legend = TRUE,#指定是否显示图例。
         legend_breaks=c(-5,0,5),#指定图例中显示的数据范围为-10到10。
         fontsize_row = 10,#分别指定行标签和列标签的字体大小。
         fontsize_col = 10)

# 绘制生存曲线 ------------------------------------------------------------------
data_use_surv <- data_use[,c(36,1,2)]
summary(data_use_surv$time)#检查时间的单位
time_year <- (data_use_surv$time)/365#换算成以年为单位
data_use_surv$time_year <- time_year#添加一列时间（年）
data_use_surv$status <- ifelse(data_use_surv$status == "Dead",1,0) #死亡为1，存活为0
table(data_use_surv$status)#查看生存状态
colnames(data_use_surv)#查看行名
fit<-survfit(Surv(time_year,status)~Immune_group,data = data_use_surv)#拟合生存函数

ggsurvplot(fit, 
           data = data_use_surv,#你的数据集
           surv.median.line = "hv", #添加中位生存期
           size = 1, # 改变线条大小
           cex.lab=2,
           break.time.by = 1, #横坐标时间间隔
           xlim = c(0,6),
           axis.title.x =element_text(size=5), 
           axis.title.y = element_text(size=5),
           palette = c("#F04972","#7BB0E0"),#自定义调色板
           conf.int = TRUE, #添加置信区间
           pval = TRUE, #添加p值
           risk.table = TRUE, #添加风险表
           xlab = "Follow-up years",  
           ylab="Survival probability ",
           risk.table.col = "strata",#按组划分的风险表颜色 
           legend.labs =  c("Immune score low","Immune score high"), #更改图例标签
           risk.table.height = 0.3,  
           ggtheme = theme_bw() 
)
