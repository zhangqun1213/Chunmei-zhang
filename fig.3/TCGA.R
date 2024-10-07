setwd("D:/R/TCGA临床数据")

# 下载所需要的R包 ----------------------------------------------------------------
# 检查compareGroup是否安装
if (!requireNamespace("compareGroups", quietly = TRUE)) {
  install.packages("compareGroups")
}

# 加载所需包 -------------------------------------------------------------------
library(compareGroups)

# 载入示例数据 ------------------------------------------------------------------
load("TCGA-ESCA_clinical.rdata")

identical(rownames(clin_info),colnames(mrna_expr_counts))
# 提取你要的数据所在的列 -------------------------------------------------------------
colnames(clin_info)
# 根据临床信息列名，来判断你自己要提取的数据
# xx = clin_info$xx (向量名称 = 你的临床数据框$数据框中的列)

# ①提取病人的ID
ID = clin_info$barcode
# ②提取病人的生存状态
status = clin_info$vital_status
# ③提取病人的生存时间
time1 = clin_info$days_to_last_follow_up
time2 = clin_info$days_to_death
# 两个时间，生存分析以days_to_death死亡为准，如果为na，则选择days_to_last_follow_up
# ④提取病人的性别
gender = clin_info$gender
# ⑤提取病人的年龄
age = clin_info$age_at_index
# ⑥提取病人的分期信息
stage = clin_info$ajcc_pathologic_stage

stage_T = clin_info$ajcc_pathologic_t
stage_N = clin_info$ajcc_pathologic_n
stage_M = clin_info$ajcc_pathologic_m

# 合并生存信息 ------------------------------------------------------------------
mydata <- data.frame(ID,
                     status,
                     time1,
                     time2,
                     gender,
                     age,
                     stage,
                     stage_T,
                     stage_N,
                     stage_M)
data <- mydata #习惯性备份

# 开始处理筛选获得的临床数据 -----------------------------------------------------------
# ①先把行名换成样本名称
rownames(data) <- data$ID
data <- subset(data, select = -ID)#删除ID这一列
# ②生存结局变量处理（转换成因子变量/数值变量）
table(data$status)
data$status <- as.factor(data$status)
data$status <- ifelse(data$status == "Alive", 0, 1) # 转成数值型
# ③生存时间处理（合并两个生存时间）
data$time <- ifelse(is.na(data$time2), data$time1, data$time2)
data <- subset(data, select = -c(time1, time2))#删除没用的时间列
# ④性别（二分类变量）的处理
table(data$gender)
data$gender <- as.factor(data$gender)
data$gender <- ifelse(data$gender == "female", 0, 1) # 转成数值型
# ⑤年龄（连续型变量）的处理
summary(data$age)
data$age <- as.numeric(data$age)
# 假如有缺失，最常见的就是均值/中位数填补
data[1, 5] <- NA #[行，列]
data$age[is.na(data$age)] <- 61
# ⑥分期（有序多分类/无序多分类变量）的处理
table(data$stage)
# 首先判断是否真正缺失（有些有T N M分期，但是stage缺失，可以自行单个填补
data[8, 4] <- "我在这里" #[行，列]
# 遇到缺失值，最简单的办法，就是删除缺失的样本
newdata <- data[!is.na(data$stage), ]# 生成一个新的数据框，删除缺失的样本
# 假如你明确知道想要填补的某一个值，比如M分期中，未知的都可以填Mx
data$stage_M[is.na(data$stage_M)] <- "Mx"
data$stage[is.na(data$stage)] <- "stage x"

# N分类转为少一点的分类
table(data$stage)
data$stage_copy <- ifelse(data$stage %in% c("Stage I", "Stage IA", "Stage IB"), "Stage I",
                              ifelse(data$stage %in% c("Stage II", "Stage IIA", "Stage IIB"), "Stage II",
                                     ifelse(data$stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), "Stage III",
                                            ifelse(data$stage %in% c("Stage IV", "Stage IVA"), "Stage IV",
                                                   ifelse(data$stage == "stage x", "stage x", NA)))))
table(data$stage_copy)

# 数据保存与三线表制作 ----------------------------------------------------
colnames(data)
data <- data[,c("status","time","gender","age","stage_copy",
                     "stage_T","stage_N","stage_M")]
write.csv(data,"临床数据.csv")
# 简易三线表
tab_1 <- descrTable(status ~ ., # 指定描述status分组下，所有的变量
                   data = data)
print(tab_1)
# 表格输出，保存为word版本
export2word(tab_1, file='临床数据表格.docx') 