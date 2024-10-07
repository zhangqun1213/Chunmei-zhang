
# 加载当前工作路径 ----------------------------------------------------------------
setwd("D:/bilibiliR/45、免疫浸润相关性棒棒糖图")

# 下载所需要的R包 ----------------------------------------------------------------
# 检查 ggplot2包是否已经安装
if (!require("ggplot2", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggplot2")
}

# 加载所需包 -------------------------------------------------------------------------
library(ggplot2)

# 加载示例数据 -------------------------------------------------------------------
data <- read.csv("correlation_result.csv", row.names = 1,header = T, sep = ",", stringsAsFactors = F)

# 设定目标基因
target_gene <- "SHISA3"
mydata <- data[data$Gene == target_gene,]# 选择目标基因的数据
head(mydata)

# 数据处理 -----------------------------------------------------------------------
# 使用reorder函数重新排列y轴
mydata$im_cell <- reorder(mydata$im_cell, mydata$Cor)

# 开始绘图 --------------------------------------------------------------------
plot <- ggplot(mydata, aes(x=Cor, y=im_cell)) +   # 绘制散点图
        geom_segment(aes(x = 0, xend = Cor, y = im_cell, yend = im_cell), color="black") + # 添加线段图层
        geom_point(aes(size=abs(Cor), colour=p.value), alpha=0.5) + # 添加点图层 
        scale_colour_gradient(low = "#339D5A", high = "#EDB306") + # 设置颜色渐变 
        #ED063F #339D5A #EDB306 #备选颜色参数
        # 设置点的大小
        scale_size_continuous(range = c(2, 10)) + # 设置点的大小范围
        theme_minimal() +  # 设置主题 
        # theme(legend.position = "none") +  # 去掉图例
        # 调整坐标轴
        # scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) + # 设置x轴
        # scale_y_discrete(limits = rev(levels(mydata$im_cell))) + # 设置y轴
        # 调整文本样式
        theme(axis.line = element_line(size = 1.0),
              axis.text = element_text(size = 12, face = "bold"),
              axis.text.y = element_text(hjust = 1),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 将标题居中
              axis.title = element_text(size = 14, face = "bold")) +
       # 添加标签  
        labs(title = "目标基因SHISA3与免疫细胞相关性图", xlab="Cor", ylab="im_cell")

plot # 显示图形

# 覆盖不同的主题
plot + theme_grey() + ggtitle("theme_grey()")
plot + theme_bw() + ggtitle("theme_bw()")
plot + theme_linedraw() + ggtitle("theme_linedraw()")
plot + theme_light() + ggtitle("theme_light()")
plot + theme_dark() + ggtitle("theme_dark()")
plot + theme_classic() + ggtitle("theme_classic()")

# 保存图片
ggsave("target_cor.png", width = 8, height = 6, dpi = 300)


# Tips:
# 要快速掌握用ggplot绘图可以试试2023-08-16发布的ggplotAssist教程绘图 -----------------------------
# 从github上下载ggplotAssist
devtools::install_github("cardiomoon/ggplotAssist") 
# 安装完后关闭当前Rstudio，重新打开
library(ggplotAssist)