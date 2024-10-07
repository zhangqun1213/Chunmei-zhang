# 载入工作路径
setwd("D:/R/GO富集圈图")


# 下载所需要的R包 ----------------------------------------------------------------
# 检查 GOplot 包是否已经安装
if (!require("GOplot", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("GOplot")
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 检查 clusterProfiler 包是否已经安装
if (!require("clusterProfiler", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  BiocManager::install("clusterProfiler")
}
# 检查 org.Hs.eg.db 包是否已经安装
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  BiocManager::install("org.Hs.eg.db")
}

# 加载R包 --------------------------------------------------------------------
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(GOplot)#可视化

# 导入差异分析的结果 ---------------------------------------------------------------
DEG <- read.table("DEG_result.txt",header = T)

Threshold <- 0.0005 #设定阈值

DEG_use <- DEG[DEG$adj.P.Val < Threshold,] #过滤步骤

colnames(DEG_use) # 查看列名，下方选择列的依据
colnames(DEG_use)[1] <- "SYMBOL"# 修改列名，为下一步合并做准备
DEG_use <- DEG_use[,c("SYMBOL","logFC")] #只要基因名称和差异倍数两列




# 读取数据：基因差异分析的结果 ----------------------------------------------------------
Genes <- bitr(DEG_use$SYMBOL,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

data_use <- merge(DEG_use,Genes,by = "SYMBOL") #按照SYMBOL列合并数据，删除没有匹配到ENTREZID的基因行

# GO富集分析与整理 -------------------------------------------------------------
GO <- enrichGO(gene = data_use$ENTREZID, # 输入基因的 ENTREZID
               OrgDb = org.Hs.eg.db, # 注释信息来自 org.Hs.eg.db 数据库
               keyType = "ENTREZID", # 指定使用的基因标识符类型
               ont = "all", # 本次富集分析所使用的本体类型: "BP" (生物学过程)、"CC" (细胞组分)、"MF" (分子功能) 或 "all"（所有三个）
               pAdjustMethod = "BH", # 用于多重检验校正的 p 值校正方法 (Benjamini-Hochberg 方法)
               pvalueCutoff = 1, # p 值的阈值 (p 值低于此阈值的基因将被视为显著)
               qvalueCutoff = 1, # q 值的阈值 (q 值低于此阈值的基因将被视为显著)
               minGSSize = 5, # 基因集的最小大小 (基因的数量)
               maxGSSize = 5000, # 基因集的最大大小
               readable = TRUE) # 是否将 ENTREZID 转换为基因符号以提高可读性


# 将 GO 结果转换为数据框格式
GO_result <- as.data.frame(GO) 
# 选择显著富集功能通路，根据给定的阈值（Threshold）
GO_result <- GO_result[(GO_result$pvalue < 0.05 & GO_result$p.adjust < 0.05),] 
# 创建一个数据框，包含富集分析结果的GO类型
go_result <- data.frame(Category = GO_result$ONTOLOGY, 
                        ID = GO_result$ID, # 包含富集分析结果的通路 ID
                        Term = GO_result$Description, # 包含富集分析结果的通路描述
                        Genes = gsub("/", ", ", GO_result$geneID), # 包含富集分析结果中的基因 ID，多个 ID 用逗号分隔
                        adj_pval = GO_result$p.adjust) # 包含富集分析结果中校正后的 p 值
# 创建一个包含基因 ID 和差异表达值的数据框
genelist <- data.frame(ID = data_use$SYMBOL, logFC = data_use$logFC) 
# 将 genelist 数据框的行名设置为 ID
row.names(genelist)=genelist[,1] 
# 生成一个通路富集分析结果的可视化图表
circ <- circle_dat(go_result, genelist) 
head(circ)

# 可视化准备及绘图 ----------------------------------------------------------------
nrow(go_result) #查看一共有多少条通路
term_num = 3 # 设定图中要展示的GO的条目数量  

nrow(genelist) #查看一共有多少个基因
gene_num = nrow(genelist)  # 设定图中要展示的基因数量上限 

# 生成一个通路富集分析结果的和弦图表
chord <- chord_dat(circ, #一个包含通路富集分析结果的数据框
                   genelist[1:gene_num,], #一个包含基因 ID 和差异表达值的数据框，其中只包含前 gene_num 个基因
                   go_result$Term[1:term_num] # 一个包含通路描述的向量，其中只包含前 term_num 个通路
                   )

#开始绘图（生成PDF清晰度更高）

# 一、GO富集圈图
pdf(file="GO富集圈图.pdf",width = 12,height = 12)
GOChord(chord, 
        title = "GOcircos",
        # 标题
        space = 0.01,
        # 设置基因之间的间距为0.01
        gene.order = 'logFC',    
        # 根据logFC的值对基因进行排序
        gene.size = 2,
        # 基因标签的大小
        nlfc = 1,
        # 定义logFC列数（默认值=1）
        lfc.col=c('#F5053C','#EDF5F7','cyan1'),
        # 指定的logFC的填充颜色：c（低值的颜色，中点的颜色，高值的颜色）
        lfc.min = -3,
        # 指定logFC比例的最小值（默认值= -3）
        lfc.max = 3,
        # 指定logFC比例的最大值（默认值= 3）
        gene.space = 0.2,       
        # 基因离圆圈距离
        border.size = 0.2,  
        # 定义功能区边框的大小
        process.label = 12,
        # 图例的大小
        limit = c(1,5)
        # 控制着和弦图中GO条目与基因的数量,限制输出过多或过少的关系,保证可视化效果清晰易读。
        # 第一个值为1,表示每个基因至少被分配一个GO条目
        # 第二个值为5,表示每个展示的GO条目至少关联5个基因。
        )       
dev.off()

# 二、GO圆形树状图
pdf(file="GO圆形树状图.pdf",width = 10,height = 8)
GOCluster(circ, 
          metric = "euclidean",
          # 设置距离度量方法为欧氏距离
           clust.by='logFC',
          # 根据logFC值进行聚类
          lfc.col=c('#F5053C','#EDF5F7','cyan1'),
          # 指定的logFC的填充颜色：c（低值的颜色，中点的颜色，高值的颜色）
          as.character(go_result[1:term_num,3])
          # 选择go_result中的前termNum个GO术语
          )
dev.off()

# 三、太极图
GOCircle(circ,nsub = 8)

go_result$ID[1:20]
Select_GO_Term <- c('GO:0042113',
                    'GO:0042100',
                    'GO:0048731',
                    'GO:0030183',
                    'GO:0003341')
GOCircle(circ,nsub = Select_GO_Term)
# - 最外面显示GO条目,外圈表示基因（上调/下调）

# 四、气泡图
GOBubble(circ,
         labels = 9,#图中标签（GO：...）显示的阈值，对-log（p.adj）>6的term进行标注，数值越大，显示越少
         display = 'multiple',#在同一图上展示多个变量的气泡
         colour = c('#F5053C','#F9F871','cyan1'),#设置颜色
         bg.col = T # 是否添加背景颜色
         )
GOBubble(circ
         #,table.legend = F # 是否显示右侧GO具体条目
         )
#Z-score 的计算方式是通过基因在通路中的表达水平和通路中所有基因的表达水平的比较进行标准化。