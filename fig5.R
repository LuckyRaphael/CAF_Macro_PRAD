rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig5/')
# 加载必要的R包
library(data.table)   # 高效数据读取
library(GSVA)        # 基因集变异分析
library(tidyverse)   # 数据整理和可视化
library(survival)    # 生存分析
library(survminer)   # 生存分析可视化
library(ggsurvfit)   # 生存曲线绘制
library(gridExtra)   # 图形排版
library(dplyr)       # 数据处理
library(tidyr)       # 数据整理
library(stringr)     # 字符串处理
library(patchwork)   # 图形组合
####转录组生存分析####
exp <- read.csv('../results/TCGA_PRAD_TPM.csv',
                row.names = 1, check.names = FALSE)
table(substr(colnames(exp), 14, 15))
Group <- ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group <- factor(Group,levels = c("normal","tumor"))
table(Group)
# 筛选肿瘤样本
exp <- exp[, Group == "tumor"]

# 读取细胞亚群marker基因
genesets_mye <- as.data.frame(fread(
  input = "../results/Mye_allmarkers.csv",
  header = TRUE,data.table=FALSE)) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster)  # 选择gene和cluster列

genesets_fib <- as.data.frame(fread(
  input = "../results/CAF_allmarkers.csv",
  header = TRUE,data.table=FALSE)) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 100, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster)  # 选择gene和cluster列

genesets <- rbind(genesets_fib, genesets_mye)
colnames(genesets)
genesets <- split(genesets$gene, genesets$cluster)
# 重新排序, 将CAF排在前面
genesets <- genesets[order(!grepl("CAF", names(genesets), ignore.case = TRUE))]
names(genesets)
# 设置选择的细胞亚群名称
select_celltype<- c("POSTN+CAF", "Macro_APOE")
# 提取选择的细胞亚群对应的marker基因
genesets <- genesets[names(genesets) %in% select_celltype] 

gc()
# GSVA打分
gsva_res <- gsva(expr = as.matrix(exp),
                 gset.idx.list = genesets,
                 mx.diff = T, # 数据为正态分布则T，双峰则F
                 kcdf = "Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                 min.sz = 2, max.sz = 100,
                 verbose = T,
                 parallel.sz = 15,
                 method = "ssgsea") # 并行线程数目
gsva_res[1:nrow(gsva_res), 1:5]
gsva_res <- as.data.frame(t(gsva_res))
gsva_res[1:8,1:ncol(gsva_res)]
colnames(gsva_res)




###########################################################
theta1 <- read.csv("../Fig2/theta.csv", check.names = FALSE)
colnames(theta1)
colnames(theta1)[1] <- 'X'
theta1 <- tibble::column_to_rownames(theta1,'X')

theta2 <- read.csv("../results/theta_mye.csv", check.names = FALSE)
colnames(theta2)
colnames(theta2)[1] <- 'X'
theta2 <- tibble::column_to_rownames(theta2,'X')
identical(rownames(theta1),rownames(theta2))
gsva_res <- as.data.frame(cbind(theta1,theta2))
colnames(gsva_res)





### 读取Xena生存数据
# tcga_suv <- fread(
#   input = "../Fig2/XENA/TCGA_PRAD_surv.csv",
#   header = TRUE,       # 第一行是列名
#   data.table = FALSE
# )
# 
# # 数据处理
# tcga_suv <- tcga_suv %>%
#   tibble::column_to_rownames(var = "V1") # 将V1列转换为行名
# tcga_suv[1:6, 1:ncol(tcga_suv)]
# rownames(tcga_suv) <- paste0(rownames(tcga_suv), "-01A") # 标准化样本名格式

load('../Fig2/otherbulk/TCGA_PRAD.RData')
tcga_suv <-TCGA_PRAD[,c(1:2)]
rownames(tcga_suv) <- paste0(rownames(tcga_suv), "A") # 标准化样本名格式
# 数据清洗
table(tcga_suv$OS, useNA = "ifany")
tcga_suv$OS.time <- as.numeric(tcga_suv$OS.time)  # 转换生存时间为数值型
table(tcga_suv$OS.time, useNA = "ifany")

# 去除缺失值和生存时间为0的样本
tcga_suv <- tcga_suv %>% drop_na(OS.time)
tcga_suv <- tcga_suv[tcga_suv$OS.time > 0, ]
tcga_suv$OS.time <- as.numeric(tcga_suv$OS.time)


# tcga_suv <- tcga_suv[tcga_suv$OS.time > 30, ]
# 获取共同样本(生存数据和GSVA结果的交集)
ss <- intersect(rownames(tcga_suv), rownames(gsva_res))
tcga_suv <- tcga_suv[ss, ]
gsva_res <- gsva_res[ss, ]

# 检查样本一致性并合并数据
all(rownames(tcga_suv) == rownames(gsva_res))
rt <- cbind(tcga_suv, gsva_res)
colnames(rt)

# 数据预处理
table(rt$OS, useNA = "ifany")
range(rt$OS.time)
# rt$OS.time <- rt$OS.time/365  # 将生存时间从天转换为年
range(rt$OS.time)

# 获取需要分组的变量名(排除生存相关列)
dput(setdiff(colnames(rt), c("OS", "OS.time")))

#### 最佳截断值分组 ####
variables_to_analyze <- setdiff(colnames(rt), c("OS", "OS.time"))
# 使用survminer包寻找最佳截断值
best_thresholds <- surv_cutpoint(
  rt,
  time = "OS.time",
  event = "OS",
  variables = variables_to_analyze,
  minprop = 0.1,      # 最小比例限制
  progressbar = TRUE  # 显示进度条
)

# 根据最佳截断值进行分组
for (var in variables_to_analyze) {
  cutoff <- unname(best_thresholds[[var]][["estimate"]])
  rt[[var]] <- ifelse(rt[[var]] > cutoff, "High", "Low") 
}


# #### 中位值分组 ####
# variables_to_analyze <- setdiff(colnames(rt), c("OS", "OS.time"))
# # 根据中位值进行分组
# for (var in variables_to_analyze) {
#   median_value <- median(rt[[var]], na.rm = TRUE)  # 计算中位值
#   rt[[var]] <- ifelse(rt[[var]] > median_value, "High", "Low")
# }



# 查看分组结果
head(rt[, grep("+", colnames(rt)), drop = FALSE])
# 重命名指定列
names(rt)[names(rt) == "Macro_APOE"] <- "APOE+Macro"
# 设置选择的细胞亚群名称
select_celltype<- c("POSTN+CAF", "APOE+Macro")
select_celltype
# 自动提取marker名称
markers <- sapply(strsplit(select_celltype, "\\+"), function(x) x[1])
markers
# 动态创建组合分组变量
rt <- rt %>%
  mutate(Group = str_c(markers[1], "+", .[[select_celltype[1]]], "_",
                       markers[2], "+", .[[select_celltype[2]]]))
table(rt$Group)

# 设置因子水平
rt[[select_celltype[1]]] <- factor(rt[[select_celltype[1]]], levels = c('High','Low'))
rt[[select_celltype[2]]] <- factor(rt[[select_celltype[2]]], levels = c('High','Low'))
rt <- rt %>%
  mutate(Group = factor(Group, levels = c(
    paste0(markers[1], "+High_", markers[2], "+High"),  # High High
    paste0(markers[1], "+High_", markers[2], "+Low"),   # High Low  
    paste0(markers[1], "+Low_", markers[2], "+High"),   # Low High
    paste0(markers[1], "+Low_", markers[2], "+Low")     # Low Low
  )))

rt1 <- rt %>%
  filter(Group %in% c(
    paste0(markers[1], "+High_", markers[2], "+High"),  # High High
    paste0(markers[1], "+Low_", markers[2], "+Low")     # Low Low
  ))%>%
  mutate(Group = factor(Group, levels = c(
    paste0(markers[1], "+High_", markers[2], "+High"),  # High High
    paste0(markers[1], "+Low_", markers[2], "+Low")     # Low Low
  )))  # 重新设置因子水平

table(rt1$Group)
####Fig6A及S4D KMplot####
source('../resource/ggsurvplotKM.R')
# 两组的KM曲线
mycolor <- c("#E31A1C","#377EB8")
p1 <- ggsurvplotKM(rt[, c("OS.time", "OS", select_celltype[1])],
                   lables = c('High','Low'),
                   risk.table = FALSE,
                   title = paste0(markers[1], "+", strsplit(select_celltype[1], "\\+")[[1]][2]))
p1
dev.off()

p2 <-  ggsurvplotKM(rt[, c("OS.time", "OS", select_celltype[2])],
                    lables = c('High','Low'),
                    risk.table = FALSE,
                    title = paste0(markers[2], "+", strsplit(select_celltype[2], "\\+")[[1]][2]))
p2
dev.off()
# 四组的KM曲线
mycolor <- c("#c71000", "#8a4198")
p3 <- ggsurvplotKM(rt1[, c("OS.time", "OS", "Group")],
                   lables = levels(rt1$Group),
                   risk.table = FALSE,
                   title = paste0(select_celltype[1], "_", select_celltype[2]))
p3
dev.off()

# 使用 patchwork 包将所有图形排列成 3 列
p <- p1 + p2 + p3 + 
  plot_layout(ncol = 3)

# 显示最终图形
print(p)
dev.off()

# 保存图片
ggsave("./Fig6A_tcga_km.pdf", 
       p, width = 18, height = 6, dpi = 300)


