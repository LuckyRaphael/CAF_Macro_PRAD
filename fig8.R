rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig8/')
####Fig4####
rm(list = ls())
gc()
library(GSVA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(stringr)
####绘制相关性散点图####
# 基因集准备
# 读取巨噬细胞(macro)特征基因集
genesets_macro=as.data.frame(fread("../results/Mye_allmarkers.csv", 
                                   header = TRUE, data.table = FALSE)) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster)

# 读取成纤维细胞(fib)特征基因集
genesets_fib=as.data.frame(fread("../results/CAF_allmarkers.csv", 
                                 header = TRUE, data.table = FALSE)) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 100, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster)

# 合并两个基因集
genesets <- rbind(genesets_macro, genesets_fib)

# 查看基因集的列名
colnames(genesets)

# 将基因集转换为GSVA需要的格式(以cluster分组的基因列表)
genesets <- split(genesets$gene, genesets$cluster)

# 替换基因集名称中的特殊字符为"_"
names(genesets) <- gsub("[-+]", "_", names(genesets))

# 查看基因集名称
names(genesets)

# 设置选择分析的细胞亚群名称
select_celltype<- c("POSTN_CAF", "Macro_APOE") #!!!

# 提取选择的细胞亚群对应的marker基因
genesets<- genesets[names(genesets) %in% select_celltype] 



####循环绘制bulk数据相关性散点图####
# 设置数据目录
data_dir <- "../Fig2/otherbulk/"
rdata_files <- list.files(path = data_dir, 
                          pattern = "\\.RData$", 
                          full.names = TRUE)
rdata_files

# 初始化图片列表
plot_list <- list()

# 循环处理每个RData文件
for (file in rdata_files) {
  # 加载数据
  var_name <- sub("\\.RData$", "", basename(file))
  cat("Processing:", file, "\n")
  load(file)
  exp <- get(var_name)
  
  exp <- na.omit(exp)
  # GSVA分析
  GSVA_res <- gsva(expr = as.matrix(t(exp[, -c(1:2)])),
                   gset.idx.list = genesets,
                   mx.diff = TRUE,
                   kcdf = "Gaussian",
                   parallel.sz = 15,
                   method = "ssgsea")
  
  GSVA_res <- as.data.frame(t(GSVA_res))
  
  # 绘图
  current_plot <- ggplot(GSVA_res, aes_string(x = select_celltype[1], y = select_celltype[2])) +
    geom_point(col = "black") +
    geom_smooth(method = lm, se = TRUE, na.rm = TRUE, fullrange = TRUE, size = 2, col = "darkred") +
    stat_cor(method = "spearman", digits = 3, size = 5) +
    ggtitle(var_name) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype = "solid"),
          axis.line = element_blank(),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(color = "black", size = 14))
  
  plot_list[[var_name]] <- current_plot
}

length(plot_list)
# 拼图
combined_plot <- wrap_plots(plot_list, ncol = 3)
# 显示
print(combined_plot)
dev.off()
# 保存图片
ggsave('./multi_cohort_corplot.pdf', combined_plot, height = 9, width =9)





