rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig6/')
####空间转录组分析####
####绘制Fig4EFG&S3E####
# devtools::install_local("../SpaGene-master.zip")
library(Seurat)
library(ggplot2)
library(scRNAtoolVis)
library(RColorBrewer)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(SingleCellExperiment)
library(BayesSpace)
library(GSVA)
# 创建主文件夹
main_fig_dir <- "./stplot/"
if (!dir.exists(main_fig_dir)) {
  dir.create(main_fig_dir, recursive = TRUE)
}

## 处理函数
## Seurat deal example one ST
process_spatial_data <- function(st, 
                                 assay = "Spatial", 
                                 npcs = 20,
                                 resolution = 0.5,
                                 ndims = 30) {
  # 1. 数据标准化
  st <- SCTransform(st, assay = assay, verbose = TRUE)
  # 2. 降维分析
  st <- RunPCA(st, assay = "SCT", verbose = FALSE)
  # 3. 可视化PCA结果
  print(ElbowPlot(st, ndims = ndims, reduction = "pca"))
  dev.off()
  # 4. 聚类分析
  st <- FindNeighbors(st, reduction = "pca", dims = 1:npcs)
  st <- FindClusters(st, resolution = resolution, verbose = FALSE)
  # 打印聚类结果
  print(table(st@meta.data$seurat_clusters))
  # 5. UMAP可视化
  st <- RunUMAP(st, reduction = "pca", dims = 1:npcs)
  # 返回处理后的对象
  return(st)
}

#查看空转样本所在文件夹
dput(list.files(path = "../data/sp/ST_all/",
                full.names = FALSE, recursive = FALSE))

# 获取所有样本
samples <- list.files(path = "../data/sp/ST_all/",
                      full.names = FALSE, recursive = FALSE)

# 调色板
color <- unique(c('#2867a9','#ca3135','#f57d34','#8469ae','#469d47',
                  '#936962','#9EDDCE','#E3CCE2','#F9C1AC','#FCD6A8','#FBF9D4',
                  '#B8D7EA','#F8BBB5','#BCE5DF','#DEDBED',"#93CA87",
                  "#ED614A","#F4D020",brewer.pal(n = 9, name = "Set1"), 
                  brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 12, name = "Set3")))

####读取基因集及细胞亚群marker基因打分####
## 读取HALLMARK基因集
gmt <- read.gmt('../resource/h.all.v2024.1.Hs.symbols.gmt')
## 去掉HALLMARK_前缀并创建基因集列表
gmt <- split(gmt$gene, gsub("^HALLMARK_", "", gmt$term))
# 设置选择的基因集名称
select_gmt <- c("HYPOXIA")
# 提取选择的基因集对应的marker基因
gmt <- gmt[names(gmt) %in% select_gmt] 

################################################################################
# 读取细胞亚群marker基因集
genesets_main <- as.data.frame(fread(
  input = "../results/Maintype_allmarkers.csv",
  header = TRUE,data.table=FALSE)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 100, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster)  # 选择gene和cluster列

genesets_mye <- as.data.frame(fread(
  input = "../results/Mye_allmarkers.csv",
  header = TRUE,data.table=FALSE)) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 100, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster) # 选择gene和cluster列

genesets_fib <- as.data.frame(fread(
  input = "../results/CAF_allmarkers.csv",
  header = TRUE,data.table=FALSE)) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 100, wt = avg_log2FC) %>%
  dplyr::select(gene, cluster) # 选择gene和cluster列

geneset <- rbind(genesets_main, genesets_fib, genesets_mye)
colnames(geneset)
geneset <- split(geneset$gene, geneset$cluster)
names(geneset)
# 设置选择的细胞亚群名称
select_celltype<- c("Epithelia", "T", "POSTN+CAF", "Macro_APOE")
# 提取选择的细胞亚群对应的marker基因
geneset <- geneset[names(geneset) %in% select_celltype] 
names(geneset)

# 将gmt和geneset合并
genesets <- c(gmt, geneset)
names(genesets)
names(genesets) <- gsub("/", "", names(genesets))
names(genesets)

