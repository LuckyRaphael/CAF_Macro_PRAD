rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig7/')
####泛癌单细胞数据读取####
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(ggplot2)
library(Seurat)
library(tidydr)
scRNA <- qs::qread("/home/w282308/JIC_CRC/data/pan/scRNA_pancancer_anno.qs")
dim(scRNA)
table(scRNA$Sample)
length(unique(scRNA$Sample))
table(scRNA$tissue)
table(scRNA$Datasets)
table(scRNA$Cancer_type)
sample_stats <- scRNA@meta.data %>%
  distinct(Sample, .keep_all = TRUE) %>%
  group_by(Datasets) %>%
  summarise(
    Total_Samples = n(),
    Normal_Samples = sum(tissue == "Normal"),
    Tumor_Samples = sum(tissue == "Tumor")
  )

print(sample_stats)

table(scRNA$Cancer_type,scRNA$Sample)
table(scRNA$celltype)
dput(unique(scRNA$celltype))
cols <- c(
  '#005f73','#099396','#ca6702','#ed9b00','#a12c2f','#bb3f02',
  '#b15d2e','#e82036','#237bba','#40ac53',
  '#f1c742','#7498b5','#e7e7e7','#d66969',
  '#e5f499','#3188bd','#fee18b','#67c1a6',
  '#9e0242','#d6404e','#fcae62','#6050a5','#aadda3','#e5f499','#3188bd','#fee18b','#67c2a5',
  "#DC050C", "#FB8072", "#1965B0", "#196580", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#F7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#Ba5ce3", "#808000",
  "#aeae5C", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
length(cols)
Maincol<- cols[1:length(unique(scRNA$celltype))] #自定义大群颜色
dput(as.character(unique(scRNA$celltype)))
#将需要的颜色与cell type对应
names(Maincol) <- c("T/NK cells", "Endothelial cells", "Myeloid", "Prolif", "B cells", 
                    "Plasma cells", "Mast cells", "Epithelial cells", "Fibroblasts"
)

pdf(file = "./pancancer_dimplotplot.pdf", width = 7, height = 5)
#umap图
DimPlot(scRNA, reduction = "umap",label=T,
        group.by = 'celltype',cols = Maincol,
        raster=TRUE,pt.size = 1)& 
  NoAxes()  &          # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "") &  # 修改labs标签
  tidydr::theme_dr(xlength = 0.2, # x轴箭头的长度比例
                   ylength = 0.2, # y轴箭头的长度比例
                   arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                                       type = "closed")) &  # 箭头类型：实心闭合箭头
  # NoLegend() &  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        panel.grid = element_blank())

dev.off()
 


pdf(file = "./pancancer_dimplotplot_Datasets.pdf", width = 7, height = 5)
#umap图
DimPlot(scRNA, reduction = "umap",label=F,
        group.by = 'Datasets',cols = cols[6:22],
        raster=TRUE,pt.size = 1)& 
  NoAxes()  &          # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "") &  # 修改labs标签
  tidydr::theme_dr(xlength = 0.2, # x轴箭头的长度比例
                   ylength = 0.2, # y轴箭头的长度比例
                   arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                                       type = "closed")) &  # 箭头类型：实心闭合箭头
  # NoLegend() &  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        panel.grid = element_blank())

dev.off()



