rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_endo/')
#--------------------------内皮细胞亚群分析--------------------------
library(Seurat)
library(data.table)
library(tibble)
library(ggplot2)
## 读取内皮细胞注释结果
scRNA <- qs::qread("../results/Endo_anno.qs") # 读取单细胞数据
colnames(scRNA@meta.data)
#--------------------------S9a 内皮细胞亚群umap图--------------------------
dim(scRNA)
length(unique(scRNA$orig.ident))
length(unique(scRNA$celltype))
colnames(scRNA@meta.data)

## 首先设置细胞身份标识
Idents(scRNA) <- 'celltype' #!!!
table(scRNA$celltype)

## 定义颜色向量
colors <- c('#eebad2','#d5adcf','#c4d691','#d1d2e9',
            '#e8bf83','#b0c3e0','#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935','#343031',
            '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31'
)
colors <- c('#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935','#343031',
            '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31')
maincol<- colors[1:length(unique(scRNA$celltype))]

dput(as.character(unique(scRNA$celltype)))
# 将需要的颜色与celltype对应
names(maincol) <- c("VenECs_KLK3", "CapECs_RGCC", "ArtECs_ENPP2", "VenECs_CXCL8", 
                    "CapECs_ESM1", "VenECs_NOP53", "VenECs_SELE", "ECs_FABP4")
pdf(file = "./S9a_dimplot_Endo.pdf", width = 8, height = 7)
DimPlot(scRNA,label = TRUE,group.by = 'celltype', #!!!
        cols = maincol, # 使用预定义颜色向量
        pt.size = 1,
        raster = FALSE,
        reduction = "umap") +
  NoAxes()  +           # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "Endothelial cell") +  # 修改labs标签
  theme_dr(xlength = 0.2, # x轴箭头的长度比例
           ylength = 0.2, # y轴箭头的长度比例
           arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                               type = "closed")) +  # 箭头类型：实心闭合箭头
  # NoLegend() +  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        panel.grid = element_blank())
dev.off()

