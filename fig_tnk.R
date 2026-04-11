rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_tnk/')
#--------------------------T细胞亚群分析--------------------------
library(Seurat)
library(data.table)
library(tibble)
library(ggplot2)
## 读取T细胞注释结果
scRNA <- qs::qread("../results/T_NK_cell.qs") # 读取单细胞数据
#--------------------------S8a T细胞亚群umap图--------------------------
dim(scRNA)
length(unique(scRNA$orig.ident))
length(unique(scRNA$celltype))

## 首先设置细胞身份标识
Idents(scRNA) <- 'celltype' #!!!
table(scRNA$celltype)

### 定义颜色向量
colors <- c('#cf5d51','#dd8349','#ed9baf','#2d6089',
            '#4d85a0','#a1d0d3','#c1d7e8','#6dadc3',
            '#698a5c','#b2d081','#e1ba32','#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935','#343031',
            '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31',
            "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#CAB2D6","#8DD3C7", "#FB8072"
            
)
maincol<- colors[1:length(unique(scRNA$celltype))]

dput(as.character(unique(scRNA$celltype)))
# 将需要的颜色与celltype对应
names(maincol) <- c("CD8_C3_ISG15", "CD4_C3_FOXP3", "CD4_C4_CD40LG", "CD4_C1_CCR7", 
                    "CD8_C1_GZMK", "CD8_C4_NDUFAF8", "Prolif_STMN1", "NK_C1_FGFBP2", 
                    "CD8_C2_CCL4", "NK_C2_TYROBP", "CD4_C6_CLU", "CD4_C7_CXCR4", 
                    "CD4_C5_CXCL13", "CD4_C2_TCF7")
pdf(file = "./S8a_dimplot_T.pdf", width = 8, height = 7)
DimPlot(scRNA,group.by = 'celltype', #!!!
        cols = maincol, # 使用预定义颜色向量
        label = FALSE,
        pt.size = 1.5,
        raster = TRUE,
        reduction = "umap") +
  NoAxes()  +           # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "T/NK cell") +  # 修改labs标签
  theme_dr(xlength = 0.2, # x轴箭头的长度比例
           ylength = 0.2, # y轴箭头的长度比例
           arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                               type = "closed")) +  # 箭头类型：实心闭合箭头
  # NoLegend() +  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        panel.grid = element_blank())
dev.off()



#--------------------------S8b T细胞亚群特定基因Featureplot--------------------------
select_markers <- c('CD3E','CD4','CD8A','GNLY')
# mycol <- c("lightgray","#f32a1f")
# # 蓝色系配色
# mycol <- c("#bebebe", "#0f5688")
# 红色系配色
mycol <- c("#bebebe", "#7c2124")
# mycol <- c("#bebebe", "#953653")
pdf(file = "./S8b_featureplot_T.pdf", width = 8, height = 7)
FeaturePlot(scRNA, 
            order = TRUE,
            features = select_markers, 
            col = mycol,
            # min.cutoff = 'q30',
            raster = TRUE,
            ncol = 2,
            pt.size=2, # 如果raster的话需要点调大
            reduction = "umap") & NoLegend() &
  NoAxes()          # 隐藏默认坐标轴标签和刻度线
dev.off()

