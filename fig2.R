rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig2/')
library(Seurat)
library(scRNAtoolVis)
library(patchwork)
#### 读取注释好的成纤维细胞亚群#####
scRNA <- qs::qread("../results/Fib_anno.qs")
table(scRNA$celltype)

dim(scRNA)
colnames(scRNA@meta.data)
scRNA$cell_subtypes <- scRNA$celltype
table(scRNA$seurat_clusters)
# UMAP降维可视化
# n.neighbor，越大越连续，越小越分散  mindist越大越cluster越松散
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:10,
                 n.neighbors = 40, min.dist = 0.5,
                 reduction.name = "umap"
)
table(scRNA$seurat_clusters)
# qs::qsave(scRNA,"../results/Fib_anno.qs")

table(scRNA$cell_subtypes)
# 自定义颜色
colors<- unique(c('#0073c2','#b40f1f',
                  '#306311','#df772b','#e184ca','#de554c','#e8e684',
                  '#b2dfe8','#667570','#5ec9b8','#314e9b','#616e88',
                  '#7623f1','#eaeaea',
                  '#c5cdf7','#7092d2','#e69fc4','#b6854d','#8dd1c6','#fffeb3','#bcb9d8','#80b2d4','#feb662',
                  '#b4df65','#d9d9d9','#be7fbf','#cbecc4','#ffee72','#e6bdcb'))

# 绘制seurat_clusters的umap图
pdf(file = "./S2A_dimplot_fibro.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'seurat_clusters', #!!!
        cols = colors, # 使用预定义颜色向量
        pt.size = 2,
        raster = TRUE,
        reduction = "umap") +
  NoAxes()  +           # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "") +  # 修改labs标签
  theme_dr(xlength = 0.2, # x轴箭头的长度比例
           ylength = 0.2, # y轴箭头的长度比例
           arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                               type = "closed")) +  # 箭头类型：实心闭合箭头
  # NoLegend() +  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        # 设置图例字体大小
        legend.text = element_text(size = 12),  # 图例项文字大小
        legend.title = element_text(size = 14), # 图例标题文字大小
        panel.grid = element_blank())
dev.off()



tissuecol<- c( '#9e0242','#fcae62','#d6404e','#6050a5','#a72b24','#aadda3','#e5f499','#3188bd','#fee18b','#67c2a5',
               '#00b4a9','#049964','#0968b3','#04c0ec','#94288e','#f47a1f','#f59f83','#eb1275')
tissuecol<- tissuecol[1:length(unique(scRNA$tissue))] #自定义大群颜色
dput(as.character(unique(scRNA$tissue)))
#将需要的颜色与cell type对应
names(tissuecol) <- c("Primary", "mLN", "CRPC", "Normal", "mCRPC", "LN")
# 绘制tissue的umap图,根据组织类型split
pdf(file = "./S2A_dimplot_fibro_tissue.pdf", width = 15, height = 4)
DimPlot(scRNA,label = FALSE,group.by = 'tissue', #!!!
        split.by = 'tissue', #split
        cols = tissuecol, 
        pt.size = 4,
        raster = TRUE,
        reduction = "umap") +
  NoAxes()  +           # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "") +  # 修改labs标签
  theme_dr(xlength = 0.2, # x轴箭头的长度比例
           ylength = 0.2, # y轴箭头的长度比例
           arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                               type = "closed")) +  # 箭头类型：实心闭合箭头
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        panel.grid = element_blank())
dev.off()



####细胞周期评分####
library(Seurat)
library(ggplot2)
library(reshape2)
# 获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(scRNA))

# 获取S期相关基因
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(scRNA))

# 细胞周期阶段评分
scRNA <- CellCycleScoring(scRNA, g2m.features=g2m_genes, s.features=s_genes)

colnames(scRNA@meta.data)
table(scRNA$Phase)
table(scRNA$cell_subtypes, scRNA$Phase)
df <- as.data.frame(table(scRNA$cell_subtypes, scRNA$Phase))
colnames(df) <- c("Subtype", "Phase", "Count")
# 排序子群
df$Subtype <- factor(df$Subtype, levels = unique(df$Subtype))

# 绘制细胞比例雷达图
pdf(file = "./S2C_cellcycle.pdf", width = 8, height = 6)
ggplot(df, aes(x = Subtype, y = Count, fill = Phase)) +
  geom_bar(stat = "identity") +
  coord_polar(start = 0) +
  theme_minimal() +
  scale_fill_manual(values = c("G1" = '#fffeb3', "G2M" = '#bcb9d8', "S" = '#80b2d4')) +
  theme(
        axis.text = element_text(color = "black"),  # 所有坐标轴文字为黑色，大小12
        axis.text.x = element_text(size = 12, margin = margin(t = 10, b = 10)),  # 增加上下边距
        axis.text.y = element_text(size = 8),  # x轴文字大小
        axis.title = element_text(color = "black", size = 12),  # 坐标轴标题
        legend.text = element_text(color = "black", size = 12),  # 图例文字
        legend.title = element_text(color = "black", size = 14),  # 图例标题
        legend.position = "right") +
  ylab("") + xlab("")
dev.off()



colors<- unique(c('#eebad2','#d5adcf','#c4d691','#d1d2e9',
                  '#e8bf83','#b0c3e0','#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935','#343031',
                  '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31',
                  '#b2dfe8','#667570','#5ec9b8','#314e9b','#616e88',
                  '#c5cdf7','#7092d2','#e69fc4','#b6854d','#8dd1c6','#fffeb3','#bcb9d8','#80b2d4','#feb662',
                  '#b4df65','#d9d9d9','#be7fbf','#cbecc4','#ffee72','#e6bdcb'))
####绘制Fig2A Dimplot####
CAFcol<- colors[1:length(unique(scRNA$cell_subtypes))]
dput(as.character(unique(scRNA$cell_subtypes)))
#将需要的颜色与cell type对应
names(CAFcol) <- c("NDUFA4L2+CAF", "CCL2+CAF", "TAGLN+CAF", "CXCL8+CAF", "RUNX2+CAF", 
                   "APOD+CAF", "POSTN+CAF", "MMP2+CAF", "MCAM+CAF")

pdf(file = "./Fig2A_dimplotplot_fibrocelltype.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'cell_subtypes', #!!!
        cols = CAFcol, # 使用预定义颜色向量
        pt.size = 3,
        raster = TRUE,
        reduction = "umap") +
  NoAxes()  +           # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "") +  # 修改labs标签
  theme_dr(xlength = 0.2, # x轴箭头的长度比例
           ylength = 0.2, # y轴箭头的长度比例
           arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                               type = "closed")) +  # 箭头类型：实心闭合箭头
  # NoLegend() +  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        # 设置图例字体大小
        legend.text = element_text(size = 12),  # 图例项文字大小
        legend.title = element_text(size = 14), # 图例标题文字大小
        panel.grid = element_blank())

dev.off()





