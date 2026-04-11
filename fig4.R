rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig4/')
####髓系细胞亚群注释结束####
library(Seurat)
library(scRNAtoolVis)
library(patchwork)
library(tidydr)
# 读取注释好的髓系细胞亚群
scRNA <- qs::qread("../results/Mye_anno.qs")
colnames(scRNA@meta.data)
scRNA$cell_subtypes <- scRNA$celltype


####绘制Fig5B 髓系亚群Dimplot####
cols <- unique(c(
  '#93bdd1','#ea7e38','#edaf74','#2190b4','#add26b',
  '#c8e3c0','#efec96','#bc917d','#d6a6a5','#d8585c',
  '#f5d2df','#e88e90',
  '#8dd1c6','#fffeb3','#bcb9d8','#80b2d4','#feb662',
  '#b4df65','#d9d9d9','#be7fbf','#cbecc4','#ffee72','#e6bdcb',
  
  '#306311','#df772b','#e184ca','#de554c','#e8e684',
  '#b2dfe8','#667570','#5ec9b8','#314e9b','#616e88',
  '#7623f1','#eaeaea',
  
  "#DC050C", "#FB8072", "#1965B0", "#196580", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#F7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#Ba5ce3", "#808000",
  "#aeae5C", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"))
Myecol<- cols[1:length(unique(scRNA$cell_subtypes))]
dput(as.character(unique(scRNA$cell_subtypes)))
#将需要的颜色与cell type对应
names(Myecol) <- c("Macro_APOE", "Macro_CD83", "cDC2_FCER1A", "Prolif_MKI67", 
                   "Macro_CXCL10", "Mono_FCGR3A", "Mono_S100A12", "Neutro_FCGR3B", 
                   "cDC2_CD1C", "Neutro_MPO", "Macro_CCL3")


pdf(file = "./Fig5B_dimplot.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'cell_subtypes', #!!!
        cols = Myecol, # 使用预定义颜色向量
        pt.size = 2.5,
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









####绘制Marker基因Dotplot####
library(RColorBrewer)
# 提取数据自行绘图
colnames(scRNA@meta.data)
Idents(scRNA) <- 'cell_subtypes'

# 已有的marker基因集
featureSets <- list(
  Macro_APOE = c( 'APOC1', 'GPNMB', 'LGMN', 'SLC40A1','APOE'),
  Macro_CD83 = c('BAG3', 'IER5', 'RGS1', 'IL1B', 'CD83'),
  Mono_S100A12 = c( 'S100A8', 'S100A9', 'VCAN', 'EREG','S100A12'),
  cDC2_FCER1A = c('BIRC3', 'PPA1', 'TXN', 'CPVL', 'FCER1A'),
  Neutro_FCGR3B = c( 'CXCR2', 'CSF3R', 'FPR1', 'TREM1','FCGR3B'),
  Mono_FCGR3A = c('SMIM25', 'LILRB2', 'CDKN1C', 'LST1', 'FCGR3A'),
  Prolif_MKI67 = c('TUBA1B', 'STMN1', 'TOP2A', 'CENPF','MKI67'),
  Neutro_MPO = c('ELANE', 'PRTN3', 'AZU1', 'CTSG','MPO'),
  Macro_CXCL10 = c('ISG15', 'GBP1', 'MX1', 'GBP4','CXCL10'),
  cDC2_CD1C = c( 'ARL4C', 'EZR', 'JAML', 'RUNX3','CD1C'),
  Macro_CCL3 = c( 'GABARAP', 'TYMP', 'UCP2', 'PRELID1','CCL3')
)
# 合并所有基因为一个向量
all_genes <- unlist(featureSets)
# 使用Seurat自带气泡图函数绘图然后提取绘图数据
p=DotPlot(scRNA,
          features = featureSets) + RotatedAxis()+
  theme(axis.text.x = element_text(size = 8))
p
dev.off()


library(magrittr)
result_data <- p$data %>% # 提取气泡图数据
  dplyr::select(feature.groups, features.plot, id, avg.exp.scaled, pct.exp) %>%
  dplyr::rename(Function = feature.groups, 
                Gene = features.plot, 
                Group = id, 
                AvgExpr = avg.exp.scaled, 
                PctExpr = pct.exp) %>%
  dplyr::mutate(
    Function = factor(Function, levels = names(featureSets)), #  确保功能和基因有正确的顺序
    Gene = factor(Gene, levels = rev(all_genes)),
    Group = factor(Group, levels = names(featureSets)) # 调整细胞类型的顺序
  )


pdf(file = "./Mye_markergene_dotplot.pdf", width = 12, height = 5)
# 分面
ggplot(result_data, aes(y = Group, x = Gene, color = AvgExpr, size = PctExpr)) +
  geom_point() +
  scale_color_gradientn(
    # colours = rev(brewer.pal(11, "Spectral")),
    # colours = rev(brewer.pal(11, "RdBu")),      # 红蓝发散（常用）
    # colours = rev(brewer.pal(11, "RdYlBu")),    # 红黄蓝发散
    # colours = rev(brewer.pal(11, "RdYlGn")),    # 红黄绿发散
    # colours = rev(brewer.pal(11, "PRGn")),      # 紫绿发散
    colours = rev(brewer.pal(11, "PiYG")),      # 粉绿发散
    # colours = rev(brewer.pal(11, "BrBG")),     # 棕绿发散
    guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"),
    name = "Relative\nExpression"
  ) +
  scale_size_continuous(name = "Percentage\nExpressed", range = c(0, 6)) +
  # facet_grid(. ~ Function, scales = "free_x", space = "free_x") +  # 横向分面
  theme_minimal() +
  xlab("") + 
  ylab("")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12,color = "black", vjust = 0.5),  # 调整x轴文本
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),  # 去除坐标轴线
    ## 第二种背景
    panel.border = element_rect(fill=NA,color="#35A79D", size=1, linetype="solid"), # 添加边框
    panel.background = element_rect(fill = "#F1F6FC"),# "white"设置主题(背景)
    # strip.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = c( "lightblue")), 
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white")   # 设置绘图区域背景为白色
  )
dev.off()



####绘制S7A Featureplot####
JIC_markers <- c(
  "CD163","MRC1", "CCR5", # Macrophage
  "CD14","FCN1", # Monocyte
  "IDO1","CLEC9A", # cDC1
  "CD1C","CLEC10A", # cDC2
  "LILRA4","JCHAIN", # pDC
  "LAMP3", # mDC
  "CSF3R","S100A8","S100A9", # Neutrophils
  "TOP2A","MKI67","STMN1" # Proliferation
)

pdf(file = "./S7A_featureplot.pdf", width = 15, height = 9)
FeaturePlot(scRNA, raster=T, order=T, 
            features=JIC_markers, 
            col = c("#e0f2f7", "#b3d8e6", "#80bed9", "#4aa3c8", "#1f8abf", "#1570a6", "#0f5688"),
            # col=c("lightgray","#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),# c("lightgray","#953553")
            pt.size=3, reduction="umap",
            ncol=6)& 
  theme_void() &  # 移除背景和网格
  theme(
    legend.frame = element_rect(colour = "black"),        
    legend.ticks = element_line(colour = "black", linewidth = 0),
    legend.key.width = unit(0.3, "cm"),        
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(color = 'black', face = "bold", size = 8)
  ) &
  NoAxes()  # 移除坐标轴
dev.off()





