rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_b/')
#--------------------------B细胞亚群分析--------------------------
library(Seurat)
library(data.table)
library(tibble)
library(ggplot2)
## 读取B细胞注释结果
scRNA <- qs::qread("../results/B_cell.qs") # 读取单细胞数据
#--------------------------S10a B细胞亚群umap图--------------------------
dim(scRNA)
length(unique(scRNA$orig.ident))
length(unique(scRNA$celltype))

## 首先设置细胞身份标识
Idents(scRNA) <- 'celltype' #!!!
table(scRNA$celltype)

## 定义颜色向量
colors <- c('#eebad2','#d5adcf','#c4d691','#d1d2e9',
            '#e8bf83','#b0c3e0','#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935','#343031',
            '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31',"#9E9E9E","#D64550",
            "#FF977E","#A43B76",
            "#9A64A0","#750985",
            "#5ECBC8", "#28788D", "#37A794", 
            "#4C5D8A", "#039BE5", "#BF360C"
)

maincol <- colors[1:length(unique(scRNA$celltype))]

dput(as.character(unique(scRNA$celltype)))
# 将需要的颜色与celltype对应
names(maincol) <- c("IgG_Plasma", "Naive_B", "Memory_B", "IgA_Plasma", "GC_B", 
                    "Prolif_B")
pdf(file = "./dimplot_Bsubtype.pdf", width = 8, height = 7)
DimPlot(scRNA,label = TRUE,group.by = 'celltype', #!!!
        cols = maincol, # 使用预定义颜色向量
        pt.size = 1,
        raster = FALSE,
        reduction = "umap") +
  NoAxes()  +           # 隐藏默认坐标轴标签和刻度线
  labs(x = "UMAP1", y = "UMAP2", 
       title = "B cell") +  # 修改labs标签
  theme_dr(xlength = 0.2, # x轴箭头的长度比例
           ylength = 0.2, # y轴箭头的长度比例
           arrow = grid::arrow(length = unit(0.1, "inches"), # 箭头尺寸：0.1英寸
                               type = "closed")) +  # 箭头类型：实心闭合箭头
  # NoLegend() +  # 移除图例
  theme(aspect.ratio = 1,  # 设置图形宽高比为1:1（正方形）
        panel.grid = element_blank())
dev.off()



#--------------------------S10b B细胞亚群在不同分组间的比例饼图--------------------------
# 加载必要的包
library(ggplot2)
library(dplyr)
# 绘制数据集细胞数量比例饼图
library(RColorBrewer)
colnames(scRNA@meta.data)
scRNA$tissue[is.na(scRNA$tissue)] <- "NA"
scRNA$group <- scRNA$tissue
table(scRNA$group)
# 更简洁的版本
groups <- unique(scRNA@meta.data$group)
groups <- groups[groups != "NA"]
groups <- groups[groups != "Noinfo"]

for (group_name in groups) {
  # 筛选数据
  group_data <- scRNA@meta.data %>% filter(group == group_name)
  celltype_table <- table(group_data$celltype)
  
  mynames <- names(celltype_table)
  myratio <- as.numeric(celltype_table)
  pielabel <- paste0(mynames, "(", round(myratio/sum(myratio)*100, 2), "%)")
  current_maincol <- maincol[mynames]
  
  # 保存PDF
  pdf(paste0("S10b_pie_chart_", group_name, ".pdf"), width = 9, height = 7)
  pie(myratio, labels = pielabel, radius = 1, clockwise = FALSE,
      main = paste("", group_name), 
      border = "white", col = current_maincol)
  dev.off()
}





#--------------------------S10C B细胞亚群在不同样本中的桑基比例柱状图--------------------------
## scRNAtoolVis
library(scRNAtoolVis)
colnames(scRNA@meta.data)
p <- cellRatioPlot(object = scRNA,
                   sample.name = "tissue", #!!!
                   celltype.name = "celltype", 
                   fill.col = maincol, # 修改颜色!!!
                   flow.alpha = 0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 倾斜45度并右对齐
p
ggsave(filename = "./s10c2_tissueratio_barplot.pdf", p, width = 6, height = 6, dpi = 300) 
dev.off()






