rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig1/')
library(Seurat)
library(tidydr)
## 读取单细胞数据
scRNA <- qs::qread("../results/StepF.All_Cells_Seurat_Object.qs")
dim(scRNA)
colnames(scRNA@meta.data)
table(scRNA$ct.sub)
table(scRNA$group)
table(scRNA$data.sets)
table(scRNA$group2)
table(scRNA$orig.ident)
table(scRNA$sampleID)
table(scRNA$patientID)
length(unique(scRNA$orig.ident))
length(unique(scRNA$sampleID))
length(unique(scRNA$patientID))
scRNA$Sample <- scRNA$sampleID
scRNA$tissue <- scRNA$group
scRNA$Datasets <- scRNA$data.sets
scRNA$celltype <- scRNA$ct.sub
#写出合并的单细胞数据
qs::qsave(scRNA,"../results/StepF.All_Cells_Seurat_Object.qs")


# DimPlot(scRNA, group.by = "ct.sub", label = T, pt.size = 2, raster = T)&  # + ggsci::scale_color_lancet()
#   theme_dr(arrow = grid::arrow(length = unit(0, "inches")))&
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# dev.off()
# 
# 
# DimPlot(scRNA, group.by = "data.sets", label = T, pt.size = 2, raster = T)&  # + ggsci::scale_color_lancet()
#   theme_dr(arrow = grid::arrow(length = unit(0, "inches")))&
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# dev.off()



rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig1/')
####获取细胞大群的marker基因便于空转打分####
scRNA <- qs::qread("../results/StepF.All_Cells_Seurat_Object.qs")
dim(scRNA)
length(unique(scRNA$Sample))
future::plan("multicore", workers = 10)
library(dplyr)
library(data.table)
colnames(scRNA@meta.data)
scRNA@active.ident<-factor(scRNA$celltype)
badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),
                       x=rownames(scRNA), perl=T, value=T))
markers <- FindAllMarkers(
  subset(x = scRNA, downsample = 500),
  # scRNA,
  features=rownames(scRNA)[-which(rownames(scRNA) %in% badgenelt)],
  only.pos = TRUE, 
  min.pct = 0.25, logfc.threshold = 0.25)
markers$pct_diff<-markers$pct.1-markers$pct.2
# 获取每个细胞亚群Top50marker基因
markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) ->top50 
# 保存全部marker基因
saveRDS(markers,"../results/Maintype_CellType_diffgenes.RDS") 
fwrite(markers , "../results/Maintype_allmarkers.csv")
# 保存Top50 Marker基因
fwrite(top50 , "../results/Maintype_top50markers.csv")




table(scRNA$celltype)
####单细胞大群注释完成####
# 直接保存髓系跟成纤维细胞亚群
scRNA_subset <- subset(scRNA, subset = celltype == "Epithelia")
qs::qsave(scRNA_subset,"../results/scRNA_epi_raw.qs")
scRNA_subset <- subset(scRNA, subset = celltype == "Fib")
qs::qsave(scRNA_subset,"../results/scRNA_fibro_raw.qs")
scRNA_subset <- subset(scRNA, subset = celltype == "Endo")
qs::qsave(scRNA_subset,"../results/scRNA_endo_raw.qs")
scRNA_subset <- subset(scRNA, subset = celltype == "Mye")
qs::qsave(scRNA_subset,"../results/scRNA_mye_raw.qs")
scRNA_subset <- subset(scRNA, subset = celltype == "B")
qs::qsave(scRNA_subset,"../results/scRNA_B_raw.qs")
scRNA_subset <- subset(scRNA, subset = celltype == "T")
qs::qsave(scRNA_subset,"../results/scRNA_T_raw.qs")














####各细胞类型亚群注释结束####
rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig1/')
####细胞注释完成后绘制Fig1####
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
library(scRNAtoolVis)
library(RColorBrewer)
scRNA <- qs::qread("../results/StepF.All_Cells_Seurat_Object.qs")
colnames(scRNA@meta.data)
dim(scRNA)
table(scRNA$Datasets)
length(unique(scRNA$Sample))


# 聚类
# scRNA <- FindNeighbors(scRNA, dims = 20, reduction = 'harmony') 
# scRNA <- FindClusters(scRNA, resolution = 0.4)
# table(scRNA@meta.data$seurat_clusters)
# 自定义颜色
colors<- unique(c('#639dd0','#e1959a','#ffc79a','#7cc1f3','#e1eddc',
                  '#e84a34','#4dbbd6','#3c5487','#f29b80','#8591b4','#808080',
                  '#f09c00','#016170','#b11f13','#016170','#a7d3b4','#e6daa5',
                  '#005f73','#099396','#ca6702','#ed9b00','#a12c2f','#bb3f02',
                  '#b15d2e','#e82036','#237bba','#40ac53',
                  '#f1c742','#7498b5','#e7e7e7','#d66969',
                  '#e5f499','#3188bd','#fee18b','#67c1a6',
                  '#9e0242','#d6404e','#fcae62','#6050a5','#aadda3','#e5f499','#3188bd','#fee18b','#67c2a5',
                  '#00b4a9','#049964','#0968b3','#04c0ec','#94288e','#f47a1f','#a72b24','#f59f83','#eb1275',
                  '#e8733d','#981118','#000000','#e9c56d','#006a81','#199169',
                  '#d37db5','#ad60a5','#7e2e92','#e775ad','#01afb2','#86d3e4','#4c92ce','#a0c3e7','#6d65ae',
                  '#850663','#ffbde5','#7860f5','#2000fc','#9dfffd','#01afb2','#46da44','#b8712e','#fe691e',
                  '#c5cdf7','#7092d2','#e69fc4','#b6854d','#8dd1c6','#fffeb3','#bcb9d8','#80b2d4','#feb662',
                  '#b4df65','#d9d9d9','#be7fbf','#cbecc4','#ffee72','#e6bdcb'))


pdf(file = "./dimplotplot_datasets.pdf", width = 8, height = 6)
## 绘制seurat_clusters的umap图
DimPlot(scRNA, label = FALSE, group.by = 'Datasets', #!!!
        cols = colors, # 使用预定义颜色向量
        pt.size = 1.5,
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
        panel.grid = element_blank())
dev.off()




tissuecol<- c( '#9e0242','#fcae62','#d6404e','#6050a5','#a72b24','#aadda3','#e5f499','#3188bd','#fee18b','#67c2a5',
               '#00b4a9','#049964','#0968b3','#04c0ec','#94288e','#f47a1f','#f59f83','#eb1275')
tissuecol<- tissuecol[1:length(unique(scRNA$tissue))] #自定义大群颜色
dput(as.character(unique(scRNA$tissue)))
#将需要的颜色与cell type对应
names(tissuecol) <- c("Primary", "mLN", "CRPC", "Normal", "mCRPC", "LN")
pdf(file = "./dimplotplot_tissue.pdf", width = 8, height = 6)
## 绘制seurat_clusters的umap图
DimPlot(scRNA, label = FALSE, group.by = 'tissue', #!!!
        cols = tissuecol, # 使用预定义颜色向量
        pt.size = 1.5,
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
        panel.grid = element_blank())
dev.off()






colnames(scRNA@meta.data)
####比例饼图####
# 绘制数据集细胞数量比例饼图
library(RColorBrewer)
mynames <- table(scRNA$Datasets) %>% names()
myratio <- table(scRNA$Datasets) %>% as.numeric()
pielabel <- paste0(mynames,"(",round(myratio/sum(myratio)*100,2),"%)")
colorS <- brewer.pal(12,"Set3")
pdf(file = "./S1B_pieplot_datasets_cellnumber.pdf", width = 13, height = 9)
pie(myratio,
    labels = pielabel,
    radius =1, clockwise=F,
    main = "Cell number",  # 添加标题
    border = "white",
    col = colors)
dev.off()


# 绘制数据集样本数量比例饼图
# 统计每个数据集的唯一样本数（Sample列是样本ID）
library(dplyr)
library(RColorBrewer)
dataset_sample_counts <- scRNA@meta.data %>%
  distinct(Datasets, Sample) %>%  # 去重，确保每个样本只计一次
  dplyr::count(Datasets, name = "sample_count")  # 计算每个数据集的样本数

# 提取数据
mynames <- dataset_sample_counts$Datasets
myratio <- dataset_sample_counts$sample_count
pielabel <- paste0(mynames, " (", round(myratio/sum(myratio)*100, 2), "%)")  # 仅显示百分比
# 设置颜色
colorS <- brewer.pal(length(mynames), "Set3")  # 自动匹配颜色数量
pdf(file = "./S1C_pieplot_datasets_samplenumber.pdf", width = 13, height = 9)
pie(myratio, 
    labels = pielabel,
    radius = 1,
    clockwise = FALSE,
    main = "Sample number",  # 标题
    border = "white",
    col = colors)
dev.off()




####绘制Figure1B大群注释UMAP图####
library(scRNAtoolVis)
table(scRNA$celltype)
colors<- unique(c('#da1921','#5c2d8b','#187e8c','#64abc9','#ffca56','#343031',
                  '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31',
                  '#e54a32','#00a087','#f29b80','#91d1c1','#4cbbd6','#395388','#8491b5',
                  '#b2dfe8','#667570','#5ec9b8','#314e9b','#616e88',
                  '#c5cdf7','#7092d2','#e69fc4','#b6854d','#8dd1c6','#fffeb3','#bcb9d8','#80b2d4','#feb662',
                  '#b4df65','#d9d9d9','#be7fbf','#cbecc4','#ffee72','#e6bdcb'))
# ??clusterCornerAxes
Maincol<- colors[1:length(unique(scRNA$celltype))] #自定义大群颜色
dput(as.character(unique(scRNA$celltype)))
#将需要的颜色与cell type对应
names(Maincol) <- c("Epithelia", "T", "Mye", "Fib", "Endo", "Mast", "B")

pdf(file = "./dimplotplot_celltype.pdf", width = 8, height = 6)
## 绘制seurat_clusters的umap图
DimPlot(scRNA, label = FALSE, group.by = 'celltype', #!!!
        cols = Maincol, # 使用预定义颜色向量
        pt.size = 1.5,
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
        panel.grid = element_blank())
dev.off()





