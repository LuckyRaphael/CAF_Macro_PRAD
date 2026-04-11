rm(list = ls())
gc()
setwd('./')
#### CytoTRACE1 ####
rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig3/')
# devtools::install_local("../CytoTRACE_0.3.3.tar.gz")
# install.packages("../monocle.tar", repos = NULL, type = "source")
library(Seurat)
library(CytoTRACE)
library(data.table)
library(tibble)
library(ggplot2)
scRNA <- qs::qread("../results/Fib_anno.qs")
scRNA$cell_subtypes <- scRNA$celltype
exp1 <- as.matrix(scRNA@assays$RNA@counts) #V4
# exp1 <- as.matrix(GetAssayData(scRNA,layer="counts")) #V5
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,] #过滤基因
exp1[1:4,1:4]
results <- CytoTRACE(exp1,ncores = 20)
phenot <- scRNA$cell_subtypes
phenot <- as.character(phenot)
names(phenot) <- rownames(scRNA@meta.data)
emb <- scRNA@reductions[["umap"]]@cell.embeddings
if (!dir.exists('./CytoTRACE1/')) {
  dir.create('./CytoTRACE1/')
}
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './CytoTRACE1/')
dev.off()
plotCytoGenes(results, numOfGenes = 30, outputDir = './CytoTRACE1/')



rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig3/')
library(Seurat)
library(monocle)
library(tidydr)
library(viridis)
# 读取注释好的成纤维细胞亚群
scRNA <- qs::qread("../results/Fib_anno.qs")
scRNA$cell_subtypes <- scRNA$celltype
dim(scRNA)
top50 <- read.csv('../results/CAF_top50markers.csv')
# 自定义颜色
colors<- unique(c('#eebad2','#d5adcf','#c4d691','#d1d2e9',
                  '#e8bf83','#b0c3e0','#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935','#343031',
                  '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31',
                  '#b2dfe8','#667570','#5ec9b8','#314e9b','#616e88',
                  '#c5cdf7','#7092d2','#e69fc4','#b6854d','#8dd1c6','#fffeb3','#bcb9d8','#80b2d4','#feb662',
                  '#b4df65','#d9d9d9','#be7fbf','#cbecc4','#ffee72','#e6bdcb'))
CAFcol<- colors[1:length(unique(scRNA$cell_subtypes))]
dput(unique(scRNA$cell_subtypes))
#将需要的颜色与cell type对应
names(CAFcol) <- c("NDUFA4L2+CAF", "CCL2+CAF", "TAGLN+CAF", "CXCL8+CAF", "RUNX2+CAF", 
                   "APOD+CAF", "POSTN+CAF", "MMP2+CAF", "MCAM+CAF")
####Monocle2成纤维细胞亚群拟时序####
scRNA$celltype <- scRNA$cell_subtypes # 添加一列celltype,便于后续作图
dim(scRNA)
# ## 控制细胞数,内存大可不进行此操作
# set.seed(123456)  # 确保结果可重复
# ## 随机抽取xx个细胞
# a <- sample(1:ncol(scRNA),15000)
# scRNA <- scRNA[,a]
# dim(scRNA)



####第一种####
dat <- Seurat::as.CellDataSet(scRNA)
dat <- estimateSizeFactors(dat)
dat <- detectGenes(dat, min_expr = 10)
fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(dat),num_cells_expressed >= 10))
dat <- reduceDimension(dat,
                       max_components = 2,
                       norm_method = 'log',
                       num_dim = 20,
                       reduction_method = 'tSNE', # c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree")
                       verbose = T,
                       residualModelFormulaStr = "~orig.ident", #去除样本影响
                       check_duplicates=F)
dat <- clusterCells(dat,verbose = F)
clustering_DEG_genes <- differentialGeneTest(dat[expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 20)

# dat <- setOrderingFilter(dat,
#                          ordering_genes = clustering_DEG_genes$gene_short_name[clustering_DEG_genes$use_for_ordering=="TRUE"])
dat <- setOrderingFilter(dat, ordering_genes = unique(top50$gene))
# dat <- setOrderingFilter(dat,
#                          ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000])
dat <- reduceDimension(dat, method = 'DDRTree')
dat <- orderCells(dat)



# ####第二种####
# ### 1. 导入数据，创建对象，预处理
# dat <- Seurat::as.CellDataSet(scRNA)
# #数据预处理，估计size factor和离散度，类似归一化，标准化
# #运行时间较长，可修改核心数
# dat <- estimateSizeFactors(dat)
# dat <- estimateDispersions(dat, cores=30, relative_expr = TRUE)
# dat <- detectGenes(dat, min_expr = 1)#计算每个基因在多少细胞中表达0.1
# 
# ### 2. 用seurat筛选高变基因-常用
# order.genes <- VariableFeatures(FindVariableFeatures(scRNA, assay = "RNA", nfeatures = 1000), assay = "RNA")
# 
# # ### 2. 选择选择排序基因，monocle方法，2000左右
# # disp_table <- dispersionTable(dat)
# # order.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >=
# #                         1 * dispersion_fit) %>% pull(gene_id) %>% as.character() ###0.1可更改
# 
# #标记这些基因，指明基因用于后续的聚类/排序
# dat <- setOrderingFilter(dat,order.genes)
# plot_ordering_genes(dat)
# 
# ### 3. 降维排序，推断轨迹，并按照拟时序给细胞排序
# dat <- reduceDimension(dat, max_components = 2, 
#                        num_dim = 20,
#                        residualModelFormulaStr = "~orig.ident", 
#                        reduction_method = 'DDRTree')
# 
# # residualModelFormulaStr减少其他因素的影响，比如不同样本、不同批次
# dat <- orderCells(dat)#运行时间较长

#############################################################################
# 可以人工设置起始点
## ordering cells by assigning root nodes
GM_state <- function(cds, starting_point, cluster){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds@phenoData@data[,cluster])[,starting_point]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
table(dat$celltype)
root_start = GM_state(cds = dat,starting_point = "POSTN+CAF",cluster = "celltype")
root_start
dat <- monocle::orderCells(dat, root_state = root_start)


# # 查看当前的State分布
# table(pData(dat)$State)
# # 删除State为xx的细胞
# dat <- dat[, pData(dat)$State != 1]
# # 确认删除后的细胞数量
# table(pData(dat)$State)


saveRDS(dat,"../results/CAF_monocle_top50_2.RDS")
dat <- readRDS('../results/CAF_monocle_top50_2.RDS')
# 细胞类型分面拟时序图
p <- plot_cell_trajectory(dat, color_by = "celltype",
                          cell_link_size = 1, theta=0.8,
                          show_backbone=FALSE, cell_size = 0.75)+
  facet_wrap(~celltype, nrow = 3)+
  theme_void()+
  theme(legend.position = "none")+
  scale_color_manual(values=CAFcol)
p
dev.off()
ggsave('./S3A_CAF_CellType_wrap_top50.pdf', p, width = 9, height = 9)

# 细胞类型拟时序图(不分面)
p <- plot_cell_trajectory(dat, color_by = "celltype",
                          cell_link_size = 1, cell_size = 1,
                          show_tree=T,
                          show_backbone=FALSE, theta=0.8)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = CAFcol)
p
dev.off()
ggsave('./Fig3A_CAF_CellType_top50.pdf', p, width=8, height=6)

# CAF亚群随tissue类型分面
p <- plot_cell_trajectory(dat, color_by = "celltype",
                          cell_link_size = 1, theta=0.8,
                          show_backbone=FALSE, cell_size = 1)+
  facet_wrap(~tissue, nrow = 1)+
  theme_void()+
  theme(legend.position = "none")+
  scale_color_manual(values=CAFcol)
p
dev.off()
ggsave('./S3A2_CAF_tissue_wrap_top50.pdf', p, width = 12,height = 5) # useDingbat=F



## 拟时序pseudotime图
# dat$Pseudotime <- max(dat$Pseudotime) - dat$Pseudotime # 使拟时序方向相反
p <- plot_cell_trajectory(dat, color_by = "Pseudotime", cell_size = 1,
                          show_backbone=FALSE, cell_link_size = 1
)+
  scale_colour_gradientn(
    colours = (viridis(100)),
    # colours = rev(colorRampPalette(c('#fcffd4','#52bbc1','#0e206a'))(100)),
    # colours = colorRampPalette(c("#4DAC26", "white", "#D01C8B"))(100),
    # colours = colorRampPalette(rev(brewer.pal(9, "PRGn")))(100),
    # colours = colorRampPalette(c("navy","white","firebrick3"))(100),
    na.value = "transparent",
    name = "Pseudotime" ,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p
dev.off()
ggsave('./Fig3B_CAF_Pseudotime_top50.pdf', p, width = 8, height = 6) # useDingbat=F



# State
p<-plot_cell_trajectory(dat,color_by = "State",
                        cell_size = 1, cell_link_size = 1, theta=0.8,
                        show_tree=T, show_backbone=FALSE)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values=jjAnno::useMyCol('calm',n=7))
p
dev.off()
ggsave('./Fig3C_CAF_State_top50.pdf', p, width=8, height=6)

## State饼图
length(unique(dat$State))
tmp=sapply(1:length(unique(dat$State)),function(x){
  sapply(levels(factor(dat$celltype)),function(y){
    length(which(dat$celltype==y & dat$State==x))
  })
})
colnames(tmp)=paste0('State',1:length(unique(dat$State)))

plot_lists=lapply(colnames(tmp),function(x){
  df=data.frame(
    value=as.numeric(as.vector(tmp[,x])),
    group=factor(rownames(tmp),levels=levels(factor(dat$celltype)))
  )
  ggplot(df, aes(x="",y=value,fill=group, show.legend = FALSE))+
    geom_bar(stat="identity",width = 1)+
    scale_fill_manual(values=CAFcol)+labs(title=x)+
    coord_polar("y")+theme_minimal()+theme(
      axis.text = element_blank(),  # 移除所有坐标轴文本
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"))+Seurat::NoLegend()
  
})
pdf('./Fig3C_Monocle2_State_CellType_pie.pdf')
cowplot::plot_grid(plotlist=plot_lists,nrow=1)
dev.off()


