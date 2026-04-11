rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_epi/')
scRNA <- qs::qread("../results/Epi_anno.qs")
# scRNA$celltype <- scRNA$ct.sub.epi
# table(scRNA$celltype)
# scRNA$celltype <- factor(as.character(scRNA$celltype))
# table(scRNA$celltype)
# qs::qsave(scRNA,"../results/Epi_anno.qs")
colnames(scRNA@meta.data)
DimPlot(scRNA, group.by = "celltype", reduction = "umap", label = T)
dev.off()
scRNA$cell_subtypes <- scRNA$celltype
table(scRNA$cell_subtypes)
Epicol<- unique(c('#f79f7c','#9eafd3','#e58dbc','#a8d156','#be7042',
                  '#b53365','#e68122','#efb51b','#6fa8d0','#117e7a','#8660b1',
                  '#d71921','#232a6a','#1a8a40','#892390','#f27e26'))[1:length(unique(scRNA$cell_subtypes))]
dput(as.character(unique(scRNA$cell_subtypes)))
#将需要的颜色与cell type对应
names(Epicol) <- c("Luminal", "Club", "Basal", "LPCs", "Hillock", "NE")
# 小提琴图箱线图
library(SCP)
pdf("./cnv_epiclusters_boxplot.pdf", width=8, height=6)
FeatureStatPlot(
  srt = scRNA, group.by = "cell_subtypes",
  palcolor = Epicol, 
  bg_palcolor = Epicol,
  ncol = 1,
  multiplegroup_comparisons = TRUE,
  multiple_method = "kruskal.test",
  # sig_label = 'p.format',
  stat.by = 'totalCNV', add_box = TRUE,
  bg.by = 'cell_subtypes'
) & xlab('') & 
  ylab('totalCNV') &
  NoLegend()
dev.off()

table(scRNA$tissue)
scRNA$tissue <- factor(scRNA$tissue, levels = c("Primary", "CRPC", "mCRPC"))
groupcol <- c("#92C5DE","#F4A582","darkred") # c( "#8DD3C7","#BC80BD","lightblue")
names(groupcol) <- c("Primary", "CRPC", "mCRPC")
pdf("./cnv_epitissue_boxplot.pdf", width=6, height=6)
FeatureStatPlot(
  srt = scRNA, group.by = "tissue",
  palcolor = groupcol, 
  bg_palcolor = groupcol,
  ncol = 1,
  multiplegroup_comparisons = TRUE,
  multiple_method = "kruskal.test",
  # sig_label = 'p.format',
  stat.by = 'totalCNV', add_box = TRUE,
  bg.by = 'tissue'
) & xlab('') & 
  ylab('totalCNV') &
  NoLegend()
dev.off()


pdf(file = "./epi_tissue_umap.pdf", width = 6, height = 5)
DimPlot(scRNA, label = FALSE, group.by = 'tissue', #!!!
        cols = groupcol, # 使用预定义颜色向量
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



select_gene <- 'totalCNV'
p <- ggdraw() &
  draw_plot(FeaturePlot(scRNA, raster=T, order=T, 
                        # raster.dpi=c(2048,2048), 
                        min.cutoff = "q10",   # 最小值
                        max.cutoff = 3000,
                        features=select_gene, col=c( "#ADD8E6", "#6A0DAD"),
                        pt.size=2, reduction="umap", 
                        ncol=1)&
              NoAxes(),scale = 0.9)&
  #绘制箭头
  draw_plot(ggplot(tibble(group = c("UMAP1", "UMAP2"),
                          x = c(0, 0), xend = c(1, 0),
                          y = c(0, 0), yend = c(0, 1),
                          lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                          angle = c(0, 90))) +
              geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
                           arrow = arrow(angle = 20, type = "closed", length = unit(0.05, "npc")),
                           size = 0.8, lineend = "round") +
              geom_text(aes(lx, ly, label = group, angle = angle), size = 4) +
              theme_void() +
              coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1)),
            x = 0.05,y = 0.05,width = 0.2,height = 0.2)
p
dev.off()
ggsave('./epi_cnv_featureplot.pdf', p, width=7, height=6)



