################################################################################
# Section 1: Data Loading & Metadata Setup
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_epi/')

library(Seurat)
library(tidydr)
library(SCP)
library(ggplot2)
library(cowplot)
library(dplyr)

# Load annotated epithelial Seurat object
scRNA <- qs::qread("../results/Epi_anno.qs")

# Inspect metadata and initial UMAP
colnames(scRNA@meta.data)
# DimPlot(scRNA, group.by = "celltype", reduction = "umap", label = T)

# Standardize subtype column
scRNA$cell_subtypes <- scRNA$celltype
table(scRNA$cell_subtypes)

# Define color palette for epithelial subtypes
Epicol <- unique(c('#f79f7c','#9eafd3','#e58dbc','#a8d156','#be7042',
                   '#b53365','#e68122','#efb51b','#6fa8d0','#117e7a','#8660b1',
                   '#d71921','#232a6a','#1a8a40','#892390','#f27e26'))[1:length(unique(scRNA$cell_subtypes))]

# Map colors to specific cell types
names(Epicol) <- c("Luminal", "Club", "Basal", "LPCs", "Hillock", "NE")

################################################################################
# Section 2: CNV Statistical Analysis (Boxplots/Violin Plots)
################################################################################

# 1. totalCNV by Cell Subtype
pdf("./cnv_epiclusters_boxplot.pdf", width=8, height=6)
FeatureStatPlot(
  srt = scRNA, 
  group.by = "cell_subtypes",
  palcolor = Epicol, 
  bg_palcolor = Epicol,
  ncol = 1,
  multiplegroup_comparisons = TRUE,
  multiple_method = "kruskal.test",
  stat.by = 'totalCNV', 
  add_box = TRUE,
  bg.by = 'cell_subtypes'
) & xlab('') & ylab('totalCNV') & NoLegend()
dev.off()

# 2. totalCNV by Tissue Type
# Set factor levels for clinical progression
scRNA$tissue <- factor(scRNA$tissue, levels = c("Primary", "CRPC", "mCRPC"))
groupcol <- c("#92C5DE", "#F4A582", "darkred")
names(groupcol) <- c("Primary", "CRPC", "mCRPC")

pdf("./cnv_epitissue_boxplot.pdf", width=6, height=6)
FeatureStatPlot(
  srt = scRNA, 
  group.by = "tissue",
  palcolor = groupcol, 
  bg_palcolor = groupcol,
  ncol = 1,
  multiplegroup_comparisons = TRUE,
  multiple_method = "kruskal.test",
  stat.by = 'totalCNV', 
  add_box = TRUE,
  bg.by = 'tissue'
) & xlab('') & ylab('totalCNV') & NoLegend()
dev.off()

################################################################################
# Section 3: UMAP Visualizations
################################################################################

# UMAP colored by Tissue
pdf(file = "./epi_tissue_umap.pdf", width = 6, height = 5)
DimPlot(scRNA, label = FALSE, group.by = 'tissue', 
        cols = groupcol, pt.size = 1.5, raster = TRUE, reduction = "umap") +
  NoAxes() + # Hide default axes and ticks
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

# FeaturePlot for totalCNV with custom corner axes
select_gene <- 'totalCNV'

p <- ggdraw() &
  draw_plot(
    FeaturePlot(scRNA, raster=T, order=T, 
                min.cutoff = "q10", 
                max.cutoff = 3000,
                features = select_gene, 
                cols = c("#ADD8E6", "#6A0DAD"), # Light blue to dark purple
                pt.size = 2, reduction = "umap") & 
      NoAxes(), 
    scale = 0.9
  ) &
  # Add custom UMAP arrow axes
  draw_plot(
    ggplot(tibble(group = c("UMAP1", "UMAP2"),
                  x = c(0, 0), xend = c(1, 0),
                  y = c(0, 0), yend = c(0, 1),
                  lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                  angle = c(0, 90))) +
      geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
                   arrow = arrow(angle = 20, type = "closed", length = unit(0.05, "npc")),
                   linewidth = 0.8, lineend = "round") +
      geom_text(aes(lx, ly, label = group, angle = angle), size = 4) +
      theme_void() +
      coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1)),
    x = 0.05, y = 0.05, width = 0.2, height = 0.2
  )

ggsave('./epi_cnv_featureplot.pdf', p, width=7, height=6)
