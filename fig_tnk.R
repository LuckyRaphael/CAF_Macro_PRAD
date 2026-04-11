################################################################################
# Section 1: Environment Setup and Data Loading
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_tnk/')

library(Seurat)
library(data.table)
library(tibble)
library(ggplot2)
library(tidydr)

# Load annotated T/NK cell Seurat object
scRNA <- qs::qread("../results/T_NK_cell.qs")

################################################################################
# Section 2: Supplemental Figure S8a - T/NK Cell Subtype UMAP
################################################################################

# Set identity to celltype
Idents(scRNA) <- 'celltype'

# Define custom color palette for T/NK sub-clusters
colors <- c(
  '#cf5d51','#dd8349','#ed9baf','#2d6089','#4d85a0','#a1d0d3',
  '#c1d7e8','#6dadc3','#698a5c','#b2d081','#e1ba32','#5c2d8b',
  '#187e8c','#64abc9','#ffca56','#de3935','#343031'
)

# Map colors to the 14 identified clusters
maincol <- colors[1:length(unique(scRNA$celltype))]
names(maincol) <- c(
  "CD8_C3_ISG15", "CD4_C3_FOXP3", "CD4_C4_CD40LG", "CD4_C1_CCR7", 
  "CD8_C1_GZMK", "CD8_C4_NDUFAF8", "Prolif_STMN1", "NK_C1_FGFBP2", 
  "CD8_C2_CCL4", "NK_C2_TYROBP", "CD4_C6_CLU", "CD4_C7_CXCR4", 
  "CD4_C5_CXCL13", "CD4_C2_TCF7"
)

# Plot UMAP
pdf(file = "./S8a_dimplot_T.pdf", width = 8, height = 7)
DimPlot(scRNA, group.by = 'celltype', 
        cols = maincol, 
        label = FALSE, 
        pt.size = 1.5, 
        raster = TRUE, 
        reduction = "umap") +
  NoAxes() + 
  labs(x = "UMAP1", y = "UMAP2", title = "T/NK cell") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

################################################################################
# Section 3: Supplemental Figure S8b - Canonical Marker FeaturePlot
################################################################################

# Markers: CD3E (T-cell), CD4 (Helper), CD8A (Cytotoxic), GNLY (NK/Cytotoxicity)
select_markers <- c('CD3E','CD4','CD8A','GNLY')

# Define Red-themed gradient for expression
# Using gray for low expression and a deep red (#7c2124) for high expression
mycol <- c("#bebebe", "#7c2124")

pdf(file = "./S8b_featureplot_T.pdf", width = 8, height = 7)
FeaturePlot(scRNA, 
            order = TRUE, 
            features = select_markers, 
            cols = mycol, 
            raster = TRUE, 
            ncol = 2, 
            pt.size = 2, 
            reduction = "umap") & 
  NoLegend() & 
  NoAxes()
dev.off()
