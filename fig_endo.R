################################################################################
# Section 1: Environment Setup and Data Loading
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig_endo/')

library(Seurat)
library(data.table)
library(tibble)
library(ggplot2)
library(tidydr) # Required for theme_dr

# Load annotated Endothelial cell Seurat object
scRNA <- qs::qread("../results/Endo_anno.qs")

################################################################################
# Section 2: Supplemental Figure S9a - Endothelial Subtype UMAP
################################################################################

# Set active identity to celltype for consistency
Idents(scRNA) <- 'celltype'

# Define a distinct color palette for EC subtypes
# Optimized to ensure clear differentiation between Venous, Capillary, and Arterial ECs
colors <- c(
  '#5c2d8b', '#187e8c', '#64abc9', '#ffca56', 
  '#de3935', '#343031', '#424b00', '#376696'
)

# Map colors to specific EC sub-clusters
maincol <- colors[1:length(unique(scRNA$celltype))]
names(maincol) <- c(
  "VenECs_KLK3", "CapECs_RGCC", "ArtECs_ENPP2", "VenECs_CXCL8", 
  "CapECs_ESM1", "VenECs_NOP53", "VenECs_SELE", "ECs_FABP4"
)

# Generate UMAP plot
pdf(file = "./S9a_dimplot_Endo.pdf", width = 8, height = 7)
DimPlot(
  scRNA, 
  label = TRUE, 
  group.by = 'celltype', 
  cols = maincol, 
  pt.size = 1, 
  raster = FALSE, 
  reduction = "umap"
) +
  NoAxes() + 
  labs(
    x = "UMAP1", 
    y = "UMAP2", 
    title = "Endothelial Cell Clusters"
  ) +
  # Custom axis arrows for a cleaner aesthetic
  theme_dr(
    xlength = 0.2, 
    ylength = 0.2, 
    arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    aspect.ratio = 1, 
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()
