################################################################################
# Section 1: Environment Setup and Data Loading
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig7/')

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidydr)

# Load annotated pan-cancer single-cell object
scRNA <- qs::qread("/home/w282308/JIC_CRC/data/pan/scRNA_pancancer_anno.qs")

# Data Overview
dim(scRNA)
table(scRNA$Cancer_type)

# Calculate sample-level statistics per Dataset
sample_stats <- scRNA@meta.data %>%
  distinct(Sample, .keep_all = TRUE) %>%
  group_by(Datasets) %>%
  summarise(
    Total_Samples = n(),
    Normal_Samples = sum(tissue == "Normal"),
    Tumor_Samples = sum(tissue == "Tumor")
  )

print(sample_stats)

################################################################################
# Section 2: Color Palette and Metadata Mapping
################################################################################

# Custom master color palette
cols <- c(
  '#005f73','#099396','#ca6702','#ed9b00','#a12c2f','#bb3f02',
  '#b15d2e','#e82036','#237bba','#40ac53','#f1c742','#7498b5',
  '#e7e7e7','#d66969','#e5f499','#3188bd','#fee18b','#67c1a6',
  '#9e0242','#d6404e','#fcae62','#6050a5','#aadda3','#DC050C'
)

# Map colors to major cell types
Maincol <- cols[1:length(unique(scRNA$celltype))]
names(Maincol) <- c(
  "T/NK cells", "Endothelial cells", "Myeloid", "Prolif", "B cells", 
  "Plasma cells", "Mast cells", "Epithelial cells", "Fibroblasts"
)

################################################################################
# Section 3: Dimensional Reduction Visualization (UMAP)
################################################################################

# 1. UMAP Plot by Cell Type
pdf(file = "./pancancer_dimplotplot.pdf", width = 7, height = 5)
DimPlot(scRNA, reduction = "umap", label = TRUE, 
        group.by = 'celltype', cols = Maincol, 
        raster = TRUE, pt.size = 1) & 
  NoAxes() & 
  labs(x = "UMAP1", y = "UMAP2", title = "") & 
  tidydr::theme_dr(
    xlength = 0.2, 
    ylength = 0.2, 
    arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")
  ) & 
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

# 2. UMAP Plot by Dataset (Batch/Source)
pdf(file = "./pancancer_dimplotplot_Datasets.pdf", width = 7, height = 5)
DimPlot(scRNA, reduction = "umap", label = FALSE, 
        group.by = 'Datasets', cols = cols[6:22], 
        raster = TRUE, pt.size = 1) & 
  NoAxes() & 
  labs(x = "UMAP1", y = "UMAP2", title = "") & 
  tidydr::theme_dr(
    xlength = 0.2, 
    ylength = 0.2, 
    arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")
  ) & 
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()
