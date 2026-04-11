################################################################################
# Section 1: Environment Setup & Data Loading
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig4/')

library(Seurat)
library(scRNAtoolVis)
library(patchwork)
library(tidydr)
library(RColorBrewer)
library(magrittr)
library(ggplot2)

# Load annotated Myeloid cell sub-clusters
scRNA <- qs::qread("../results/Mye_anno.qs")
scRNA$cell_subtypes <- scRNA$celltype

################################################################################
# Section 2: Myeloid Sub-cluster Visualization (UMAP)
################################################################################

# Define comprehensive color palette
cols <- unique(c(
  '#93bdd1','#ea7e38','#edaf74','#2190b4','#add26b',
  '#c8e3c0','#efec96','#bc917d','#d6a6a5','#d8585c',
  '#f5d2df','#e88e90','#8dd1c6','#fffeb3','#bcb9d8',
  '#80b2d4','#feb662','#b4df65','#d9d9d9','#be7fbf',
  '#cbecc4','#ffee72','#e6bdcb','#306311','#df772b',
  '#e184ca','#de554c','#e8e684','#b2dfe8','#667570',
  '#5ec9b8','#314e9b','#616e88','#7623f1','#eaeaea'
))

Myecol <- cols[1:length(unique(scRNA$cell_subtypes))]

# Map colors to specific cell types
names(Myecol) <- c(
  "Macro_APOE", "Macro_CD83", "cDC2_FCER1A", "Prolif_MKI67", 
  "Macro_CXCL10", "Mono_FCGR3A", "Mono_S100A12", "Neutro_FCGR3B", 
  "cDC2_CD1C", "Neutro_MPO", "Macro_CCL3"
)

# Plot Figure 5B: Myeloid UMAP
pdf(file = "./Fig5B_dimplot.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'cell_subtypes', 
        cols = Myecol, pt.size = 2.5, raster = TRUE, reduction = "umap") +
  NoAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank())
dev.off()

################################################################################
# Section 3: Marker Gene Visualization (DotPlot)
################################################################################

Idents(scRNA) <- 'cell_subtypes'

# Define canonical marker gene sets
featureSets <- list(
  Macro_APOE    = c('APOC1', 'GPNMB', 'LGMN', 'SLC40A1', 'APOE'),
  Macro_CD83    = c('BAG3', 'IER5', 'RGS1', 'IL1B', 'CD83'),
  Mono_S100A12  = c('S100A8', 'S100A9', 'VCAN', 'EREG', 'S100A12'),
  cDC2_FCER1A   = c('BIRC3', 'PPA1', 'TXN', 'CPVL', 'FCER1A'),
  Neutro_FCGR3B = c('CXCR2', 'CSF3R', 'FPR1', 'TREM1', 'FCGR3B'),
  Mono_FCGR3A   = c('SMIM25', 'LILRB2', 'CDKN1C', 'LST1', 'FCGR3A'),
  Prolif_MKI67  = c('TUBA1B', 'STMN1', 'TOP2A', 'CENPF', 'MKI67'),
  Neutro_MPO    = c('ELANE', 'PRTN3', 'AZU1', 'CTSG', 'MPO'),
  Macro_CXCL10  = c('ISG15', 'GBP1', 'MX1', 'GBP4', 'CXCL10'),
  cDC2_CD1C     = c('ARL4C', 'EZR', 'JAML', 'RUNX3', 'CD1C'),
  Macro_CCL3    = c('GABARAP', 'TYMP', 'UCP2', 'PRELID1', 'CCL3')
)

all_genes <- unlist(featureSets)

# Generate DotPlot and extract data for custom ggplot
p_dot <- DotPlot(scRNA, features = featureSets) + RotatedAxis()

result_data <- p_dot$data %>%
  dplyr::select(feature.groups, features.plot, id, avg.exp.scaled, pct.exp) %>%
  dplyr::rename(Function = feature.groups, Gene = features.plot, Group = id, 
                AvgExpr = avg.exp.scaled, PctExpr = pct.exp) %>%
  dplyr::mutate(
    Function = factor(Function, levels = names(featureSets)),
    Gene = factor(Gene, levels = rev(all_genes)),
    Group = factor(Group, levels = names(featureSets))
  )

# Plot custom formatted DotPlot
pdf(file = "./Mye_markergene_dotplot.pdf", width = 12, height = 5)
ggplot(result_data, aes(y = Group, x = Gene, color = AvgExpr, size = PctExpr)) +
  geom_point() +
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "PiYG")), # Pink-Green divergent palette
    guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"),
    name = "Relative\nExpression"
  ) +
  scale_size_continuous(name = "Percentage\nExpressed", range = c(0, 6)) +
  theme_minimal() +
  xlab("") + ylab("") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = "black", vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "#35A79D", linewidth = 1),
    panel.background = element_rect(fill = "#F1F6FC"),
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white")
  )
dev.off()

################################################################################
# Section 4: Supplemental FeaturePlots (S7A)
################################################################################

JIC_markers <- c(
  "CD163","MRC1", "CCR5",       # Macrophage
  "CD14","FCN1",               # Monocyte
  "IDO1","CLEC9A",             # cDC1
  "CD1C","CLEC10A",            # cDC2
  "LILRA4","JCHAIN",           # pDC
  "LAMP3",                     # mDC
  "CSF3R","S100A8","S100A9",   # Neutrophils
  "TOP2A","MKI67","STMN1"      # Proliferation
)

pdf(file = "./S7A_featureplot.pdf", width = 15, height = 9)
FeaturePlot(scRNA, raster = TRUE, order = TRUE, 
            features = JIC_markers, 
            cols = c("#e0f2f7", "#b3d8e6", "#80bed9", "#4aa3c8", "#1f8abf", "#1570a6", "#0f5688"),
            pt.size = 3, reduction = "umap", ncol = 6) & 
  theme_void() & 
  theme(
    legend.frame = element_rect(colour = "black"),        
    legend.ticks = element_blank(),
    legend.key.width = unit(0.3, "cm"),        
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(color = 'black', face = "bold", size = 8)
  ) &
  NoAxes()
dev.off()


