################################################################################
# Section 1: Environment Setup and Data Loading
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig2/')

library(Seurat)
library(scRNAtoolVis)
library(patchwork)
library(ggplot2)
library(reshape2)

# Load annotated Fibroblast subtype data
scRNA <- qs::qread("../results/Fib_anno.qs")
table(scRNA$celltype)

# Synchronize subtype naming
scRNA$cell_subtypes <- scRNA$celltype

################################################################################
# Section 2: Dimensional Reduction (UMAP Tuning)
################################################################################

# Re-run UMAP for visualization optimization
# Note: higher n.neighbors = more global/continuous, lower = more local/discrete
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:10,
                 n.neighbors = 40, min.dist = 0.5,
                 reduction.name = "umap")

# Define primary cluster color palette
cluster_colors <- unique(c('#0073c2','#b40f1f','#306311','#df772b','#e184ca',
                           '#de554c','#e8e684','#b2dfe8','#667570','#5ec9b8',
                           '#314e9b','#616e88','#7623f1','#eaeaea'))

# Plot UMAP by Seurat Clusters
pdf(file = "./S2A_dimplot_fibro.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'seurat_clusters', 
        cols = cluster_colors, pt.size = 2, raster = TRUE, reduction = "umap") +
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
# Section 3: Visualization by Tissue Type
################################################################################

tissuecol <- c('#9e0242','#fcae62','#d6404e','#6050a5','#a72b24','#aadda3',
               '#e5f499','#3188bd','#fee18b','#67c2a5','#00b4a9','#049964')
tissuecol <- tissuecol[1:length(unique(scRNA$tissue))]
names(tissuecol) <- c("Primary", "mLN", "CRPC", "Normal", "mCRPC", "LN")

# Split UMAP by Tissue for distribution comparison
pdf(file = "./S2A_dimplot_fibro_tissue.pdf", width = 15, height = 4)
DimPlot(scRNA, label = FALSE, group.by = 'tissue', 
        split.by = 'tissue', cols = tissuecol, pt.size = 4, 
        raster = TRUE, reduction = "umap") +
  NoAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

################################################################################
# Section 4: Cell Cycle Analysis
################################################################################

# Identify G2M and S phase genes available in the dataset
g2m_genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(scRNA))
s_genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(scRNA))

# Perform scoring
scRNA <- CellCycleScoring(scRNA, g2m.features = g2m_genes, s.features = s_genes)

# Prepare data for circular bar plot (Radar-style)
df <- as.data.frame(table(scRNA$cell_subtypes, scRNA$Phase))
colnames(df) <- c("Subtype", "Phase", "Count")
df$Subtype <- factor(df$Subtype, levels = unique(df$Subtype))

pdf(file = "./S2C_cellcycle.pdf", width = 8, height = 6)
ggplot(df, aes(x = Subtype, y = Count, fill = Phase)) +
  geom_bar(stat = "identity") +
  coord_polar(start = 0) +
  theme_minimal() +
  scale_fill_manual(values = c("G1" = '#fffeb3', "G2M" = '#bcb9d8', "S" = '#80b2d4')) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12, margin = margin(t = 10, b = 10)),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(color = "black", size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "right") +
  ylab("") + xlab("")
dev.off()

################################################################################
# Section 5: Final CAF Subtype UMAP (Figure 2A)
################################################################################

caf_palette <- c('#eebad2','#d5adcf','#c4d691','#d1d2e9','#e8bf83',
                 '#b0c3e0','#5c2d8b','#187e8c','#64abc9')

names(caf_palette) <- c("NDUFA4L2+CAF", "CCL2+CAF", "TAGLN+CAF", "CXCL8+CAF", "RUNX2+CAF", 
                        "APOD+CAF", "POSTN+CAF", "MMP2+CAF", "MCAM+CAF")

pdf(file = "./Fig2A_dimplotplot_fibrocelltype.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'cell_subtypes', 
        cols = caf_palette, pt.size = 3, raster = TRUE, reduction = "umap") +
  NoAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank())
dev.off()
