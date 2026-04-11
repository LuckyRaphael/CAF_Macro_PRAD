################################################################################
# Section 1: Initial Data Processing & Metadata Standardization
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig1/')

library(Seurat)
library(tidydr)

# Read single-cell data
scRNA <- qs::qread("../results/StepF.All_Cells_Seurat_Object.qs")

# Inspect data structure
dim(scRNA)
colnames(scRNA@meta.data)
table(scRNA$ct.sub)
table(scRNA$group)
table(scRNA$data.sets)

# Standardize metadata column names
scRNA$Sample <- scRNA$sampleID
scRNA$tissue <- scRNA$group
scRNA$Datasets <- scRNA$data.sets
scRNA$celltype <- scRNA$ct.sub

# Save updated object
qs::qsave(scRNA, "../results/StepF.All_Cells_Seurat_Object.qs")

################################################################################
# Section 2: Marker Gene Identification for Cell Type Scoring
################################################################################

# Enable parallel processing
future::plan("multicore", workers = 10)

library(dplyr)
library(data.table)

scRNA@active.ident <- factor(scRNA$celltype)

# Define list of genes to exclude (Ribosomal, Mitochondrial, Non-coding, etc.)
badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)

badgenelt <- unique(grep(pattern = paste(c(badterm), collapse = "|"), 
                         x = rownames(scRNA), perl = T, value = T))

# Identify marker genes for major cell types (downsampled for performance)
markers <- FindAllMarkers(
  subset(x = scRNA, downsample = 500),
  features = rownames(scRNA)[-which(rownames(scRNA) %in% badgenelt)],
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

markers$pct_diff <- markers$pct.1 - markers$pct.2

# Extract Top 50 markers per cluster
top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

# Export marker lists
saveRDS(markers, "../results/Maintype_CellType_diffgenes.RDS") 
fwrite(markers, "../results/Maintype_allmarkers.csv")
fwrite(top50, "../results/Maintype_top50markers.csv")

################################################################################
# Section 3: Subsetting Data for Downstream Lineage Analysis
################################################################################

# After major cell type annotation, save subsets for sub-clustering
cell_types <- c("Epithelia", "Fib", "Endo", "Mye", "B", "T")
file_names <- c("epi", "fibro", "endo", "mye", "B", "T")

for(i in seq_along(cell_types)){
  tmp_subset <- subset(scRNA, subset = celltype == cell_types[i])
  qs::qsave(tmp_subset, paste0("../results/scRNA_", file_names[i], "_raw.qs"))
}

################################################################################
# Section 4: Figure 1 Visualizations (UMAP & Proportions)
################################################################################

library(ggplot2)
library(cowplot)
library(stringr)
library(viridis)
library(scCustomize)
library(scRNAtoolVis)
library(RColorBrewer)

# Global Color Palette
colors <- unique(c('#639dd0','#e1959a','#ffc79a','#7cc1f3','#e1eddc',
                   '#e84a34','#4dbbd6','#3c5487','#f29b80','#8591b4','#808080',
                   '#f09c00','#016170','#b11f13','#016170','#a7d3b4','#e6daa5',
                   '#005f73','#099396','#ca6702','#ed9b00','#a12c2f','#bb3f02',
                   '#b15d2e','#e82036','#237bba','#40ac53',
                   '#f1c742','#7498b5','#e7e7e7','#d66969',
                   '#e5f499','#3188bd','#fee18b','#67c1a6',
                   '#9e0242','#d6404e','#fcae62','#6050a5','#aadda3','#e5f499',
                   '#3188bd','#fee18b','#67c2a5','#00b4a9','#049964','#0968b3',
                   '#04c0ec','#94288e','#f47a1f','#a72b24','#f59f83','#eb1275',
                   '#e8733d','#981118','#000000','#e9c56d','#006a81','#199169',
                   '#d37db5','#ad60a5','#7e2e92','#e775ad','#01afb2','#86d3e4',
                   '#4c92ce','#a0c3e7','#6d65ae','#850663','#ffbde5','#7860f5',
                   '#2000fc','#9dfffd','#01afb2','#46da44','#b8712e','#fe691e',
                   '#c5cdf7','#7092d2','#e69fc4','#b6854d','#8dd1c6','#fffeb3',
                   '#bcb9d8','#80b2d4','#feb662','#b4df65','#d9d9d9','#be7fbf',
                   '#cbecc4','#ffee72','#e6bdcb'))

# 1. Plot UMAP by Dataset
pdf(file = "./dimplotplot_datasets.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'Datasets', 
        cols = colors, pt.size = 1.5, raster = TRUE, reduction = "umap") +
  NoAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

# 2. Plot UMAP by Tissue Type
tissuecol <- c('#9e0242','#fcae62','#d6404e','#6050a5','#a72b24','#aadda3',
               '#e5f499','#3188bd','#fee18b','#67c2a5','#00b4a9','#049964',
               '#0968b3','#04c0ec','#94288e','#f47a1f','#f59f83','#eb1275')
tissuecol <- tissuecol[1:length(unique(scRNA$tissue))]
names(tissuecol) <- c("Primary", "mLN", "CRPC", "Normal", "mCRPC", "LN")

pdf(file = "./dimplotplot_tissue.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'tissue', 
        cols = tissuecol, pt.size = 1.5, raster = TRUE, reduction = "umap") +
  NoAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

# 3. Pie Chart: Cell count distribution per Dataset
mynames <- table(scRNA$Datasets) %>% names()
myratio <- table(scRNA$Datasets) %>% as.numeric()
pielabel <- paste0(mynames, " (", round(myratio/sum(myratio)*100, 2), "%)")

pdf(file = "./S1B_pieplot_datasets_cellnumber.pdf", width = 13, height = 9)
pie(myratio, labels = pielabel, radius = 1, clockwise = F,
    main = "Cell Distribution per Dataset", border = "white", col = colors)
dev.off()

# 4. Pie Chart: Sample count distribution per Dataset
dataset_sample_counts <- scRNA@meta.data %>%
  distinct(Datasets, Sample) %>%
  dplyr::count(Datasets, name = "sample_count")

mynames_s <- dataset_sample_counts$Datasets
myratio_s <- dataset_sample_counts$sample_count
pielabel_s <- paste0(mynames_s, " (", round(myratio_s/sum(myratio_s)*100, 2), "%)")

pdf(file = "./S1C_pieplot_datasets_samplenumber.pdf", width = 13, height = 9)
pie(myratio_s, labels = pielabel_s, radius = 1, clockwise = FALSE,
    main = "Sample Distribution per Dataset", border = "white", col = colors)
dev.off()

# 5. Final Major Cell Type UMAP
Main_colors <- unique(c('#da1921','#5c2d8b','#187e8c','#64abc9','#ffca56','#343031',
                        '#424b00','#376696','#a6cbaa','#273669','#a8bbd3','#407f31',
                        '#e54a32','#00a087','#f29b80','#91d1c1','#4cbbd6','#395388'))

Maincol <- Main_colors[1:length(unique(scRNA$celltype))]
names(Maincol) <- c("Epithelia", "T", "Mye", "Fib", "Endo", "Mast", "B")

pdf(file = "./dimplotplot_celltype.pdf", width = 8, height = 6)
DimPlot(scRNA, label = FALSE, group.by = 'celltype', 
        cols = Maincol, pt.size = 1.5, raster = TRUE, reduction = "umap") +
  NoAxes() +
  labs(x = "UMAP1", y = "UMAP2", title = "Global Cell Types") +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()




