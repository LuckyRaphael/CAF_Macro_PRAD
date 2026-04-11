################################################################################
# Section 1: Environment Setup and Library Loading
################################################################################
rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig6/')

library(Seurat)
library(ggplot2)
library(scRNAtoolVis)
library(RColorBrewer)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(SingleCellExperiment)
library(BayesSpace)
library(GSVA)
library(dplyr)

# Create output directory for spatial plots
main_fig_dir <- "./stplot/"
if (!dir.exists(main_fig_dir)) {
  dir.create(main_fig_dir, recursive = TRUE)
}

################################################################################
# Section 2: Spatial Data Processing Function
################################################################################

# Function to standardize, reduce dimensions, and cluster a single ST sample
process_spatial_data <- function(st, 
                                 assay = "Spatial", 
                                 npcs = 20,
                                 resolution = 0.5,
                                 ndims = 30) {
  # 1. Normalization using SCTransform
  st <- SCTransform(st, assay = assay, verbose = TRUE)
  
  # 2. PCA Dimensional Reduction
  st <- RunPCA(st, assay = "SCT", verbose = FALSE)
  print(ElbowPlot(st, ndims = ndims, reduction = "pca"))
  
  # 3. Clustering
  st <- FindNeighbors(st, reduction = "pca", dims = 1:npcs)
  st <- FindClusters(st, resolution = resolution, verbose = FALSE)
  print(table(st@meta.data$seurat_clusters))
  
  # 4. UMAP Embedding
  st <- RunUMAP(st, reduction = "pca", dims = 1:npcs)
  
  return(st)
}

# List available ST samples
samples <- list.files(path = "../data/sp/ST_all/", full.names = FALSE, recursive = FALSE)

# Define custom color palette
st_colors <- unique(c('#2867a9','#ca3135','#f57d34','#8469ae','#469d47',
                    '#936962','#9EDDCE','#E3CCE2','#F9C1AC','#FCD6A8','#FBF9D4',
                    '#B8D7EA','#F8BBB5','#BCE5DF','#DEDBED',"#93CA87",
                    "#ED614A","#F4D020", brewer.pal(9, "Set1"), 
                    brewer.pal(8, "Set2"), brewer.pal(12, "Set3")))

################################################################################
# Section 3: Gene Set Preparation (Hallmarks & Cell Type Markers)
################################################################################

## 1. Process Hallmark Pathway (Hypoxia)
gmt_raw <- read.gmt('../resource/h.all.v2024.1.Hs.symbols.gmt')
gmt_list <- split(gmt_raw$gene, gsub("^HALLMARK_", "", gmt_raw$term))
hypoxia_gmt <- gmt_list[names(gmt_list) %in% "HYPOXIA"]

## 2. Process Cell Type Markers (Top 100 markers by log2FC)
# Load markers from Main Type, Myeloid, and Fibroblast analysis
genesets_main <- fread("../results/Maintype_allmarkers.csv", data.table = FALSE) %>%
  group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% select(gene, cluster)

genesets_mye <- fread("../results/Mye_allmarkers.csv", data.table = FALSE) %>%
  group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% select(gene, cluster)

genesets_fib <- fread("../results/CAF_allmarkers.csv", data.table = FALSE) %>%
  group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% select(gene, cluster)

# Combine and filter for target cell types
combined_markers <- rbind(genesets_main, genesets_fib, genesets_mye)
marker_list <- split(combined_markers$gene, combined_markers$cluster)

select_celltype <- c("Epithelia", "T", "POSTN+CAF", "Macro_APOE")
marker_list_filtered <- marker_list[names(marker_list) %in% select_celltype]

## 3. Merge All Gene Sets for Scoring
genesets_final <- c(hypoxia_gmt, marker_list_filtered)
names(genesets_final) <- gsub("/", "", names(genesets_final)) # Clean names for GSVA

# Preview final gene sets
print(names(genesets_final))
