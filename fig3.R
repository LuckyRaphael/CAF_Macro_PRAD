################################################################################
# Section 1: CytoTRACE Analysis for Differentiation Potential
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig3/')

library(Seurat)
library(CytoTRACE)
library(data.table)
library(tibble)
library(ggplot2)

# Load annotated fibroblast data
scRNA <- qs::qread("../results/Fib_anno.qs")
scRNA$cell_subtypes <- scRNA$celltype

# Extract and filter count matrix
exp1 <- as.matrix(scRNA@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0, 1, sum) >= 5, ] # Filter genes expressed in < 5 cells

# Run CytoTRACE
results <- CytoTRACE(exp1, ncores = 20)

# Prepare metadata for plotting
phenot <- as.character(scRNA$cell_subtypes)
names(phenot) <- rownames(scRNA@meta.data)
emb <- scRNA@reductions[["umap"]]@cell.embeddings

# Output directory and plotting
if (!dir.exists('./CytoTRACE1/')) dir.create('./CytoTRACE1/')

plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './CytoTRACE1/')
plotCytoGenes(results, numOfGenes = 30, outputDir = './CytoTRACE1/')

################################################################################
# Section 2: Monocle2 Pseudotime Trajectory Construction
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig3/')

library(Seurat)
library(monocle)
library(tidydr)
library(viridis)

scRNA <- qs::qread("../results/Fib_anno.qs")
scRNA$cell_subtypes <- scRNA$celltype
top50 <- read.csv('../results/CAF_top50markers.csv')

# Define custom CAF subtype colors
CAFcol <- unique(c('#eebad2','#d5adcf','#c4d691','#d1d2e9','#e8bf83',
                   '#b0c3e0','#5c2d8b','#187e8c','#64abc9','#ffca56','#de3935'))
CAFcol <- CAFcol[1:length(unique(scRNA$cell_subtypes))]
names(CAFcol) <- c("NDUFA4L2+CAF", "CCL2+CAF", "TAGLN+CAF", "CXCL8+CAF", "RUNX2+CAF", 
                   "APOD+CAF", "POSTN+CAF", "MMP2+CAF", "MCAM+CAF")

# 1. Convert to CellDataSet and Preprocess
dat <- Seurat::as.CellDataSet(scRNA)
dat <- estimateSizeFactors(dat)
dat <- detectGenes(dat, min_expr = 10)

# Filter genes by expression frequency
fData(dat)$use_for_ordering <- fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(dat), num_cells_expressed >= 10))

# Preliminary dimension reduction
dat <- reduceDimension(dat, max_components = 2, norm_method = 'log', 
                       num_dim = 20, reduction_method = 'tSNE', 
                       residualModelFormulaStr = "~orig.ident", verbose = T)
dat <- clusterCells(dat, verbose = F)

# 2. Ordering and Trajectory Construction
# Use specific top marker genes for trajectory ordering
dat <- setOrderingFilter(dat, ordering_genes = unique(top50$gene))
dat <- reduceDimension(dat, method = 'DDRTree')
dat <- orderCells(dat)

# Define root state (Starting point of the trajectory)
GM_state <- function(cds, starting_point, cluster){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds@phenoData@data[,cluster])[,starting_point]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else { return (1) }
}

root_start <- GM_state(cds = dat, starting_point = "POSTN+CAF", cluster = "celltype")
dat <- monocle::orderCells(dat, root_state = root_start)

# Save result
saveRDS(dat, "../results/CAF_monocle_top50_2.RDS")

################################################################################
# Section 3: Figure Generation
################################################################################

# 1. Cell Type Wrapped Trajectory
p1 <- plot_cell_trajectory(dat, color_by = "celltype", cell_link_size = 1, 
                           theta = 0.8, show_backbone = FALSE, cell_size = 0.75) +
      facet_wrap(~celltype, nrow = 3) +
      theme_void() + theme(legend.position = "none") +
      scale_color_manual(values = CAFcol)
ggsave('./S3A_CAF_CellType_wrap_top50.pdf', p1, width = 9, height = 9)

# 2. Combined Cell Type Trajectory
p2 <- plot_cell_trajectory(dat, color_by = "celltype", cell_link_size = 1, 
                           cell_size = 1, show_tree = T, show_backbone = FALSE, theta = 0.8) +
      theme_dr(arrow = grid::arrow(length = unit(0, "inches"))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = CAFcol)
ggsave('./Fig3A_CAF_CellType_top50.pdf', p2, width = 8, height = 6)

# 3. Tissue-wise Wrapped Trajectory
p3 <- plot_cell_trajectory(dat, color_by = "celltype", cell_link_size = 1, 
                           theta = 0.8, show_backbone = FALSE, cell_size = 1) +
      facet_wrap(~tissue, nrow = 1) +
      theme_void() + theme(legend.position = "none") +
      scale_color_manual(values = CAFcol)
ggsave('./S3A2_CAF_tissue_wrap_top50.pdf', p3, width = 12, height = 5)

# 4. Pseudotime Trajectory
p4 <- plot_cell_trajectory(dat, color_by = "Pseudotime", cell_size = 1, 
                           show_backbone = FALSE, cell_link_size = 1) +
      scale_colour_gradientn(colours = viridis(100), name = "Pseudotime", 
                             guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      theme_dr(arrow = grid::arrow(length = unit(0, "inches"))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave('./Fig3B_CAF_Pseudotime_top50.pdf', p4, width = 8, height = 6)

# 5. State Visualization and Pie Charts
p5 <- plot_cell_trajectory(dat, color_by = "State", cell_size = 1, 
                           cell_link_size = 1, show_tree = T, show_backbone = FALSE) +
      theme_dr(arrow = grid::arrow(length = unit(0, "inches"))) +
      scale_color_manual(values = jjAnno::useMyCol('calm', n = 7))
ggsave('./Fig3C_CAF_State_top50.pdf', p5, width = 8, height = 6)

# Pie chart distribution across States
tmp <- sapply(1:length(unique(dat$State)), function(x){
  sapply(levels(factor(dat$celltype)), function(y){
    length(which(dat$celltype == y & dat$State == x))
  })
})
colnames(tmp) <- paste0('State', 1:length(unique(dat$State)))

plot_lists <- lapply(colnames(tmp), function(x){
  df <- data.frame(value = as.numeric(as.vector(tmp[,x])),
                   group = factor(rownames(tmp), levels = levels(factor(dat$celltype))))
  ggplot(df, aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = CAFcol) + labs(title = x) +
    coord_polar("y") + theme_minimal() + 
    theme(axis.text = element_blank(), axis.title = element_blank(),
          panel.grid = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 14, face = "bold")) + Seurat::NoLegend()
})

pdf('./Fig3C_Monocle2_State_CellType_pie.pdf', width = 12, height = 4)
cowplot::plot_grid(plotlist = plot_lists, nrow = 1)
dev.off()


