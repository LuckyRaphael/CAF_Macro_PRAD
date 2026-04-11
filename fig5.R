################################################################################
# Section 1: Environment Setup and Library Loading
################################################################################
rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig5/')

library(data.table)
library(GSVA)
library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)
library(patchwork)
library(stringr)

################################################################################
# Section 2: TCGA Bulk Expression Data Preprocessing
################################################################################

# Load TCGA PRAD TPM expression matrix
exp <- read.csv('../results/TCGA_PRAD_TPM.csv', row.names = 1, check.names = FALSE)

# Identify Tumor and Normal samples (TCGA barcode 14-15: <10 is Tumor)
sample_types <- ifelse(as.numeric(str_sub(colnames(exp), 14, 15)) < 10, 'tumor', 'normal')
exp <- exp[, sample_types == "tumor"]

################################################################################
# Section 3: Marker Gene Set Preparation (GSVA/ssGSEA)
################################################################################

# Load Myeloid markers (Top 30 by log2FC)
genesets_mye <- fread("../results/Mye_allmarkers.csv", data.table = FALSE) %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) %>%
  select(gene, cluster)

# Load Fibroblast (CAF) markers (Top 100 by log2FC)
genesets_fib <- fread("../results/CAF_allmarkers.csv", data.table = FALSE) %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  select(gene, cluster)

# Combine and refine gene sets
genesets_all <- rbind(genesets_fib, genesets_mye)
genesets_list <- split(genesets_all$gene, genesets_all$cluster)

# Select target subtypes for survival analysis
select_celltype <- c("POSTN+CAF", "Macro_APOE")
genesets_filtered <- genesets_list[names(genesets_list) %in% select_celltype]

# Calculate Enrichment Scores (ssGSEA)
gsva_res <- gsva(expr = as.matrix(exp), 
                 gset.idx.list = genesets_filtered, 
                 method = "ssgsea", 
                 kcdf = "Gaussian", 
                 parallel.sz = 15)

gsva_res <- as.data.frame(t(gsva_res))

################################################################################
# Section 4: Survival Data Integration
################################################################################

# Load pre-processed TCGA survival data
load('../Fig2/otherbulk/TCGA_PRAD.RData') # Contains TCGA_PRAD object
tcga_suv <- TCGA_PRAD[, c(1:2)] # Assuming columns 1:2 are OS and OS.time
rownames(tcga_suv) <- paste0(rownames(tcga_suv), "A") 

# Clean survival data
tcga_suv <- tcga_suv %>%
  mutate(OS.time = as.numeric(OS.time)) %>%
  drop_na(OS.time) %>%
  filter(OS.time > 0)

# Synchronize samples between expression and survival data
common_samples <- intersect(rownames(tcga_suv), rownames(gsva_res))
rt <- cbind(tcga_suv[common_samples, ], gsva_res[common_samples, ])

# Standardize column name for plotting
colnames(rt)[colnames(rt) == "Macro_APOE"] <- "APOE+Macro"
select_celltype[select_celltype == "Macro_APOE"] <- "APOE+Macro"

################################################################################
# Section 5: Stratification and Grouping
################################################################################

# Determine Optimal Cutpoints for High/Low grouping
variables_to_split <- select_celltype
res_cut <- surv_cutpoint(rt, time = "OS.time", event = "OS", variables = variables_to_split)

# Apply categorization
for (var in variables_to_split) {
  cutoff <- res_cut[[var]][["estimate"]]
  rt[[var]] <- ifelse(rt[[var]] > cutoff, "High", "Low")
}

# Create Combined Groups (High-High vs Low-Low)
markers <- sapply(strsplit(select_celltype, "\\+"), `[`, 1)

rt <- rt %>%
  mutate(Group = paste0(markers[1], "+", .[[select_celltype[1]]], "_", 
                        markers[2], "+", .[[select_celltype[2]]]))

# Filter for the extreme comparison groups
target_groups <- c(paste0(markers[1], "+High_", markers[2], "+High"),
                   paste0(markers[1], "+Low_", markers[2], "+Low"))

rt_filtered <- rt %>%
  filter(Group %in% target_groups) %>%
  mutate(Group = factor(Group, levels = target_groups))

################################################################################
# Section 6: Visualization (Kaplan-Meier Plots)
################################################################################

source('../resource/ggsurvplotKM.R') # Custom plotting function

# 1. Survival by Subtype 1 (POSTN+CAF)
p1 <- ggsurvplotKM(rt[, c("OS.time", "OS", select_celltype[1])],
                   lables = c('High','Low'), 
                   title = select_celltype[1])

# 2. Survival by Subtype 2 (APOE+Macro)
p2 <- ggsurvplotKM(rt[, c("OS.time", "OS", select_celltype[2])],
                   lables = c('High','Low'), 
                   title = select_celltype[2])

# 3. Combined Group Survival (Double High vs Double Low)
p3 <- ggsurvplotKM(rt_filtered[, c("OS.time", "OS", "Group")],
                   lables = levels(rt_filtered$Group), 
                   title = "Combined: CAF & Myeloid")

# Combine and save
final_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("./Fig6A_tcga_km.pdf", final_plot, width = 18, height = 6, dpi = 300)
