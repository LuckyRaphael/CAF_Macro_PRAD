################################################################################
# Section 1: Environment Setup and Gene Set Preparation
################################################################################

rm(list = ls())
gc()
setwd('/home/w282308/JIC_PRAD/Fig8/')

library(GSVA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(stringr)
library(patchwork)

# 1. Load and process Macrophage markers (Top 10)
genesets_macro <- fread("../results/Mye_allmarkers.csv", data.table = FALSE) %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  select(gene, cluster)

# 2. Load and process CAF markers (Top 100)
genesets_fib <- fread("../results/CAF_allmarkers.csv", data.table = FALSE) %>% 
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  select(gene, cluster)

# 3. Combine and Format for GSVA
genesets_all <- rbind(genesets_macro, genesets_fib)
genesets_list <- split(genesets_all$gene, genesets_all$cluster)

# Sanitize names (replace special characters with underscores for compatibility)
names(genesets_list) <- gsub("[-+]", "_", names(genesets_list))

# Filter for specific sub-clusters of interest
select_celltype <- c("POSTN_CAF", "Macro_APOE")
genesets_filtered <- genesets_list[names(genesets_list) %in% select_celltype]

################################################################################
# Section 2: Batch Processing Multi-Cohort Correlation
################################################################################

# Identify all RData cohort files
data_dir <- "../Fig2/otherbulk/"
rdata_files <- list.files(path = data_dir, pattern = "\\.RData$", full.names = TRUE)

plot_list <- list()

# Loop through each cohort
for (file in rdata_files) {
  var_name <- sub("\\.RData$", "", basename(file))
  message("Processing Cohort: ", var_name)
  
  # Load the expression matrix (assuming the object name matches the filename)
  load(file)
  exp_data <- get(var_name)
  exp_data <- na.omit(exp_data)
  
  # Perform ssGSEA scoring
  # Note: t(exp_data[, -c(1:2)]) assumes columns 1 and 2 are metadata (like Survival info)
  GSVA_res <- gsva(expr = as.matrix(t(exp_data[, -c(1:2)])),
                   gset.idx.list = genesets_filtered,
                   method = "ssgsea",
                   kcdf = "Gaussian",
                   mx.diff = TRUE,
                   parallel.sz = 15)
  
  df_plot <- as.data.frame(t(GSVA_res))
  
  # Generate Correlation Scatter Plot
  p <- ggplot(df_plot, aes_string(x = select_celltype[1], y = select_celltype[2])) +
    geom_point(color = "black", alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", fill = "pink", size = 1.2) +
    stat_cor(method = "spearman", digits = 3, size = 4.5, label.sep = "\n") +
    ggtitle(var_name) +
    labs(x = paste(select_celltype[1], "Score"), 
         y = paste(select_celltype[2], "Score")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  plot_list[[var_name]] <- p
}

################################################################################
# Section 3: Final Composition and Export
################################################################################

# Arrange all cohort plots into a grid
combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_annotation(title = 'Correlation between POSTN+ CAF and APOE+ Macro across Cohorts',
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

# Save the multi-panel plot
ggsave('./multi_cohort_corplot.pdf', combined_plot, height = 10, width = 12)

# Print to device
print(combined_plot)
