# =============================================================
# Script: 03_DESeq2_visualization.R
# Purpose: Visualization of DESeq2 results
# Input: dds_object.RData
# Output: PCA plot, volcano plots, heatmaps
# Author: Mauli Bhavsar
# Date: 03/04/2026
# =============================================================

# Load required packages
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)

# Load saved DESeq2 object
# No need to rerun DESeq2 - just load what we saved!
results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"

load(paste0(results_dir, "dds_object.RData"))

cat("DESeq2 object loaded successfully!\n")
cat("Objects available:", ls(), "\n")


# ============================================
# Plot 1: PCA Plot
# ============================================

# Extract PCA data manually for full control
pca_data <- plotPCA(vsd, 
                    intgroup = c("diagnosis", "region"),
                    returnData = TRUE)

# Check what percentage variance each PC explains
percentVar <- round(100 * attr(pca_data, "percentVar"))

cat("PC1 explains:", percentVar[1], "% variance\n")
cat("PC2 explains:", percentVar[2], "% variance\n")
print(head(pca_data))

# Now plot with full control
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                 color = diagnosis, 
                                 shape = region)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = c(
    "AD"    = "#E74C3C",
    "Old"   = "#3498DB",
    "Young" = "#2ECC71"
  )) +
  scale_shape_manual(values = c("A" = 16, "T" = 17)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of AD, Old and Young Brain Samples",
    subtitle = "Variance Stabilized Counts | n=30 samples",
    color = "Diagnosis",
    shape = "Brain Region"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(pca_plot)

ggsave(
  filename = paste0(results_dir, "PCA_plot.png"),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("PCA plot saved!\n")




# ============================================
# Plot 2: Volcano Plots
# ============================================

# Load results dataframes
res_AD_vs_Young_df <- read.csv(paste0(results_dir, "AD_vs_Young_results.csv"))
res_AD_vs_Old_df   <- read.csv(paste0(results_dir, "AD_vs_Old_results.csv"))
res_Old_vs_Young_df <- read.csv(paste0(results_dir, "Old_vs_Young_results.csv"))

# --- Volcano Plot 1: AD vs Young ---
volcano_AD_vs_Young <- EnhancedVolcano(res_AD_vs_Young_df,
                                       lab     = res_AD_vs_Young_df$gene,
                                       x       = "log2FoldChange",
                                       y       = "padj",
                                       title   = "AD vs Young",
                                       subtitle = paste0("Total DEGs: 3,929 | Up: 1,773 | Down: 2,156"),
                                       pCutoff = 0.05,
                                       FCcutoff = 1.0,
                                       pointSize = 2.0,
                                       labSize = 3.0,
                                       col = c("grey70", "grey70", "steelblue", "#E74C3C"),
                                       colAlpha = 0.6,
                                       legendPosition = "right",
                                       drawConnectors = TRUE,
                                       widthConnectors = 0.5,
                                       max.overlaps = 20
)

print(volcano_AD_vs_Young)

ggsave(
  filename = paste0(results_dir, "AD_vs_Young_volcano.png"),
  plot = volcano_AD_vs_Young,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Volcano plot 1 saved!\n")


# --- Volcano Plot 2: AD vs Old ---
volcano_AD_vs_Old <- EnhancedVolcano(res_AD_vs_Old_df,
                                     lab      = res_AD_vs_Old_df$gene,
                                     x        = "log2FoldChange",
                                     y        = "padj",
                                     title    = "AD vs Old",
                                     subtitle = paste0("Total DEGs: 684 | Up: 321 | Down: 363"),
                                     pCutoff  = 0.05,
                                     FCcutoff = 1.0,
                                     pointSize = 2.0,
                                     labSize  = 3.0,
                                     col      = c("grey70", "grey70", "steelblue", "#E74C3C"),
                                     colAlpha = 0.6,
                                     legendPosition = "right",
                                     drawConnectors = TRUE,
                                     widthConnectors = 0.5,
                                     max.overlaps = 20
)

print(volcano_AD_vs_Old)

ggsave(
  filename = paste0(results_dir, "AD_vs_Old_volcano.png"),
  plot = volcano_AD_vs_Old,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Volcano plot 2 saved!\n")

# --- Volcano Plot 3: Old vs Young ---
volcano_Old_vs_Young <- EnhancedVolcano(res_Old_vs_Young_df,
                                        lab      = res_Old_vs_Young_df$gene,
                                        x        = "log2FoldChange",
                                        y        = "padj",
                                        title    = "Old vs Young",
                                        subtitle = paste0("Total DEGs: 361 | Up: 56 | Down: 305"),
                                        pCutoff  = 0.05,
                                        FCcutoff = 1.0,
                                        pointSize = 2.0,
                                        labSize  = 3.0,
                                        col      = c("grey70", "grey70", "steelblue", "#3498DB"),
                                        colAlpha = 0.6,
                                        legendPosition = "right",
                                        drawConnectors = TRUE,
                                        widthConnectors = 0.5,
                                        max.overlaps = 20
)

print(volcano_Old_vs_Young)

ggsave(
  filename = paste0(results_dir, "Old_vs_Young_volcano.png"),
  plot = volcano_Old_vs_Young,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Volcano plot 3 saved!\n")


# ============================================
# Plot 3: Heatmaps
# ============================================

# Get variance stabilized counts for visualization
vsd_matrix <- assay(vsd)

# Create annotation dataframe for columns
annotation_col <- data.frame(
  Diagnosis = metadata$diagnosis,
  Region    = metadata$region,
  row.names = metadata$sample
)

# Define colors for annotation
annotation_colors <- list(
  Diagnosis = c(
    "AD"    = "#E74C3C",
    "Old"   = "#3498DB",
    "Young" = "#2ECC71"
  ),
  Region = c(
    "A" = "#F39C12",
    "T" = "#8E44AD"
  )
)

# Save vsd object for future use
save(vsd, vsd_matrix,
     file = paste0(results_dir, "vsd_object.RData"))

cat("VSD object saved!\n")

# --- Heatmap 1: Top 50 DEGs from AD vs Young ---

# Get top 50 most significant genes
top50_AD_vs_Young <- res_AD_vs_Young_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

# Extract their expression values
heatmap_data_1 <- vsd_matrix[top50_AD_vs_Young, ]

# Scale by row (gene) so patterns are visible
# This converts to Z-scores — how many SDs above/below mean
heatmap_data_1_scaled <- t(scale(t(heatmap_data_1)))

# Plot heatmap
pheatmap(
  heatmap_data_1_scaled,
  annotation_col  = annotation_col,
  annotation_colors = annotation_colors,
  color           = colorRampPalette(
    rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  fontsize_row    = 8,
  main            = "Top 50 DEGs: AD vs Young",
  filename        = paste0(results_dir, "heatmap_AD_vs_Young.png"),
  width           = 10,
  height          = 12
)


cat("Heatmap 1 saved!\n")

# --- Heatmap 2: Top 50 DEGs from AD vs Old ---
top50_AD_vs_Old <- res_AD_vs_Old_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

heatmap_data_2 <- vsd_matrix[top50_AD_vs_Old, ]
heatmap_data_2_scaled <- t(scale(t(heatmap_data_2)))

pheatmap(
  heatmap_data_2_scaled,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  color             = colorRampPalette(
    rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 8,
  main              = "Top 50 DEGs: AD vs Old",
  filename          = paste0(results_dir, "heatmap_AD_vs_Old.png"),
  width             = 10,
  height            = 12
)

cat("Heatmap 2 saved!\n")
