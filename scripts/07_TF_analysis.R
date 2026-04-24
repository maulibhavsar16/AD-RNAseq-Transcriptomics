# =============================================================
# Script: 07_TF_analysis.R
# Purpose: Transcription Factor Activity Inference
# Input: DESeq2 VST counts, DoRothEA regulons
# Output: TF activity scores and plots
# Author: Mauli
# Date: 18/04/2026
# =============================================================

# Install required packages
BiocManager::install("decoupleR")
BiocManager::install("OmnipathR")

# Load libraries
library(decoupleR)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

cat("Packages loaded!\n")


# =============================================================
# Step 2: Load data and get DoRothEA regulons
# =============================================================

# Set paths
results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"
tf_dir      <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\TF_analysis\\"

# Load VST normalized counts
load(paste0(results_dir, "vsd_object.RData"))
load(paste0(results_dir, "dds_object.RData"))

# Load metadata
metadata <- read.csv(
  "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\data\\metadata.csv"
)

cat("Data loaded!\n")

# Get DoRothEA regulons
# These are TF-target gene relationships
# Confidence levels A-E (A = highest confidence)
# We use A and B only — most reliable interactions
net <- get_dorothea(
  organism   = "human",
  levels     = c("A", "B")
)

cat("DoRothEA regulons loaded!\n")
cat("Number of TF-target interactions:", nrow(net), "\n")
cat("Number of unique TFs:", length(unique(net$source)), "\n")

# Preview
print(head(net))

BiocManager::install("OmnipathR")


library(decoupleR)

net <- get_dorothea(
  organism = "human",
  levels   = c("A", "B")
)

cat("DoRothEA regulons loaded!\n")
cat("Number of TF-target interactions:", nrow(net), "\n")
cat("Number of unique TFs:", length(unique(net$source)), "\n")

print(head(net))


# =============================================================
# Step 3: Run TF activity inference
# =============================================================

# Get VST matrix — genes as rows, samples as columns
# This is the format decoupleR needs
vst_matrix <- assay(vsd)

cat("Expression matrix dimensions:\n")
cat("Genes:", nrow(vst_matrix), "\n")
cat("Samples:", ncol(vst_matrix), "\n")

# Run weighted mean method
# This is the recommended method in decoupleR
# It calculates TF activity for each sample individually
tf_activities <- run_wmean(
  mat  = vst_matrix,
  net  = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  times   = 100,
  minsize = 5
)

cat("\nTF activity inference complete!\n")
cat("Dimensions:", dim(tf_activities), "\n")
print(head(tf_activities))


# =============================================================
# Step 4: Extract significant TF activities
# =============================================================

# Filter for significant statistic type
tf_scores <- tf_activities %>%
  filter(statistic == "corr_wmean")

cat("Total TF-sample combinations:", nrow(tf_scores), "\n")

# Convert to wide matrix
# Rows = TFs, Columns = samples
tf_matrix <- tf_scores %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(
    names_from  = condition,
    values_from = score
  ) %>%
  column_to_rownames("source") %>%
  as.matrix()

cat("TF activity matrix dimensions:\n")
cat("TFs:", nrow(tf_matrix), "\n")
cat("Samples:", ncol(tf_matrix), "\n")

# Add metadata for ordering
# Order samples by diagnosis
sample_order <- metadata %>%
  arrange(diagnosis) %>%
  pull(sample)

# Reorder columns
tf_matrix <- tf_matrix[, sample_order]

cat("\nSample order:\n")
print(colnames(tf_matrix))


# =============================================================
# Step 5: Find differentially active TFs between AD and Old
# =============================================================

# Get sample indices for each group
ad_samples  <- metadata$sample[metadata$diagnosis == "AD"]
old_samples <- metadata$sample[metadata$diagnosis == "Old"]
young_samples <- metadata$sample[metadata$diagnosis == "Young"]

# Calculate mean activity per group for each TF
mean_AD    <- rowMeans(tf_matrix[, ad_samples])
mean_Old   <- rowMeans(tf_matrix[, old_samples])
mean_Young <- rowMeans(tf_matrix[, young_samples])

# Calculate difference — AD vs Old
diff_AD_vs_Old   <- mean_AD - mean_Old
diff_AD_vs_Young <- mean_AD - mean_Young

# Run t-test for each TF — AD vs Old
pvals_AD_vs_Old <- apply(tf_matrix, 1, function(x){
  t.test(x[ad_samples], x[old_samples])$p.value
})

# Run t-test — AD vs Young
pvals_AD_vs_Young <- apply(tf_matrix, 1, function(x){
  t.test(x[ad_samples], x[young_samples])$p.value
})

# Combine into results dataframe
tf_results <- data.frame(
  TF               = rownames(tf_matrix),
  mean_AD          = mean_AD,
  mean_Old         = mean_Old,
  mean_Young       = mean_Young,
  diff_AD_vs_Old   = diff_AD_vs_Old,
  diff_AD_vs_Young = diff_AD_vs_Young,
  pval_AD_vs_Old   = pvals_AD_vs_Old,
  pval_AD_vs_Young = pvals_AD_vs_Young
) %>%
  mutate(
    padj_AD_vs_Old   = p.adjust(pval_AD_vs_Old, method = "BH"),
    padj_AD_vs_Young = p.adjust(pval_AD_vs_Young, method = "BH")
  ) %>%
  arrange(pval_AD_vs_Old)

cat("Top 20 differentially active TFs — AD vs Old:\n")
print(head(tf_results[,
                      c("TF", "mean_AD", "mean_Old",
                        "diff_AD_vs_Old", "pval_AD_vs_Old", "padj_AD_vs_Old")], 20))

# Save results
write.csv(tf_results,
          paste0(tf_dir, "TF_activity_results.csv"),
          row.names = FALSE)

cat("\nTF results saved!\n")


# Reload and view top 20
tf_results <- read.csv(paste0(tf_dir, "TF_activity_results.csv"))

cat("Top 20 differentially active TFs — AD vs Old:\n")
print(head(tf_results[,
                      c("TF", "mean_AD", "mean_Old",
                        "diff_AD_vs_Old", "pval_AD_vs_Old", "padj_AD_vs_Old")], 20))

# =============================================================
# Step 6: TF Activity Heatmap
# =============================================================

# Select top 25 most variable TFs across samples
tf_variance <- apply(tf_matrix, 1, var)
top25_tfs <- names(sort(tf_variance, decreasing = TRUE)[1:25])

# Subset matrix
tf_heatmap_data <- tf_matrix[top25_tfs, ]

# Scale by row for visualization
tf_scaled <- t(scale(t(tf_heatmap_data)))

# Annotation for columns
annotation_col <- data.frame(
  Diagnosis = metadata$diagnosis,
  row.names = metadata$sample
)[sample_order, , drop = FALSE]

annotation_colors <- list(
  Diagnosis = c(
    "AD"    = "#E74C3C",
    "Old"   = "#3498DB",
    "Young" = "#2ECC71"
  )
)

# Plot heatmap
pheatmap(
  tf_scaled,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  color             = colorRampPalette(
    rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  fontsize_row      = 9,
  main              = "TF Activity Scores — Top 25 Variable TFs",
  filename          = paste0(tf_dir, "TF_activity_heatmap.png"),
  width             = 10,
  height            = 10
)

cat("TF activity heatmap saved!\n")

# Dotplot of top differentially active TFs
top_tfs <- tf_results %>%
  filter(pval_AD_vs_Old < 0.05) %>%
  arrange(diff_AD_vs_Old)

ggplot(top_tfs,
       aes(x = diff_AD_vs_Old,
           y = reorder(TF, diff_AD_vs_Old),
           color = diff_AD_vs_Old,
           size  = -log10(pval_AD_vs_Old))) +
  geom_point() +
  scale_color_gradient2(
    low      = "#3498DB",
    mid      = "white",
    high     = "#E74C3C",
    midpoint = 0
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 12) +
  labs(
    title  = "Differentially Active TFs — AD vs Old",
    x      = "Activity Difference (AD - Old)",
    y      = "",
    color  = "Activity\nDifference",
    size   = "-log10(pval)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = paste0(tf_dir, "TF_activity_dotplot.png"),
  width = 10, height = 8, dpi = 300
)

cat("TF dotplot saved!\n")

# Save tf_matrix for future use
save(tf_matrix, tf_activities, tf_scores,
     file = paste0(tf_dir, "TF_objects.RData"))

cat("TF objects saved!\n")
