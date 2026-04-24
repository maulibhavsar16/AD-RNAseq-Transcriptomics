# =============================================================
# Script: 02_DESeq2_analysis.R
# Purpose: Differential expression analysis using DESeq2
# Input: GSE159699_summary_count.star.txt, metadata.csv
# Output: DESeq2 results CSV files, dds_object.RData
# Author: Mauli Bhavsar
# Date: 03/04/2026
# =============================================================

# Install Bioconductor manager first
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



# Install supporting packages from CRAN
install.packages(c(
  "tidyverse",     # Data manipulation and plotting
  "ggplot2",       # Publication quality plots
  "pheatmap",      # Heatmaps
  "RColorBrewer"   # Color palettes
))

# Install visualization packages
BiocManager::install(c(
  "EnhancedVolcano",  # Volcano plots
  "apeglm"            # Log fold change shrinkage
))

install.packages("rlang")


library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

cat("All packages loaded successfully!\n")

# ============================================
# Step 2: Load count matrix and metadata
# ============================================

# Load count matrix
counts <- read.table(
  "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\data\\raw\\GSE159699_summary_count.star.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# Load metadata
metadata <- read.csv(
  "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\data\\metadata.csv"
)

# Basic exploration
cat("Count matrix dimensions:\n")
dim(counts)

cat("\nFirst 5 genes, first 5 samples:\n")
counts[1:5, 1:5]

cat("\nMetadata:\n")
print(metadata)

cat("\nSample groups:\n")
table(metadata$diagnosis)

# Fix count matrix column names
# Remove leading X and replace dots with hyphens
colnames(counts) <- gsub("^X", "", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))

# Verify fix
cat("Count matrix column names:\n")
print(colnames(counts))

cat("\nMetadata sample names:\n")
print(metadata$sample)

# Check if they match perfectly
cat("\nDo names match?\n")
all(colnames(counts) == metadata$sample)

# ============================================
# Step 3: Create DESeq2 object
# ============================================

# Convert diagnosis and region to factors
# Factors tell DESeq2 these are categories not numbers
metadata$diagnosis <- factor(metadata$diagnosis, 
                             levels = c("Young", "Old", "AD"))
metadata$region <- factor(metadata$region)

# Set rownames of metadata to sample names
# DESeq2 requires this for matching
rownames(metadata) <- metadata$sample

# Convert counts to matrix format
counts_matrix <- as.matrix(counts)

# Verify counts are integers (DESeq2 requirement)
cat("Are counts integers?\n")
is.integer(counts_matrix)

# Convert to integer if needed
storage.mode(counts_matrix) <- "integer"

# Create DESeq2 object
# design = ~ region + diagnosis means:
# account for brain region differences, then test for diagnosis effect
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = metadata,
  design    = ~ region + diagnosis
)

cat("DESeq2 object created!\n")
print(dds)

# ============================================
# Step 4: Pre-filtering
# ============================================

cat("Genes before filtering:", nrow(dds), "\n")

# Keep only genes with at least 10 counts total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")
cat("Genes removed:", sum(!keep), "\n")

# ============================================
# Step 5: Run DESeq2
# ============================================

cat("Running DESeq2 - this may take a few minutes...\n")

dds <- DESeq(dds)

cat("DESeq2 complete!\n")



# ============================================
# Step 6: Extract DESeq2 results
# ============================================

# Comparison 1: AD vs Young (most dramatic difference)
res_AD_vs_Young <- results(dds,
                           contrast = c("diagnosis", "AD", "Young"),
                           alpha = 0.05)

# Comparison 2: AD vs Old (AD-specific changes)
res_AD_vs_Old <- results(dds,
                         contrast = c("diagnosis", "AD", "Old"),
                         alpha = 0.05)

# Comparison 3: Old vs Young (normal aging)
res_Old_vs_Young <- results(dds,
                            contrast = c("diagnosis", "Old", "Young"),
                            alpha = 0.05)

# Summary of each comparison
cat("=== AD vs Young ===\n")
summary(res_AD_vs_Young)

cat("\n=== AD vs Old ===\n")
summary(res_AD_vs_Old)

cat("\n=== Old vs Young ===\n")
summary(res_Old_vs_Young)


# ============================================
# Save DESeq2 results to CSV files
# ============================================

results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"

# Convert results to dataframes and add gene names
res_AD_vs_Young_df <- as.data.frame(res_AD_vs_Young) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

res_AD_vs_Old_df <- as.data.frame(res_AD_vs_Old) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

res_Old_vs_Young_df <- as.data.frame(res_Old_vs_Young) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Save to CSV
write.csv(res_AD_vs_Young_df,
          paste0(results_dir, "AD_vs_Young_results.csv"),
          row.names = FALSE)

write.csv(res_AD_vs_Old_df,
          paste0(results_dir, "AD_vs_Old_results.csv"),
          row.names = FALSE)

write.csv(res_Old_vs_Young_df,
          paste0(results_dir, "Old_vs_Young_results.csv"),
          row.names = FALSE)

# Save DESeq2 object for later use
save(dds, res_AD_vs_Young, res_AD_vs_Old, res_Old_vs_Young,
     file = paste0(results_dir, "dds_object.RData"))

cat("All results saved successfully!\n")