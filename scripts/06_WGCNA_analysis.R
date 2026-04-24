# =============================================================
# Script: 06_WGCNA_analysis.R
# Purpose: Weighted Gene Co-expression Network Analysis
# Input: DESeq2 VST normalized counts, metadata
# Output: WGCNA modules, trait correlations, hub genes
# Author: Your Name
# Date: 2025
# =============================================================

#Step 1- Install and load WGCNA packages

# Install WGCNA
install.packages("WGCNA")

# Load libraries
library(WGCNA)
library(DESeq2)
library(tidyverse)

# Critical WGCNA setting — always include this
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

cat("WGCNA loaded successfully!\n")


# Install missing dependency
BiocManager::install("impute")

# Also install another common WGCNA dependency
BiocManager::install("preprocessCore")

# Now load WGCNA
library(WGCNA)
library(DESeq2)
library(tidyverse)


cat("WGCNA loaded successfully!\n")


# Load saved objects
results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"
wgcna_dir   <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\WGCNA\\"

# Load VST normalized counts and metadata
load(paste0(results_dir, "vsd_object.RData"))
load(paste0(results_dir, "dds_object.RData"))

# Get VST matrix — transposed for WGCNA
# WGCNA needs samples as rows, genes as columns
datExpr <- t(assay(vsd))

cat("Expression matrix dimensions:\n")
cat("Samples:", nrow(datExpr), "\n")
cat("Genes:", ncol(datExpr), "\n")

# Load metadata
metadata <- read.csv(
  "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\data\\metadata.csv"
)

cat("\nMetadata loaded:", nrow(metadata), "samples\n")



# ============================================
# Step 3: Quality check for WGCNA
# ============================================

# Check for genes with too many missing values
# Check quality
gsg <- goodSamplesGenes(datExpr, verbose = 3)
cat("All genes and samples pass quality check:", gsg$allOK, "\n")

if(!gsg$allOK){
  if(sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:",
                     paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if(sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

cat("Samples:", nrow(datExpr), "\n")
cat("Genes:", ncol(datExpr), "\n")

# Plot in RStudio AND save to file
sampleTree <- hclust(dist(datExpr), method = "average")

# This shows in RStudio Plots panel
par(cex = 0.8, mar = c(0, 4, 2, 0))
plot(sampleTree,
     main = "Sample Clustering to Detect Outliers",
     sub  = "", xlab = "")

# Save to WGCNA folder
dev.copy(png, paste0(wgcna_dir, "sample_clustering_tree.png"),
         width = 1200, height = 600)
dev.off()

cat("Plot shown and saved!\n")


# ============================================
# Step 4: Filter to top variable genes
# ============================================

# Calculate variance for each gene
gene_variances <- apply(datExpr, 2, var)

# Keep top 5000 most variable genes
# Standard practice for WGCNA on laptops
n_genes <- 5000
top_genes <- order(gene_variances, decreasing = TRUE)[1:n_genes]
datExpr_filtered <- datExpr[, top_genes]

cat("Genes before filtering:", ncol(datExpr), "\n")
cat("Genes after filtering:", ncol(datExpr_filtered), "\n")
cat("Samples retained:", nrow(datExpr_filtered), "\n")



# ============================================
# Step 5: Soft thresholding power selection
# ============================================

# Test powers from 1 to 20
powers <- c(1:20)

sft <- pickSoftThreshold(
  datExpr_filtered,
  powerVector  = powers,
  verbose      = 5,
  networkType  = "signed"
)

cat("Suggested soft threshold power:", sft$powerEstimate, "\n")

# Plot results in RStudio
par(mfrow = c(1, 2))

# Plot 1 — Scale free topology fit
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R²)",
     type = "n",
     main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.8, col = "blue", lty = 2)

# Plot 2 — Mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, col = "red")

# Save plot
dev.copy(png,
         paste0(wgcna_dir, "soft_threshold_selection.png"),
         width = 1200, height = 600)
dev.off()

cat("Soft threshold plot saved!\n")


# ============================================
# Step 6: Build co-expression network
# ============================================
# Fix — override the cor function conflict
cor <- WGCNA::cor

# Now rebuild the network
net <- blockwiseModules(
  datExpr_filtered,
  power             = 9,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

cat("Network built successfully!\n")
cat("Number of modules found:", length(table(net$colors)), "\n")
cat("\nModule sizes:\n")
print(table(net$colors))


# ============================================
# Step 7: Plot module dendrogram
# ============================================

# Plot in RStudio
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors,
  "Module Colors",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene Dendrogram and Module Colors"
)

# Save
dev.copy(png,
         paste0(wgcna_dir, "module_dendrogram.png"),
         width = 1200, height = 600)
dev.off()

cat("Dendrogram saved!\n")


# ============================================
# Step 8: Module-Trait Correlation
# ============================================

# Create numeric trait matrix
# WGCNA needs numbers not categories
traits <- data.frame(
  row.names  = metadata$sample,
  AD         = as.numeric(metadata$diagnosis == "AD"),
  Old        = as.numeric(metadata$diagnosis == "Old"),
  Young      = as.numeric(metadata$diagnosis == "Young"),
  Region_A   = as.numeric(metadata$region == "A"),
  Region_T   = as.numeric(metadata$region == "T")
)

cat("Trait matrix:\n")
print(head(traits))

# Get module eigengenes
# Eigengene = first principal component of each module
# Represents the overall expression pattern of that module
MEs <- moduleEigengenes(datExpr_filtered, net$colors)$eigengenes
MEs <- orderMEs(MEs)

cat("\nModule eigengenes dimensions:", dim(MEs), "\n")
cat("Modules:", colnames(MEs), "\n")

# Calculate correlations between modules and traits
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_filtered))

cat("\nModule-trait correlations calculated!\n")


# ============================================
# Step 9: Module-Trait Correlation Heatmap
# ============================================

# Create text matrix showing correlation and p-value
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Plot heatmap in RStudio
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix       = moduleTraitCor,
  xLabels      = colnames(traits),
  yLabels      = colnames(MEs),
  ySymbols     = colnames(MEs),
  colorLabels  = FALSE,
  colors       = blueWhiteRed(50),
  textMatrix   = textMatrix,
  setStdMargins = FALSE,
  cex.text     = 0.8,
  zlim         = c(-1, 1),
  main         = "Module-Trait Correlation Heatmap"
)

# Save
dev.copy(png,
         paste0(wgcna_dir, "module_trait_heatmap.png"),
         width = 900, height = 700)
dev.off()

cat("Module-trait heatmap saved!\n")


# ============================================
# Step 10: Hub genes in Yellow module
# ============================================

# Calculate module membership (kME)
# kME = correlation of each gene with module eigengene
# Higher kME = more central to the module = hub gene candidate

kME <- signedKME(datExpr_filtered, MEs)

# Focus on yellow module
yellow_genes <- names(net$colors)[net$colors == "yellow"]

# Get kME for yellow module specifically
yellow_kME <- data.frame(
  gene = yellow_genes,
  kME  = kME$kMEyellow[match(yellow_genes, rownames(kME))]
) %>%
  arrange(desc(kME))

cat("Yellow module — Top 20 hub genes:\n")
print(head(yellow_kME, 20))

# Save hub genes
write.csv(yellow_kME,
          paste0(wgcna_dir, "yellow_module_hub_genes.csv"),
          row.names = FALSE)

# Also get hub genes for all significant modules
brown_genes <- names(net$colors)[net$colors == "brown"]
brown_kME <- data.frame(
  gene = brown_genes,
  kME  = kME$kMEbrown[match(brown_genes, rownames(kME))]
) %>%
  arrange(desc(kME))

cat("\nBrown module — Top 20 hub genes:\n")
print(head(brown_kME, 20))

write.csv(brown_kME,
          paste0(wgcna_dir, "brown_module_hub_genes.csv"),
          row.names = FALSE)

cat("\nHub genes saved!\n")


# ============================================
# Step 11: Visualize module membership
# ============================================

# Plot kME distribution for yellow module
par(mfrow = c(1,2))

# Yellow module kME distribution
hist(yellow_kME$kME,
     main = "Yellow Module\nGene-Module Membership",
     xlab = "Module Membership (kME)",
     col  = "gold",
     breaks = 20)
abline(v = 0.8, col = "red", lty = 2)

# Brown module kME distribution
hist(brown_kME$kME,
     main = "Brown Module\nGene-Module Membership",
     xlab = "Module Membership (kME)",
     col  = "brown",
     breaks = 20)
abline(v = 0.8, col = "red", lty = 2)

# Save
dev.copy(png,
         paste0(wgcna_dir, "module_membership_distributions.png"),
         width = 1000, height = 500)
dev.off()

cat("Module membership plots saved!\n")

# Save complete WGCNA objects for reproducibility
save(net, MEs, moduleTraitCor, moduleTraitPvalue,
     kME, datExpr_filtered, traits,
     file = paste0(wgcna_dir, "WGCNA_objects.RData"))

cat("WGCNA objects saved!\n")