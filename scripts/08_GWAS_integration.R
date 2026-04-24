# =============================================================
# Script: 08_GWAS_integration.R
# Purpose: Integrate transcriptomic findings with AD GWAS data
# Input: DESeq2 results, AD GWAS genes from GWAS Catalog
# Output: Overlap analysis, prioritized gene list
# Author: Mauli
# Date: 22/04/2026
# =============================================================

# Load libraries
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)


# Set paths
results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"
gwas_dir    <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\GWAS\\"

# Load DESeq2 results
res_AD_vs_Young_df <- read.csv(paste0(results_dir, "AD_vs_Young_results.csv"))
res_AD_vs_Old_df   <- read.csv(paste0(results_dir, "AD_vs_Old_results.csv"))

cat("Libraries loaded and data ready!\n")

# =============================================================
# Step 2: Get AD GWAS genes from GWAS Catalog
# =============================================================

# Connect to Ensembl biomaRt
mart <- useMart("ensembl",
                dataset = "hsapiens_gene_ensembl")

cat("Connected to Ensembl!\n")

# Known AD GWAS genes from major studies
# Source: GWAS Catalog + Lambert et al. 2013 + Bellenguez et al. 2022
# These are the most replicated AD risk genes
ad_gwas_genes <- c(
  # Most established AD risk genes
  "APOE", "BIN1", "CLU", "PICALM", "CR1",
  "MS4A6A", "MS4A4E", "CD33", "EPHA1", "CD2AP",
  # More recent GWAS hits
  "ABCA7", "SORL1", "FERMT2", "SLC24A4", "RIN3",
  "INPP5D", "MEF2C", "NME8", "ZCWPW1", "CELF1",
  "NYAP1", "PTK2B", "CASS4", "TRIP4", "ECHDC3",
  "ADAMTS4", "ADAMTS1", "ACE", "BICKEL", "PLCG2",
  "ABI3", "TREM2", "LILRB2", "AGRN", "APP",
  "PSEN1", "PSEN2", "MAPT", "GRN", "TOMM40",
  "PVRL2", "NECTIN2", "APOC1", "APOC4", "BCAM",
  "BIN1", "BRIDGE1", "CAND1", "KAT8", "RABEP1"
)

# Remove duplicates
ad_gwas_genes <- unique(ad_gwas_genes)

cat("Total AD GWAS genes:", length(ad_gwas_genes), "\n")
cat("\nAD GWAS gene list:\n")
print(ad_gwas_genes)


# =============================================================
# Reload all required objects
# =============================================================

library(tidyverse)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)

# Set paths
results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"
gwas_dir    <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\GWAS\\"

# Load DESeq2 results
res_AD_vs_Young_df <- read.csv(paste0(results_dir, "AD_vs_Young_results.csv"))
res_AD_vs_Old_df   <- read.csv(paste0(results_dir, "AD_vs_Old_results.csv"))

# Reload GWAS genes
ad_gwas_genes <- c(
  "APOE", "BIN1", "CLU", "PICALM", "CR1",
  "MS4A6A", "MS4A4E", "CD33", "EPHA1", "CD2AP",
  "ABCA7", "SORL1", "FERMT2", "SLC24A4", "RIN3",
  "INPP5D", "MEF2C", "NME8", "ZCWPW1", "CELF1",
  "NYAP1", "PTK2B", "CASS4", "TRIP4", "ECHDC3",
  "ADAMTS4", "ADAMTS1", "ACE", "BICKEL", "PLCG2",
  "ABI3", "TREM2", "LILRB2", "AGRN", "APP",
  "PSEN1", "PSEN2", "MAPT", "GRN", "TOMM40",
  "PVRL2", "NECTIN2", "APOC1", "APOC4", "BCAM",
  "BRIDGE1", "CAND1", "KAT8", "RABEP1"
)
ad_gwas_genes <- unique(ad_gwas_genes)

cat("Everything reloaded!\n")
cat("AD vs Young results:", nrow(res_AD_vs_Young_df), "genes\n")
cat("AD vs Old results:", nrow(res_AD_vs_Old_df), "genes\n")
cat("GWAS genes:", length(ad_gwas_genes), "\n")


# =============================================================
# Step 3: Find overlap with our DEGs
# =============================================================

# Get our significant DEGs
sig_AD_vs_Young <- res_AD_vs_Young_df %>%
  filter(!is.na(padj), padj < 0.05,
         abs(log2FoldChange) > 0.5) %>%
  pull(gene)

sig_AD_vs_Old <- res_AD_vs_Old_df %>%
  filter(!is.na(padj), padj < 0.05,
         abs(log2FoldChange) > 0.5) %>%
  pull(gene)

cat("DEGs AD vs Young:", length(sig_AD_vs_Young), "\n")
cat("DEGs AD vs Old:", length(sig_AD_vs_Old), "\n")
cat("AD GWAS genes:", length(ad_gwas_genes), "\n")

# Find overlaps
overlap_Young <- intersect(sig_AD_vs_Young, ad_gwas_genes)
overlap_Old   <- intersect(sig_AD_vs_Old, ad_gwas_genes)
overlap_both  <- intersect(overlap_Young, overlap_Old)

cat("\n=== OVERLAP RESULTS ===\n")
cat("GWAS genes in AD vs Young DEGs:", length(overlap_Young), "\n")
cat("GWAS genes in AD vs Old DEGs:", length(overlap_Old), "\n")
cat("GWAS genes in BOTH comparisons:", length(overlap_both), "\n")

cat("\nGWAS genes found in AD vs Young:\n")
print(overlap_Young)

cat("\nGWAS genes found in AD vs Old:\n")
print(overlap_Old)

cat("\nGWAS genes found in BOTH comparisons:\n")
print(overlap_both)

# Get their expression details
gwas_deg_details <- res_AD_vs_Young_df %>%
  filter(gene %in% overlap_Young) %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  arrange(padj)

cat("\nExpression details of GWAS-DEG overlapping genes:\n")
print(gwas_deg_details)

# Save overlap results
write.csv(gwas_deg_details,
          paste0(gwas_dir, "GWAS_DEG_overlap.csv"),
          row.names = FALSE)

cat("\nOverlap results saved!\n")


# =============================================================
# Step 4: Test statistical significance of overlap
# =============================================================

# Hypergeometric test
# Is finding 8 GWAS genes in our DEG list more than expected by chance?

# Parameters
total_genes    <- nrow(res_AD_vs_Young_df)  # total genes tested
total_gwas     <- length(ad_gwas_genes)      # total GWAS genes
total_degs     <- length(sig_AD_vs_Young)    # our DEGs
overlap_count  <- length(overlap_Young)      # overlap found

# Hypergeometric test
pval_overlap <- phyper(
  overlap_count - 1,  # observed overlaps minus 1
  total_gwas,         # GWAS genes in universe
  total_genes - total_gwas,  # non-GWAS genes
  total_degs,         # our DEG list size
  lower.tail = FALSE  # probability of >= overlap_count
)

# Expected overlap by chance
expected_overlap <- (total_gwas / total_genes) * total_degs

cat("=== Overlap Statistical Test ===\n")
cat("Total genes tested:", total_genes, "\n")
cat("Total GWAS genes:", total_gwas, "\n")
cat("Our DEGs:", total_degs, "\n")
cat("Observed overlap:", overlap_count, "\n")
cat("Expected overlap by chance:", round(expected_overlap, 2), "\n")
cat("Enrichment fold:", round(overlap_count/expected_overlap, 2), "\n")
cat("Hypergeometric p-value:", pval_overlap, "\n")


# =============================================================
# Step 5: Visualize GWAS-DEG overlap
# =============================================================

# Bar plot of overlapping genes with expression values
gwas_plot_data <- gwas_deg_details %>%
  mutate(
    direction = ifelse(log2FoldChange > 0,
                       "Upregulated in AD",
                       "Downregulated in AD"),
    gene = reorder(gene, log2FoldChange)
  )

ggplot(gwas_plot_data,
       aes(x = log2FoldChange,
           y = gene,
           fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(
    "Upregulated in AD"   = "#E74C3C",
    "Downregulated in AD" = "#3498DB"
  )) +
  theme_classic(base_size = 12) +
  labs(
    title    = "AD GWAS Genes Found in Our DEG List",
    subtitle = "AD vs Young | 8 genes with both transcriptomic and genetic evidence",
    x        = "log2 Fold Change (AD vs Young)",
    y        = "",
    fill     = ""
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  geom_text(aes(label = signif(padj, 2)),
            hjust = ifelse(gwas_plot_data$log2FoldChange > 0,
                           -0.2, 1.2),
            size = 3.5)

ggsave(
  filename = paste0(gwas_dir, "GWAS_DEG_overlap_plot.png"),
  width = 10, height = 6, dpi = 300
)

cat("GWAS overlap plot saved!\n")

# Summary table with biological annotations
gwas_summary <- data.frame(
  Gene = c("TOMM40", "ADAMTS4", "APP", "ABCA7",
           "BIN1", "CASS4", "BCAM", "PSEN1"),
  log2FC = c(-0.56, -2.35, -0.81, 0.81,
             -0.85, 0.95, 0.72, -0.84),
  padj = c(0.0001, 0.0009, 0.0013, 0.0157,
           0.0249, 0.0337, 0.0406, 0.0441),
  GWAS_rank = c("Strong", "Moderate", "Causal",
                "Strong", "Strong", "Moderate",
                "Moderate", "Causal"),
  Biological_role = c(
    "Mitochondrial import — connects to OxPhos failure",
    "ECM metalloprotease — matrix remodeling",
    "Amyloid precursor — central AD gene",
    "Lipid transport — microglial function",
    "Endocytosis — tau pathology",
    "Cytoskeletal signaling",
    "Cell adhesion",
    "Gamma-secretase — amyloid processing"
  )
)

write.csv(gwas_summary,
          paste0(gwas_dir, "GWAS_summary_table.csv"),
          row.names = FALSE)

cat("Summary table saved!\n")
print(gwas_summary)


