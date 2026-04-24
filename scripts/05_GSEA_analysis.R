# =============================================================
# Script: 05_GSEA_analysis.R
# Purpose: Gene Set Enrichment Analysis (GSEA)
# Input: DESeq2 results CSV files
# Output: GSEA results and plots
# Author: Mauli Bhavsar
# Date: 03/04/2026
# =============================================================

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(ggplot2)

# Set paths
results_dir <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\DESeq2\\"
enrich_dir  <- "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\results\\enrichment\\"

# Load DESeq2 results
res_AD_vs_Young_df  <- read.csv(paste0(results_dir, "AD_vs_Young_results.csv"))
res_AD_vs_Old_df    <- read.csv(paste0(results_dir, "AD_vs_Old_results.csv"))

cat("Results loaded successfully!\n")


#STEP-1 Prepare MSigDB Gene Set
# ============================================
# GSEA — Gene Set Enrichment Analysis
# ============================================
# Check column names first
hallmarks_raw <- msigdbr(species = "Homo sapiens",
                         collection = "H")
cat("Column names:\n")
print(colnames(hallmarks_raw))

# Check available subcollections
all_sets <- msigdbr(species = "Homo sapiens",
                    collection = "C2")
cat("Available subcollections:\n")
print(unique(all_sets$gs_subcollection))


# Get KEGG gene sets using correct subcollection name
kegg_sets <- msigdbr(species = "Homo sapiens",
                     collection = "C2",
                     subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  mutate(ncbi_gene = as.character(ncbi_gene))

cat("Hallmark gene sets:", length(unique(hallmarks$gs_name)), "\n")
cat("KEGG gene sets:", length(unique(kegg_sets$gs_name)), "\n")

# Convert to list format
hallmarks_list <- split(hallmarks$ncbi_gene,
                        hallmarks$gs_name)
kegg_list <- split(kegg_sets$ncbi_gene,
                   kegg_sets$gs_name)


# Prepare ranked gene lists
# AD vs Young
ranked_Young <- res_AD_vs_Young_df %>%
  filter(!is.na(log2FoldChange)) %>%
  filter(!is.na(padj)) %>%
  arrange(desc(log2FoldChange))

ranked_Young_entrez <- bitr(
  ranked_Young$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
) %>%
  left_join(
    ranked_Young %>% dplyr::select(gene, log2FoldChange),
    by = c("SYMBOL" = "gene")
  ) %>%
  arrange(desc(log2FoldChange))

ranks_Young <- ranked_Young_entrez$log2FoldChange
names(ranks_Young) <- ranked_Young_entrez$ENTREZID

# AD vs Old
ranked_Old <- res_AD_vs_Old_df %>%
  filter(!is.na(log2FoldChange)) %>%
  filter(!is.na(padj)) %>%
  arrange(desc(log2FoldChange))

ranked_Old_entrez <- bitr(
  ranked_Old$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
) %>%
  left_join(
    ranked_Old %>% dplyr::select(gene, log2FoldChange),
    by = c("SYMBOL" = "gene")
  ) %>%
  arrange(desc(log2FoldChange))

ranks_Old <- ranked_Old_entrez$log2FoldChange
names(ranks_Old) <- ranked_Old_entrez$ENTREZID

cat("Ranked genes AD vs Young:", length(ranks_Young), "\n")
cat("Ranked genes AD vs Old:", length(ranks_Old), "\n")

# Run GSEA — Hallmarks
set.seed(42)
gsea_hallmarks_Young <- fgsea(
  pathways = hallmarks_list,
  stats    = ranks_Young,
  minSize  = 15,
  maxSize  = 500
)

set.seed(42)
gsea_hallmarks_Old <- fgsea(
  pathways = hallmarks_list,
  stats    = ranks_Old,
  minSize  = 15,
  maxSize  = 500
)

# Run GSEA — KEGG
set.seed(42)
gsea_kegg_Young <- fgsea(
  pathways = kegg_list,
  stats    = ranks_Young,
  minSize  = 15,
  maxSize  = 500
)

set.seed(42)
gsea_kegg_Old <- fgsea(
  pathways = kegg_list,
  stats    = ranks_Old,
  minSize  = 15,
  maxSize  = 500
)

cat("\nGSEA Hallmarks significant — AD vs Young:",
    sum(gsea_hallmarks_Young$padj < 0.05, na.rm = TRUE), "\n")
cat("GSEA Hallmarks significant — AD vs Old:",
    sum(gsea_hallmarks_Old$padj < 0.05, na.rm = TRUE), "\n")
cat("GSEA KEGG significant — AD vs Young:",
    sum(gsea_kegg_Young$padj < 0.05, na.rm = TRUE), "\n")
cat("GSEA KEGG significant — AD vs Old:",
    sum(gsea_kegg_Old$padj < 0.05, na.rm = TRUE), "\n")



# ============================================
# View top GSEA results
# ============================================

# Sort by NES (Normalized Enrichment Score)
# Positive NES = upregulated in AD
# Negative NES = downregulated in AD

cat("=== Top Hallmarks — AD vs Young ===\n")
gsea_hallmarks_Young %>%
  filter(padj < 0.05) %>%
  arrange(NES) %>%
  dplyr::select(pathway, NES, pval, padj) %>%
  print()

cat("\n=== Top Hallmarks — AD vs Old ===\n")
gsea_hallmarks_Old %>%
  filter(padj < 0.05) %>%
  arrange(NES) %>%
  dplyr::select(pathway, NES, pval, padj) %>%
  print()

cat("\n=== Top KEGG — AD vs Young (top 15) ===\n")
gsea_kegg_Young %>%
  filter(padj < 0.05) %>%
  arrange(NES) %>%
  dplyr::select(pathway, NES, pval, padj) %>%
  head(15) %>%
  print()

cat("\n=== Top KEGG — AD vs Old (top 15) ===\n")
gsea_kegg_Old %>%
  filter(padj < 0.05) %>%
  arrange(NES) %>%
  dplyr::select(pathway, NES, pval, padj) %>%
  head(15) %>%
  print()



# ============================================
# GSEA Visualization
# ============================================

# Plot 1 — Hallmarks bubble plot — AD vs Young
hallmarks_Young_sig <- gsea_hallmarks_Young %>%
  filter(padj < 0.05) %>%
  arrange(NES)

hallmarks_Young_sig$pathway <- gsub(
  "HALLMARK_", "",
  hallmarks_Young_sig$pathway
)

ggplot(hallmarks_Young_sig,
       aes(x = NES,
           y = reorder(pathway, NES),
           size = -log10(padj),
           color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(
    low  = "#3498DB",
    mid  = "white",
    high = "#E74C3C",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(3, 10)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 11) +
  labs(
    title    = "GSEA Hallmarks — AD vs Young",
    x        = "Normalized Enrichment Score (NES)",
    y        = "",
    size     = "-log10(padj)",
    color    = "NES"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  filename = paste0(enrich_dir, "GSEA_hallmarks_Young_bubble.png"),
  width = 12, height = 10, dpi = 300
)

# Plot 2 — Hallmarks bubble plot — AD vs Old
hallmarks_Old_sig <- gsea_hallmarks_Old %>%
  filter(padj < 0.05) %>%
  arrange(NES)

hallmarks_Old_sig$pathway <- gsub(
  "HALLMARK_", "",
  hallmarks_Old_sig$pathway
)

ggplot(hallmarks_Old_sig,
       aes(x = NES,
           y = reorder(pathway, NES),
           size = -log10(padj),
           color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(
    low  = "#3498DB",
    mid  = "white",
    high = "#E74C3C",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(3, 10)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 11) +
  labs(
    title = "GSEA Hallmarks — AD vs Old",
    x     = "Normalized Enrichment Score (NES)",
    y     = "",
    size  = "-log10(padj)",
    color = "NES"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  filename = paste0(enrich_dir, "GSEA_hallmarks_Old_bubble.png"),
  width = 12, height = 10, dpi = 300
)

cat("GSEA bubble plots saved!\n")

# Plot 3 — Classic GSEA enrichment plot for top pathway
# Oxidative phosphorylation — our most important finding
plotEnrichment(
  hallmarks_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
  ranks_Young
) +
  labs(
    title = "GSEA Enrichment Plot",
    subtitle = "HALLMARK_OXIDATIVE_PHOSPHORYLATION — AD vs Young"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  filename = paste0(enrich_dir, "GSEA_OxPhos_enrichment_plot.png"),
  width = 10, height = 6, dpi = 300
)

cat("Enrichment plot saved!\n")

# Save all GSEA results
# Function to convert list columns to character
convert_lists <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.list(x)) {
      sapply(x, paste, collapse = ";")
    } else {
      x
    }
  })
  return(df)
}

# Apply conversion + save
write.csv(
  convert_lists(as.data.frame(gsea_hallmarks_Young)),
  paste0(enrich_dir, "GSEA_hallmarks_Young_results.csv"),
  row.names = FALSE
)

write.csv(
  convert_lists(as.data.frame(gsea_hallmarks_Old)),
  paste0(enrich_dir, "GSEA_hallmarks_Old_results.csv"),
  row.names = FALSE
)

write.csv(
  convert_lists(as.data.frame(gsea_kegg_Young)),
  paste0(enrich_dir, "GSEA_KEGG_Young_results.csv"),
  row.names = FALSE
)

write.csv(
  convert_lists(as.data.frame(gsea_kegg_Old)),
  paste0(enrich_dir, "GSEA_KEGG_Old_results.csv"),
  row.names = FALSE
)

cat("All GSEA results saved!\n")