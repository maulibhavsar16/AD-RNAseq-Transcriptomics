# =============================================================
# Script: 04_enrichment_analysis.R
# Purpose: Functional enrichment analysis (ORA and GSEA)
# Input: DESeq2 results CSV files
# Output: Enrichment results and plots
# Author: Mauli Bhavsar
# Date: 03/04/2026
# =============================================================


# Install enrichment packages
BiocManager::install(c(
    "clusterProfiler",   # ORA and GSEA
    "org.Hs.eg.db",      # Human gene annotation database
    "enrichplot",        # Visualization of enrichment results
    "fgsea"              # Fast GSEA implementation
))

# Install MSigDB gene sets
BiocManager::install("msigdbr")

# Install from CRAN
install.packages("ggupset")  # Required for some enrichplot visualizations

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
res_Old_vs_Young_df <- read.csv(paste0(results_dir, "Old_vs_Young_results.csv"))

cat("Results loaded successfully!\n")


# ============================================
# Step 2: Prepare gene lists
# ============================================

# --- For ORA: we need significant DEGs only ---

# AD vs Old significant genes (padj < 0.05, |log2FC| > 0.5)
# We focus on AD vs Old — our most biologically meaningful comparison
sig_AD_vs_Old <- res_AD_vs_Old_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 0.5)

# Split into up and down
up_AD_vs_Old   <- sig_AD_vs_Old %>% filter(log2FoldChange > 0) %>% pull(gene)
down_AD_vs_Old <- sig_AD_vs_Old %>% filter(log2FoldChange < 0) %>% pull(gene)
all_AD_vs_Old  <- sig_AD_vs_Old %>% pull(gene)

cat("AD vs Old significant DEGs:", nrow(sig_AD_vs_Old), "\n")
cat("Upregulated:", length(up_AD_vs_Old), "\n")
cat("Downregulated:", length(down_AD_vs_Old), "\n")

# --- For GSEA: we need ALL genes ranked by log2FC ---
# Remove NA values and create named vector
ranked_genes_AD_vs_Old <- res_AD_vs_Old_df %>%
  filter(!is.na(log2FoldChange)) %>%
  filter(!is.na(padj)) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange)

# Create named numeric vector — required format for GSEA
gene_ranks <- ranked_genes_AD_vs_Old$log2FoldChange
names(gene_ranks) <- ranked_genes_AD_vs_Old$gene

cat("\nTotal genes for GSEA ranking:", length(gene_ranks), "\n")
cat("Top 5 upregulated:\n")
print(head(gene_ranks, 5))
cat("Top 5 downregulated:\n")
print(tail(gene_ranks, 5))


# ============================================
# Step 3: Convert gene symbols to Entrez IDs
# ============================================

# Convert significant DEGs to Entrez IDs
sig_entrez <- bitr(
  all_AD_vs_Old,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

cat("Genes successfully converted:", nrow(sig_entrez), "\n")
cat("Genes lost in conversion:", length(all_AD_vs_Old) - nrow(sig_entrez), "\n")
print(head(sig_entrez))

# Convert full ranked list for GSEA
ranked_entrez <- bitr(
  names(gene_ranks),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# Merge with fold changes
ranked_entrez <- ranked_entrez %>%
  left_join(
    data.frame(SYMBOL = names(gene_ranks),
               log2FC = gene_ranks),
    by = "SYMBOL"
  ) %>%
  arrange(desc(log2FC))

# Create named vector with Entrez IDs
gene_ranks_entrez <- ranked_entrez$log2FC
names(gene_ranks_entrez) <- ranked_entrez$ENTREZID

cat("\nGenes in GSEA ranked list after conversion:", 
    length(gene_ranks_entrez), "\n")

# ============================================
# Step 4: GO Enrichment Analysis (ORA)
# ============================================

#THIS DIDN'T WORK (FROM HERE)
# Create background gene universe
# These are ALL genes we tested in DESeq2 — not the whole genome
universe_entrez <- bitr(
  res_AD_vs_Old_df$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

cat("Background universe size:", nrow(universe_entrez), "\n")

# Re-run GO enrichment with proper universe
go_BP <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = universe_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_MF <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = universe_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("GO Biological Process terms found:", nrow(go_BP), "\n")
cat("GO Molecular Function terms found:", nrow(go_MF), "\n")

cat("\nTop 10 GO Biological Processes:\n")
print(head(as.data.frame(go_BP)[, c("Description", "GeneRatio", "pvalue", "p.adjust")], 10))

# THIS DIDN'T WORK (UPTO HERE) SO WE DIAGNOSED WHAT WENT WRONG (SEE CHECKS BELOW)

# ============================================
# Diagnosis — check what is happening
# ============================================

# Check 1 — how many genes in sig_entrez?
cat("Significant genes with Entrez IDs:", nrow(sig_entrez), "\n")
print(head(sig_entrez))

# Check 2 — are these valid Entrez IDs?
cat("\nFirst 10 Entrez IDs:\n")
print(sig_entrez$ENTREZID[1:10])

# Check 3 — check universe size
cat("\nUniverse size:", nrow(universe_entrez), "\n")

# Check 4 — do sig genes exist in universe?
overlap <- sum(sig_entrez$ENTREZID %in% universe_entrez$ENTREZID)
cat("Sig genes found in universe:", overlap, "\n")

# Check 5 — try without universe and relaxed cutoffs
go_test <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1.0,
  qvalueCutoff  = 1.0,
  readable      = TRUE
)

cat("\nTerms found with relaxed cutoffs:", nrow(go_test), "\n")
if(nrow(go_test) > 0){
  print(head(as.data.frame(go_test)[, c("Description", "pvalue", "p.adjust")], 5))
}


# DIAGNOSIS COMPLETE, ALTERNATE GIVEN BELOW
#Solution — Use AD vs Young DEGs instead
#AD vs Young has 3,929 DEGs — much more power for enrichment analysis:

  
  
  # Convert AD vs Young significant genes
  sig_AD_vs_Young <- res_AD_vs_Young_df %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05) %>%
    filter(abs(log2FoldChange) > 0.5)
  
  sig_entrez_Young <- bitr(
    sig_AD_vs_Young$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  universe_Young <- bitr(
    res_AD_vs_Young_df$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  cat("AD vs Young significant genes:", nrow(sig_AD_vs_Young), "\n")
  cat("Converted to Entrez:", nrow(sig_entrez_Young), "\n")
  
  # Run GO BP with larger gene list
  go_BP <- enrichGO(
    gene          = sig_entrez_Young$ENTREZID,
    universe      = universe_Young$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_MF <- enrichGO(
    gene          = sig_entrez_Young$ENTREZID,
    universe      = universe_Young$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  cat("GO BP terms found:", nrow(go_BP), "\n")
  cat("GO MF terms found:", nrow(go_MF), "\n")
  
  cat("\nTop 10 GO Biological Processes:\n")
  print(head(as.data.frame(go_BP)[,
                                  c("Description", "GeneRatio", "pvalue", "p.adjust")], 10))

  
  
  # ============================================
  # Step 5: Visualize GO enrichment results
  # ============================================
  
  # --- Plot 1: Dotplot ---
  dotplot_BP <- dotplot(go_BP,
                        showCategory = 15,
                        title = "GO Biological Process — AD vs Young",
                        font.size = 10
  ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(dotplot_BP)
  
  ggsave(
    filename = paste0(enrich_dir, "GO_BP_dotplot.png"),
    plot = dotplot_BP,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  cat("GO BP dotplot saved!\n")
  
  # --- Plot 2: Barplot ---
  barplot_BP <- barplot(go_BP,
                        showCategory = 15,
                        title = "GO Biological Process — AD vs Young",
                        font.size = 10
  ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(barplot_BP)
  
  ggsave(
    filename = paste0(enrich_dir, "GO_BP_barplot.png"),
    plot = barplot_BP,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  cat("GO BP barplot saved!\n")
  
  # --- Save results to CSV ---
  write.csv(
    as.data.frame(go_BP),
    paste0(enrich_dir, "GO_BP_results.csv"),
    row.names = FALSE
  )
  
  cat("GO results saved!\n")
  
  
  
  # ============================================
  # Complete GO Enrichment — AD vs Young
  # BP already done, now MF and CC
  # ============================================
  
  # GO Molecular Function
  go_MF <- enrichGO(
    gene          = sig_entrez_Young$ENTREZID,
    universe      = universe_Young$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # GO Cellular Component
  go_CC <- enrichGO(
    gene          = sig_entrez_Young$ENTREZID,
    universe      = universe_Young$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  cat("GO MF terms found:", nrow(go_MF), "\n")
  cat("GO CC terms found:", nrow(go_CC), "\n")
  
  # Preview top results if any found
  if(nrow(go_MF) > 0){
    cat("\nTop GO MF terms:\n")
    print(head(as.data.frame(go_MF)[,
                                    c("Description", "GeneRatio", "pvalue", "p.adjust")], 10))
  }
  
  if(nrow(go_CC) > 0){
    cat("\nTop GO CC terms:\n")
    print(head(as.data.frame(go_CC)[,
                                    c("Description", "GeneRatio", "pvalue", "p.adjust")], 10))
  }
  
  
  # ============================================
  # Visualize GO CC results
  # ============================================
  
  # Dotplot for CC
  dotplot_CC <- dotplot(go_CC,
                        showCategory = 15,
                        title = "GO Cellular Component — AD vs Young",
                        font.size = 10
  ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(dotplot_CC)
  
  ggsave(
    filename = paste0(enrich_dir, "GO_CC_dotplot.png"),
    plot = dotplot_CC,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Save CC results
  write.csv(
    as.data.frame(go_CC),
    paste0(enrich_dir, "GO_CC_results.csv"),
    row.names = FALSE
  )
  
  cat("GO CC results saved!\n")
  
  # ============================================
  # GO Enrichment — AD vs Old
  # Using relaxed cutoffs due to smaller gene list
  # ============================================
  
  # Relaxed significance threshold
  go_BP_Old <- enrichGO(
    gene          = sig_entrez$ENTREZID,
    universe      = universe_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.1,
    readable      = TRUE
  )
  
  go_CC_Old <- enrichGO(
    gene          = sig_entrez$ENTREZID,
    universe      = universe_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.1,
    readable      = TRUE
  )
  
  go_MF_Old <- enrichGO(
    gene          = sig_entrez$ENTREZID,
    universe      = universe_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.1,
    readable      = TRUE
  )
  
  cat("AD vs Old GO BP terms:", nrow(go_BP_Old), "\n")
  cat("AD vs Old GO CC terms:", nrow(go_CC_Old), "\n")
  cat("AD vs Old GO MF terms:", nrow(go_MF_Old), "\n")
  
  
  # ============================================
  # KEGG Pathway Analysis
  # ============================================
  
  # AD vs Young
  kegg_Young <- enrichKEGG(
    gene         = sig_entrez_Young$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  cat("KEGG pathways found — AD vs Young:", nrow(kegg_Young), "\n")
  
  if(nrow(kegg_Young) > 0){
    cat("\nTop KEGG pathways:\n")
    print(head(as.data.frame(kegg_Young)[,
                                         c("Description", "GeneRatio", "pvalue", "p.adjust")], 15))
  }
  
  # AD vs Old — relaxed cutoffs
  kegg_Old <- enrichKEGG(
    gene         = sig_entrez$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
  
  cat("\nKEGG pathways found — AD vs Old:", nrow(kegg_Old), "\n")
  
  if(nrow(kegg_Old) > 0){
    cat("\nTop KEGG pathways AD vs Old:\n")
    print(head(as.data.frame(kegg_Old)[,
                                       c("Description", "GeneRatio", "pvalue", "p.adjust")], 15))
  }
  
  
  # ============================================
  # Visualize KEGG results
  # ============================================
  
  # Dotplot
  dotplot_kegg <- dotplot(kegg_Young,
                          showCategory = 15,
                          title = "KEGG Pathways — AD vs Young",
                          font.size = 10
  ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(dotplot_kegg)
  
  ggsave(
    filename = paste0(enrich_dir, "KEGG_dotplot_AD_vs_Young.png"),
    plot = dotplot_kegg,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Save results
  write.csv(
    as.data.frame(kegg_Young),
    paste0(enrich_dir, "KEGG_AD_vs_Young_results.csv"),
    row.names = FALSE
  )
  
  cat("KEGG AD vs Young saved!\n")
  
  