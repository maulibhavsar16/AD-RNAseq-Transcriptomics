# Transcriptomic Rewiring in Alzheimer's Disease

## Project Summary
End-to-end bulk RNA-seq analysis identifying transcriptomic signatures,
co-expression networks, transcription factor activity and genetic
overlap in Alzheimer's Disease human brain tissue.

## Dataset
- **GEO Accession:** GSE159699
- **Publication:** Nativio et al. 2020
- **Samples:** 30 postmortem lateral temporal lobe samples
- **Groups:** AD (n=12), Old Control (n=10), Young Control (n=8)

## Key Findings
1. Oxidative phosphorylation most consistently downregulated process
2. 3,929 DEGs in AD vs Young; 684 AD-specific beyond aging
3. KEGG Alzheimer's Disease pathway independently enriched (padj=8.58e-10)
4. Yellow WGCNA module correlated with AD (r=0.59); top hub: SHANK3
5. CREB1 and REST identified as key dysregulated TFs
6. 8 AD GWAS genes in DEG list including APP, PSEN1, BIN1, TOMM40

## Pipeline
FastQC/fastp → Salmon → DESeq2 → clusterProfiler/fgsea → WGCNA → decoupleR → biomaRt

## Author
Mauli Bhavsar | MSc Bioinformatics 2025
