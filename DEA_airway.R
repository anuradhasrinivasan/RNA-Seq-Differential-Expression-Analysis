# =============================================================================
# RNA-Seq Differential Expression Analysis
# Dataset: Airway Smooth Muscle Cells (Himes et al., 2014)
# Author: Anuradha Srinivasan
# Date: April 2026
# Description: Complete DEA pipeline using DESeq2 and Bioconductor
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 1: Load Required Libraries
# -----------------------------------------------------------------------------
library(DESeq2)
library(airway)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

# -----------------------------------------------------------------------------
# STEP 2: Load and Explore Dataset
# -----------------------------------------------------------------------------

# Load airway dataset
data(airway)

# Check dataset dimensions
cat("Dataset dimensions:", dim(airway), "\n")
cat("Total genes:", nrow(airway), "\n")
cat("Total samples:", ncol(airway), "\n")

# View sample information
colData(airway)

# View first few rows of count matrix
head(assay(airway))

# -----------------------------------------------------------------------------
# STEP 3: Create DESeq2 Object
# -----------------------------------------------------------------------------

# Set up DESeq2 object with cell + dex design
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# Set untreated as reference level
dds$dex <- relevel(dds$dex, ref = "untrt")

# Filter low count genes (keep genes with at least 10 counts total)
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after filtering:", nrow(dds), "\n")

# -----------------------------------------------------------------------------
# STEP 4: Run DESeq2 Analysis
# -----------------------------------------------------------------------------

# Run DESeq2
dds <- DESeq(dds)

# Get results for treated vs untreated
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Summary of results
summary(res)

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# View top 10 DEGs
head(res_ordered, 10)

# -----------------------------------------------------------------------------
# STEP 5: Add Gene Symbols to Results
# -----------------------------------------------------------------------------

# Map Ensembl IDs to gene symbols
res_ordered$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_ordered),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Map Entrez IDs for enrichment analysis
res_ordered$entrez <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_ordered),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# View results with gene symbols
head(res_ordered)

# -----------------------------------------------------------------------------
# STEP 6: Create Results Directory
# -----------------------------------------------------------------------------
dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# STEP 7: PCA Plot
# -----------------------------------------------------------------------------

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca_plot <- plotPCA(vsd, intgroup = "dex") +
  ggtitle("PCA Plot - Airway Dataset") +
  theme_bw() +
  scale_color_manual(
    values = c("untrt" = "#2196F3", "trt" = "#F44336"),
    labels = c("Untreated", "Treated"),
    name = "Treatment"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(pca_plot)
ggsave("results/PCA_airway.png", pca_plot, width = 7, height = 5, dpi = 300)
cat("PCA plot saved!\n")

# -----------------------------------------------------------------------------
# STEP 8: Volcano Plot
# -----------------------------------------------------------------------------

# Convert to dataframe for volcano plot
res_df <- as.data.frame(res_ordered)
res_df$symbol[is.na(res_df$symbol)] <- rownames(res_df)[is.na(res_df$symbol)]

# Create volcano plot
png("results/volcano_airway.png", width = 10, height = 8, units = "in", res = 300)
EnhancedVolcano(
  res_df,
  lab = res_df$symbol,
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano Plot - Dexamethasone Treated vs Untreated",
  subtitle = "Airway Smooth Muscle Cells",
  pCutoff = 0.05,
  FCcutoff = 1.5,
  pointSize = 2.0,
  labSize = 3.5,
  col = c("grey30", "#2196F3", "#4CAF50", "#F44336"),
  legendLabels = c("NS", "Log2FC", "p-value", "p-value & Log2FC"),
  legendPosition = "right",
  drawConnectors = TRUE,
  widthConnectors = 0.5
)
dev.off()
cat("Volcano plot saved!\n")

# -----------------------------------------------------------------------------
# STEP 9: Heatmap - Top 30 DEGs
# -----------------------------------------------------------------------------

# Get top 30 significant DEGs
sig_genes <- res_ordered[!is.na(res_ordered$padj) &
                           res_ordered$padj < 0.05 &
                           abs(res_ordered$log2FoldChange) > 1.5, ]

top30 <- head(sig_genes, 30)

# Extract normalized counts for top 30 genes
top30_counts <- assay(vsd)[rownames(top30), ]

# Replace Ensembl IDs with gene symbols
rownames(top30_counts) <- ifelse(
  !is.na(top30$symbol),
  top30$symbol,
  rownames(top30)
)

# Create annotation for columns
annotation_col <- data.frame(
  Treatment = colData(dds)$dex,
  row.names = colnames(top30_counts)
)

# Define annotation colors
ann_colors <- list(
  Treatment = c(untrt = "#2196F3", trt = "#F44336")
)

# Plot heatmap
png("results/heatmap_top30.png", width = 10, height = 10, units = "in", res = 300)
pheatmap(
  top30_counts,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap - Top 30 Differentially Expressed Genes",
  fontsize = 10,
  fontsize_row = 9,
  color = colorRampPalette(c("#2196F3", "white", "#F44336"))(100)
)
dev.off()
cat("Heatmap saved!\n")

# -----------------------------------------------------------------------------
# STEP 10: GO Enrichment Analysis
# -----------------------------------------------------------------------------

# Get significant DEG Entrez IDs
sig_entrez <- sig_genes$entrez[!is.na(sig_genes$entrez)]

# Run GO enrichment - Biological Process
go_results <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot GO dotplot
go_plot <- dotplot(go_results, showCategory = 15) +
  ggtitle("GO Enrichment - Biological Process") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/GO_dotplot.png", go_plot, width = 10, height = 8, dpi = 300)
cat("GO enrichment plot saved!\n")

# -----------------------------------------------------------------------------
# STEP 11: KEGG Pathway Analysis
# -----------------------------------------------------------------------------

# Run KEGG enrichment
kegg_results <- enrichKEGG(
  gene = sig_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Plot KEGG barplot
kegg_plot <- barplot(kegg_results, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment Analysis") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/KEGG_barplot.png", kegg_plot, width = 10, height = 8, dpi = 300)
cat("KEGG pathway plot saved!\n")

# -----------------------------------------------------------------------------
# STEP 12: Export Results to CSV
# -----------------------------------------------------------------------------

# Export all results
write.csv(
  as.data.frame(res_ordered),
  "results/all_DEG_results.csv",
  row.names = TRUE
)

# Export only significant DEGs
write.csv(
  as.data.frame(sig_genes),
  "results/significant_DEGs.csv",
  row.names = TRUE
)

cat("Results exported to CSV!\n")

# -----------------------------------------------------------------------------
# STEP 13: Session Info (for reproducibility)
# -----------------------------------------------------------------------------
sessionInfo()

# =============================================================================
# END OF SCRIPT
# =============================================================================
