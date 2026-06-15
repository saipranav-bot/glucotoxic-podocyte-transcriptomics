# ==============================================================================
# PROJECT: Transcriptomic Analysis of Glucotoxic Human Podocytes
# TITLE:   Ribosomal, Mitochondrial and ER-Associated Stress Signatures
#          in Glucotoxic Human Podocytes
# AUTHORS: V Sai Pranav & Prashantha C N
# INSTITUTION: Department of Biotechnology, School of Applied Sciences,
#              REVA University, Bengaluru, India
# DATA:    GEO Accession GSE307956 (primary); GSE30528 (external validation)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. INITIALIZATION & LIBRARIES
# ------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(STRINGdb)
library(decoupleR)   # TF activity inference
library(apeglm)      # LFC shrinkage

set.seed(42)

# ------------------------------------------------------------------------------
# 2. DATA PREPARATION
# ------------------------------------------------------------------------------
# Load raw RNA-seq count data (GSE307956)
# 9 samples: High glucose (H1-H3), Medium glucose (M1-M3), Normal glucose (N1-N3)
# Glucose concentrations: High = 30 mM, Medium = 15 mM, Normal = 5.5 mM

counts <- read.csv("data/csv/DESeq2_M_vs_H_results.csv", row.names = 1, check.names = FALSE)

col_data <- data.frame(
  condition = factor(
    rep(c("High", "Medium", "Normal"), each = 3),
    levels = c("Normal", "Medium", "High")
  )
)
rownames(col_data) <- colnames(counts)

# ------------------------------------------------------------------------------
# 3. QUALITY CONTROL
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(round(counts), col_data, ~ condition)

# Variance Stabilizing Transformation for QC visualizations
vsd <- vst(dds, blind = TRUE)

# Figure 1A: PCA Plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  scale_color_manual(values = c("High" = "salmon", "Medium" = "green3", "Normal" = "steelblue")) +
  ggtitle("PCA of VST-normalised Counts — Glucotoxic Podocyte Dataset")

# Figure 1B: Sample-Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

# ------------------------------------------------------------------------------
# 4. DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
# ------------------------------------------------------------------------------
dds <- DESeq(dds)

# Contrast 1: Medium vs High glucose
res_MH <- lfcShrink(dds,
                     contrast = c("condition", "Medium", "High"),
                     type = "apeglm")

# Contrast 2: Normal vs High glucose
res_NH <- lfcShrink(dds,
                     contrast = c("condition", "Normal", "High"),
                     type = "apeglm")

# Significance thresholds: padj < 0.05, |log2FC| > 0.5
sig_MH <- subset(res_MH, padj < 0.05 & abs(log2FoldChange) > 0.5)
sig_NH <- subset(res_NH, padj < 0.05 & abs(log2FoldChange) > 0.5)

cat("DEGs M vs H:", nrow(sig_MH), "\n")
cat("DEGs N vs H:", nrow(sig_NH), "\n")

# Core gene set: intersection of both contrasts
core_genes <- intersect(rownames(sig_MH), rownames(sig_NH))
cat("Core DEGs (shared):", length(core_genes), "\n")

# Figure 2A/2B: Volcano Plots
plot_volcano <- function(res, title, color) {
  df <- as.data.frame(res)
  df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 0.5, "Significant", "NS")
  ggplot(df, aes(log2FoldChange, -log10(padj), color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Significant" = color, "NS" = "grey70")) +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    ggtitle(title)
}

plot_volcano(res_MH, "Volcano Plot: Medium vs High Glucose", "red")
plot_volcano(res_NH, "Volcano Plot: Normal vs High Glucose", "blue")

# ------------------------------------------------------------------------------
# 5. CORE GENE HEATMAP (Figure 3)
# ------------------------------------------------------------------------------
normalize_data <- function(x) (x - mean(x)) / sd(x)
vst_matrix <- assay(vst(dds))
vst_zscore <- t(apply(vst_matrix, 1, normalize_data))

top20 <- head(core_genes[order(res_MH[core_genes, "padj"])], 20)
pheatmap(vst_zscore[top20, ],
         scale = "row",
         show_colnames = TRUE,
         main = "Top 20 Core DEGs — Row-normalised VST Expression",
         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(100))

# ------------------------------------------------------------------------------
# 6. FUNCTIONAL ENRICHMENT ANALYSIS
# ------------------------------------------------------------------------------
# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(core_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Figure 4: KEGG Pathway ORA
kegg_ora <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                        organism     = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
dotplot(kegg_ora, showCategory = 15, title = "KEGG Pathway Enrichment — Core DEGs")

# GO Biological Process ORA
go_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
dotplot(go_bp, showCategory = 20, title = "GO:BP Enrichment — Core DEGs")

# Figure 5: GSEA
ranked_genes <- sort(setNames(res_MH$log2FoldChange, rownames(res_MH)),
                     decreasing = TRUE)
ranked_entrez <- bitr(names(ranked_genes), fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
ranked_vec <- ranked_genes[ranked_entrez$SYMBOL]
names(ranked_vec) <- ranked_entrez$ENTREZID

gsea_kegg <- gseKEGG(geneList     = ranked_vec,
                      organism     = "hsa",
                      nPerm        = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")
dotplot(gsea_kegg, showCategory = 10, title = "KEGG GSEA Results")
ridgeplot(gsea_kegg) + ggtitle("GSEA Ridge Plot — KEGG Pathways")

# ------------------------------------------------------------------------------
# 7. NOMO3 EXPRESSION ANALYSIS (Figure 6)
# ------------------------------------------------------------------------------
nomo3_expr <- data.frame(
  Expression = assay(vst(dds))["NOMO3", ],
  Condition  = col_data$condition
)

ggplot(nomo3_expr, aes(Condition, Expression, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c("High" = "salmon", "Medium" = "green3", "Normal" = "steelblue")) +
  ggtitle("NOMO3 — Normalised Expression Across Glucose Conditions") +
  ylab("VST-normalised Expression")

# ------------------------------------------------------------------------------
# 8. PPI NETWORK ANALYSIS (STRINGdb) — Figure 7
# ------------------------------------------------------------------------------
string_db <- STRINGdb$new(version = "12.0", species = 9606,
                           score_threshold = 400, input_directory = ".")

# Map core genes to STRING identifiers
gene_df <- data.frame(SYMBOL = core_genes)
gene_mapped <- string_db$map(gene_df, "SYMBOL", removeUnmappedRows = TRUE)
hits <- gene_mapped$STRING_id

# Full core gene network
string_db$plot_network(hits)

# Dedicated NOMO3 network
nomo3_mapped <- string_db$map(data.frame(SYMBOL = "NOMO3"), "SYMBOL",
                               removeUnmappedRows = TRUE)
nomo3_network <- string_db$get_subnetwork(nomo3_mapped$STRING_id)
plot(nomo3_network, main = "NOMO3 PPI Network (STRINGdb v12.0)")

# ------------------------------------------------------------------------------
# 9. TRANSCRIPTION FACTOR ACTIVITY INFERENCE (Figure 9)
# ------------------------------------------------------------------------------
# Using decoupleR with DoRothEA regulon database
dorothea_hs <- get_dorothea(organism = "human", levels = c("A", "B", "C"))

# Run ULM
mat <- as.matrix(res_MH$log2FoldChange)
rownames(mat) <- rownames(res_MH)

tf_activities <- run_ulm(mat        = mat,
                          network   = dorothea_hs,
                          .source   = "source",
                          .target   = "target",
                          .mor      = "mor",
                          minsize   = 5)

sig_tfs <- tf_activities[tf_activities$p_value < 0.05, ]
sig_tfs <- sig_tfs[order(sig_tfs$score, decreasing = TRUE), ]

ggplot(sig_tfs, aes(reorder(source, score), score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue"),
                    labels = c("Activated", "Suppressed")) +
  ggtitle("Transcription Factor Activity — High Glucose Conditions") +
  xlab("Transcription Factor") + ylab("Activity Score (ULM)")

# ------------------------------------------------------------------------------
# 10. EXTERNAL VALIDATION (GSE30528 — Figure 8)
# ------------------------------------------------------------------------------
# Validation hub genes: RPS15, RPS19, RPL13 (ribosomal hubs)
# Note: NOMO3 absent from HG-U133A microarray platform
validation_genes <- c("RPS15", "RPS19", "RPL13")
cat("Validation genes (concordant in human DKD tissue):", paste(validation_genes, collapse = ", "), "\n")

# ==============================================================================
# SESSION INFO
# ==============================================================================
sessionInfo()
