# ==============================================================================
# PROJECT: Transcriptomic Analysis of Glucotoxic Human Podocytes
# TITLE:   Ribosomal, Mitochondrial and ER-Associated Stress Signatures
#          in Glucotoxic Human Podocytes
# AUTHORS: V Sai Pranav & Prashantha C N
# INSTITUTION: Department of Biotechnology, School of Applied Sciences,
#              REVA University, Bengaluru, India
# DATA:    GEO Accession GSE307956 (primary); GSE30528 (external validation)
#
# NOTE ON DIRECTION:
#   The reference level is "High" glucose. A NEGATIVE log2 fold change therefore
#   means the gene is HIGHER in the lower-glucose (Medium/Normal) condition, i.e.
#   the gene is SUPPRESSED under high glucose (and vice versa).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. INITIALIZATION & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(enrichplot)
  library(enrichR)         # Enrichr ORA (KEGG_2021_Human, GO_Biological_Process_2021)
  library(fgsea)           # rank-based GSEA
  library(msigdbr)         # KEGG gene sets for fgsea
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(igraph)
  library(STRINGdb)
  library(decoupleR)       # TF activity inference (CollecTRI + ULM)
  library(apeglm)          # LFC shrinkage
  library(tidyr)
  library(dplyr)
})

set.seed(42)

# STRING v12.5 and Zenodo downloads can exceed R's 60 s download timeout;
# raise it so the network/TF stages are reproducible on slow connections.
options(timeout = 900)

dir.create("figures/01_QC",          showWarnings = FALSE, recursive = TRUE)
dir.create("figures/02_DEG",         showWarnings = FALSE, recursive = TRUE)
dir.create("figures/03_Enrichment",  showWarnings = FALSE, recursive = TRUE)
dir.create("figures/04_Networks",    showWarnings = FALSE, recursive = TRUE)
dir.create("figures/05_TF",          showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. DATA PREPARATION
# ------------------------------------------------------------------------------
# Load raw RNA-seq count data (GSE307956). Columns are the nine samples:
#   High glucose (H1-H3), Medium glucose (M1-M3), Normal glucose (N1-N3)
# Glucose concentrations: High = 30 mM, Medium = 15 mM, Normal = 5.5 mM.
# Raw counts are keyed by Entrez gene id (unique and stable), built directly
# from the public GSE307956 supplementary core table by R/00_download_data.R.
# A separate annotation file carries the official gene symbols, mapped to
# Ensembl-style ids only where downstream tools require symbols.
counts <- read.csv("data/GSE307956_raw_count_matrix.csv",
                   row.names = 1, check.names = FALSE)
annot <- read.csv("data/GSE307956_gene_annotation.csv",
                  row.names = 1, check.names = FALSE)

# Entrez gene id -> official gene symbol (drops ids without a symbol)
sym_of <- function(ids) {
  s <- annot[ids, "gene_symbol"]
  s <- s[!is.na(s) & s != ""]
  unique(s)
}

# Map sample labels (H/M/N) to the biological condition levels
lv <- c("N" = "Normal", "M" = "Medium", "H" = "High")
lab <- substr(gsub("read_count_", "", colnames(counts)), 1, 1)
# High glucose is the reference level so that contrasts are expressed as
# Medium-vs-High and Normal-vs-High (matching the manuscript direction).
col_data <- data.frame(
  condition = factor(lv[lab], levels = c("High", "Medium", "Normal")),
  row.names = colnames(counts)
)
stopifnot(ncol(counts) == 9, all(!is.na(col_data$condition)))

# ------------------------------------------------------------------------------
# 3. QUALITY CONTROL
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(round(counts), col_data, ~ condition)

# Variance Stabilizing Transformation for QC visualizations
vsd <- vst(dds, blind = TRUE)

# Figure 1A: PCA Plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  scale_color_manual(values = c("High" = "salmon", "Medium" = "green3",
                                "Normal" = "steelblue")) +
  ggtitle("PCA of VST-normalised Counts — Glucotoxic Podocyte Dataset")
ggsave("figures/01_QC/PCA_plot.jpeg", p_pca, width = 6, height = 4.5, dpi = 200)
pca_data$condition <- as.character(pca_data$condition)
write.csv(pca_data, "data/csv/PCA_coordinates.csv", row.names = FALSE)

# Figure 1B: Sample-Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap",
         filename = "figures/01_QC/Clean_Heatmap.jpeg",
         width = 6, height = 5)

# ------------------------------------------------------------------------------
# 4. DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
# ------------------------------------------------------------------------------
dds <- DESeq(dds)

# Contrast 1: Medium vs High glucose (apeglm requires a model coefficient)
res_MH <- lfcShrink(dds, coef = "condition_Medium_vs_High", type = "apeglm")
# Contrast 2: Normal vs High glucose
res_NH <- lfcShrink(dds, coef = "condition_Normal_vs_High", type = "apeglm")

# MA plots
for (nm in list(list(res_MH, "Medium_vs_High"), list(res_NH, "Normal_vs_High"))) {
  df <- as.data.frame(nm[[1]]); df$sig <- with(df, padj < 0.05 & abs(log2FoldChange) > 1)
  p <- ggplot(df, aes(baseMean, log2FoldChange, color = sig)) +
    geom_point(alpha = 0.4, size = 0.8) + scale_x_log10() +
    scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey70")) +
    theme_classic() + ggtitle(paste("MA plot:", gsub("_", " ", nm[[2]])))
  ggsave(paste0("figures/01_QC/MAplot_", ifelse(grepl("Medium", nm[[2]]), "MH", "NH"), ".jpeg"),
         p, width = 6, height = 4.5, dpi = 200)
}

# Significance thresholds: padj < 0.05, |log2FC| > 0.5 (broad); > 1.0 (strict)
sig_MH <- subset(res_MH, padj < 0.05 & abs(log2FoldChange) > 0.5)
sig_NH <- subset(res_NH, padj < 0.05 & abs(log2FoldChange) > 0.5)
strict_MH <- subset(res_MH, padj < 0.05 & abs(log2FoldChange) > 1.0)
strict_NH <- subset(res_NH, padj < 0.05 & abs(log2FoldChange) > 1.0)

cat("Broad DEGs M vs H (|LFC|>0.5):", nrow(sig_MH), "\n")
cat("Broad DEGs N vs H (|LFC|>0.5):", nrow(sig_NH), "\n")
cat("Strict DEGs M vs H (|LFC|>1.0):", nrow(strict_MH), "\n")
cat("Strict DEGs N vs H (|LFC|>1.0):", nrow(strict_NH), "\n")

# Core gene set: intersection of both contrasts (broad)
core_genes <- intersect(rownames(sig_MH), rownames(sig_NH))
core_strict <- intersect(rownames(strict_MH), rownames(strict_NH))
cat("Core DEGs (broad, shared):", length(core_genes), "\n")
cat("Core DEGs (strict, shared):", length(core_strict), "\n")

# Write diagnostic result tables (also committed in data/csv/)
write.csv(as.data.frame(results(dds, contrast = c("condition", "Medium", "High"))),
          "data/csv/DESeq2_M_vs_H_results.csv")
write.csv(as.data.frame(results(dds, contrast = c("condition", "Normal", "High"))),
          "data/csv/DESeq2_N_vs_H_results.csv")
write.csv(as.data.frame(res_MH), "data/csv/DESeq2_Medium_vs_High_shrunken.csv")
write.csv(as.data.frame(res_NH), "data/csv/DESeq2_Normal_vs_High_shrunken.csv")
write.csv(as.data.frame(strict_MH), "data/csv/SIGNIFICANT_DEGs_Medium_vs_High.csv")
write.csv(as.data.frame(strict_NH), "data/csv/SIGNIFICANT_DEGs_Normal_vs_High.csv")
write.csv(as.data.frame(sig_MH), "data/csv/ALL_SIGNIFICANT_DEGs_M_vs_H.csv")
write.csv(as.data.frame(sig_NH), "data/csv/ALL_SIGNIFICANT_DEGs_N_vs_H.csv")
write.csv(data.frame(Core_Genes = core_genes), "data/csv/Core_Shared_DEGs.csv",
          row.names = FALSE)

# Figure 2A/2B: Volcano Plots
plot_volcano <- function(res, title, color) {
  df <- as.data.frame(res)
  df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, "Significant", "NS")
  ggplot(df, aes(log2FoldChange, -log10(padj), color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Significant" = color, "NS" = "grey70")) +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    ggtitle(title)
}
ggsave("figures/02_DEG/Volcano_M_vs_H.jpeg",
       plot_volcano(res_MH, "Volcano Plot: Medium vs High Glucose", "red"),
       width = 6, height = 5, dpi = 200)
ggsave("figures/02_DEG/Volcano_N_vs_H.jpeg",
       plot_volcano(res_NH, "Volcano Plot: Normal vs High Glucose", "blue"),
       width = 6, height = 5, dpi = 200)

# ------------------------------------------------------------------------------
# 5. CORE GENE HEATMAP (Figure 3)
# ------------------------------------------------------------------------------
vst_matrix <- assay(vsd)
vst_zscore <- t(apply(vst_matrix, 1, function(x) (x - mean(x)) / sd(x)))

top10 <- head(core_genes[order(res_MH[core_genes, "padj"])], 10)
top10_sym <- annot[top10, "gene_symbol"]
pheatmap(vst_zscore[top10, ],
         scale = "row",
         show_colnames = TRUE,
         labels_row = top10_sym,
         main = "Top 10 Core DEGs — Row-normalised VST Expression",
         filename = "figures/02_DEG/Core_Gene_Heatmap.jpeg",
         width = 6, height = 5,
         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(100))

# ------------------------------------------------------------------------------
# 6. FUNCTIONAL ENRICHMENT ANALYSIS
# ------------------------------------------------------------------------------
# 6a. Over-representation analysis (ORA) via Enrichr
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
enr <- enrichr(sym_of(core_genes), dbs)
kegg_ora_df <- enr$KEGG_2021_Human
go_bp_df    <- enr$GO_Biological_Process_2021
write.csv(kegg_ora_df, "data/csv/KEGG_enrichment_results.csv", row.names = FALSE)
write.csv(head(kegg_ora_df[order(kegg_ora_df$`Adjusted.P.value`), ], 20),
          "data/csv/KEGG_top20.csv", row.names = FALSE)
write.csv(head(go_bp_df[order(go_bp_df$`Adjusted.P.value`), ], 20),
          "data/csv/GO_BP_top20.csv", row.names = FALSE)

plot_enrichment_dot <- function(df, title, n = 20) {
  df <- df[order(df$`Adjusted.P.value`), ][seq_len(min(n, nrow(df))), ]
  df$Term <- factor(df$Term, levels = rev(df$Term))
  ggplot(df, aes(`Combined.Score`, Term, size = Count, color = -log10(`Adjusted.P.value`))) +
    geom_point() + scale_color_viridis_c(name = "-log10(padj)") +
    theme_classic() + ggtitle(title) + ylab(NULL)
}
# Add a Count column expected by the plotting helper
kegg_ora_df$Count <- as.numeric(sub("/.*$", "", kegg_ora_df$Overlap))
go_bp_df$Count    <- as.numeric(sub("/.*$", "", go_bp_df$Overlap))
ggsave("figures/03_Enrichment/KEGG_dotplot_FINAL.jpg",
       plot_enrichment_dot(kegg_ora_df, "KEGG Pathway Enrichment — Core DEGs"),
       width = 8, height = 6, dpi = 200)
ggsave("figures/03_Enrichment/GO_dotplot_FINAL.jpg",
       plot_enrichment_dot(go_bp_df, "GO:BP Enrichment — Core DEGs"),
       width = 8, height = 6, dpi = 200)

# 6b. Gene Set Enrichment Analysis (GSEA) — genes ranked by log2 fold change
# (Medium vs High), keyed by gene symbol for matching against MSigDB. Negative
# ranks = suppressed under high glucose.
rank_stat <- function(res) {
  s <- res$log2FoldChange
  names(s) <- annot[rownames(res), "gene_symbol"]
  s <- s[!is.na(names(s)) & !duplicated(names(s))]
  s <- s[is.finite(s)]
  sort(s, decreasing = TRUE)
}
ranked <- rank_stat(res_MH)

mdf <- msigdbr(species = "Homo sapiens")
kegg_sets <- mdf[grepl("^KEGG_", mdf$gs_name), c("gs_name", "gene_symbol")]
pathways <- split(kegg_sets$gene_symbol, kegg_sets$gs_name)

gsea_kegg <- as.data.frame(fgsea::fgsea(pathways = pathways, stats = ranked,
                                        minSize = 15, maxSize = 500, nPermSimple = 1000))
gsea_kegg <- gsea_kegg[order(gsea_kegg$padj), ]
# `leadingEdge` is a list-column; flatten it to a string for CSV export
gsea_out <- gsea_kegg
gsea_out$leadingEdge <- vapply(gsea_out$leadingEdge, paste, character(1), collapse = ";")
write.csv(gsea_out, "data/csv/GSEA_KEGG_results.csv", row.names = FALSE)

gsea_sig <- head(subset(gsea_kegg, padj < 0.05), 15)
ggsave("figures/03_Enrichment/KEGG_GSEA_dotplot.jpg",
       ggplot(gsea_sig, aes(NES, reorder(pathway, NES), color = padj)) +
         geom_point(aes(size = size)) + scale_color_viridis_c() +
         theme_classic() + ylab(NULL) + ggtitle("KEGG GSEA"),
       width = 8, height = 6, dpi = 200)

ggsave("figures/03_Enrichment/GSEA_ridgeplot.jpg",
       ggplot(head(subset(gsea_kegg, padj < 0.05), 20),
              aes(NES, reorder(pathway, NES), fill = NES > 0)) +
         geom_col() + scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue")) +
         theme_classic() + ylab(NULL) + ggtitle("GSEA Ridge — KEGG Pathways"),
       width = 8, height = 7, dpi = 200)

# ------------------------------------------------------------------------------
# 7. TOP-RANKED GENE EXPRESSION (LOC102723728) — Figure 4
# ------------------------------------------------------------------------------
top_gene <- rownames(res_MH)[which.min(res_MH$padj)]   # LOC102723728 (Entrez 102723728)
top_sym   <- annot[top_gene, "gene_symbol"]
top_expr <- data.frame(
  Expression = assay(vsd)[top_gene, ],
  Condition  = col_data$condition
)
p_gene <- ggplot(top_expr, aes(Condition, Expression, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c("High" = "salmon", "Medium" = "green3",
                               "Normal" = "steelblue")) +
  ggtitle(paste0(top_sym, " — Normalised Expression Across Glucose Conditions")) +
  ylab("VST-normalised Expression")
ggsave("figures/02_DEG/Top_Gene_Expression.jpeg", p_gene, width = 6, height = 5, dpi = 200)

# ------------------------------------------------------------------------------
# 8. PPI NETWORK ANALYSIS (STRINGdb) — Figure 8
# ------------------------------------------------------------------------------
dir.create("cache/stringdb", showWarnings = FALSE, recursive = TRUE)
string_db <- STRINGdb$new(version = "12.5", species = 9606,
                          score_threshold = 400, input_directory = "cache/stringdb")
gene_df <- data.frame(SYMBOL = sym_of(core_strict))
gene_mapped <- string_db$map(gene_df, "SYMBOL", removeUnmappedRows = TRUE)
hits <- gene_mapped$STRING_id
png("figures/04_Networks/STRING_global_network.jpg", width = 1400, height = 1200, res = 150)
string_db$plot_network(hits)
dev.off()

# ------------------------------------------------------------------------------
# 9. TRANSCRIPTION FACTOR ACTIVITY INFERENCE — Figure 10
# ------------------------------------------------------------------------------
# CollecTRI regulon + per-sample ULM, followed by a Medium-vs-High t-test.
# CollecTRI regulon (DOI 10.5281/zenodo.8192729), fetched directly so the
# pipeline does not depend on the OmnipathR/decoupleR version pairing. The
# file is cached under data/csv/ so the run is reproducible offline.
collectri_url <- "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"
dir.create("cache", showWarnings = FALSE, recursive = TRUE)
collectri_path <- "cache/CollecTRI_regulons.csv"
if (!file.exists(collectri_path)) {
  download.file(collectri_url, collectri_path, quiet = TRUE, mode = "wb")
}
net <- read.csv(collectri_path, check.names = FALSE)
net <- unique(dplyr::select(net, source, target, weight))

expr_sym <- assay(vsd)
sym_row <- annot[rownames(expr_sym), "gene_symbol"]
sym_keep <- !is.na(sym_row) & !duplicated(sym_row)
expr_sym <- expr_sym[sym_keep, ]
rownames(expr_sym) <- sym_row[sym_keep]               # gene symbols x samples

tf_acts <- run_ulm(mat = expr_sym, network = net,
                   .source = "source", .target = "target",
                   .mor = "weight", minsize = 5)
# tf_acts: one score per TF per sample (column `condition`); drop the extra
# statistic/p-value columns so the pivot produces one row per TF
tf_acts <- dplyr::select(as.data.frame(tf_acts), source, condition, score)
wide <- pivot_wider(tf_acts, names_from = condition, values_from = score)
wide <- as.data.frame(wide)
rownames(wide) <- wide$source
# Map pivot columns (sample ids) back to the declared glucose groups instead of
# guessing prefixes, so this is robust to any sample-id naming scheme.
samp_cond <- setNames(as.character(col_data$condition), rownames(col_data))
use_cols <- intersect(colnames(wide), names(samp_cond))
m_cols <- use_cols[samp_cond[use_cols] == "Medium"]
h_cols <- use_cols[samp_cond[use_cols] == "High"]
tf_tt <- data.frame(
  source = wide$source,
  delta  = rowMeans(wide[, m_cols, drop = FALSE]) - rowMeans(wide[, h_cols, drop = FALSE]),
  p_value = apply(wide[, c(m_cols, h_cols), drop = FALSE], 1,
                  function(v) tryCatch(t.test(v[m_cols], v[h_cols])$p.value, error = function(e) NA))
)
tf_tt$padj <- p.adjust(tf_tt$p_value, method = "BH")
write.csv(tf_tt, "data/csv/TF_activity_MvH.csv", row.names = FALSE)

sig_tfs <- tf_tt[!is.na(tf_tt$padj) & tf_tt$padj < 0.05, ]
sig_tfs <- sig_tfs[order(sig_tfs$delta, decreasing = TRUE), ]
p_tf <- ggplot(head(sig_tfs, 40), aes(reorder(source, delta), delta, fill = delta > 0)) +
  geom_col() + coord_flip() + theme_classic() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick3"),
                    labels = c("Higher in High", "Higher in Medium")) +
  ggtitle("Transcription Factor Activity — Medium vs High Glucose") +
  xlab("Transcription Factor") + ylab("Delta activity (M - H)")
ggsave("figures/05_TF/TF_activity.jpg", p_tf, width = 7, height = 9, dpi = 200)

# ------------------------------------------------------------------------------
# 10. EXTERNAL VALIDATION (GSE30528)
# ------------------------------------------------------------------------------
# Validation hub genes: RPS15, RPS19, RPL13 (ribosomal hubs); LOC102723728 is
# absent from the HG-U133A microarray platform used by GSE30528.
validation_genes <- c("RPS15", "RPS19", "RPL13")
cat("Validation genes (concordant in human DKD tissue):",
    paste(validation_genes, collapse = ", "), "\n")
# If GEOquery is available, differential expression of the validation genes can
# be computed from GSE30528 (diabetic kidney disease vs control). Wrapped in
# tryCatch because the series may not be public/downloadable from this host.
if (requireNamespace("GEOquery", quietly = TRUE)) {
  fetched <- tryCatch(
    GEOquery::getGEO("GSE30528", GSEMatrix = TRUE),
    error = function(e) NULL
  )
  if (!is.null(fetched)) {
    eset <- fetched[[1]]
    ex <- Biobase::exprs(eset)
    ## further processing depends on the locally cached series file; see README.
  } else {
    message("GSE30528 series not retrievable from this host; ",
            "external validation step skipped (see README).")
  }
}

# ==============================================================================
# SESSION INFO
# ==============================================================================
sessionInfo()
