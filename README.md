# Glucotoxic Podocyte Transcriptomics

Transcriptomic dissection of glucotoxicity in human podocytes using the public
RNA-seq dataset **GEO GSE307956** (High 30 mM / Medium 15 mM / Normal 5.5 mM
glucose; n = 3 per group), with external validation in **GSE30528** (diabetic
kidney disease).

High glucose **suppresses** ribosomal and oxidative-phosphorylation
programmes; the reference level in all differential-expression contrasts is
**High** glucose (a negative log2FC therefore means the gene is higher in the
lower-glucose condition).

## Pipeline

`R/glucotoxic_podocyte_analysis.R` runs end-to-end:

1. DESeq2 differential expression (Medium vs High, Normal vs High) with apeglm
   shrinkage; broad (|log2FC| > 0.5) and strict (|log2FC| > 1.0) thresholds.
2. QC: PCA, sample-distance heatmap, MA plots.
3. ORA via **Enrichr** (`KEGG_2021_Human`, `GO_Biological_Process_2021`).
4. GSEA via **fgsea** on a `sign(log2FC) * -log10(padj)` ranking.
5. PPI network via **STRINGdb** (v12.5).
6. TF activity via **decoupleR** (CollecTRI + ULM) with a Medium-vs-High t-test.

## Repository structure

```text
glucotoxic-podocyte-transcriptomics/
│
├── R/
│   └── glucotoxic_podocyte_analysis.R
│
├── data/
│   ├── GSE307956_raw_count_matrix.csv
│   └── csv/
│       ├── ALL_SIGNIFICANT_DEGs_M_vs_H.csv
│       ├── ALL_SIGNIFICANT_DEGs_N_vs_H.csv
│       ├── Core_Shared_DEGs.csv
│       ├── DESeq2_M_vs_H_results.csv
│       ├── DESeq2_Medium_vs_High_shrunken.csv
│       ├── DESeq2_N_vs_H_results.csv
│       ├── DESeq2_Normal_vs_High_shrunken.csv
│       ├── GSEA_KEGG_results.csv
│       ├── GO_BP_top20.csv
│       ├── KEGG_enrichment_results.csv
│       ├── KEGG_top20.csv
│       ├── PCA_coordinates.csv
│       ├── SIGNIFICANT_DEGs_Medium_vs_High.csv
│       ├── SIGNIFICANT_DEGs_Normal_vs_High.csv
│       └── TF_activity_MvH.csv
│
├── figures/
│   ├── 01_QC/
│   │   ├── PCA_plot.jpeg
│   │   ├── Clean_Heatmap.jpeg
│   │   ├── MAplot_MH.jpeg
│   │   └── MAplot_NH.jpeg
│   │
│   ├── 02_DEG/
│   │   ├── Volcano_M_vs_H.jpeg
│   │   ├── Volcano_N_vs_H.jpeg
│   │   ├── Core_Gene_Heatmap.jpeg
│   │   └── Top_Gene_Expression.jpeg
│   │
│   ├── 03_Enrichment/
│   │   ├── KEGG_dotplot_FINAL.jpg
│   │   ├── GO_dotplot_FINAL.jpg
│   │   ├── GO_network_FINAL.jpg
│   │   ├── KEGG_network_FINAL.jpg
│   │   ├── KEGG_GSEA_dotplot.jpg
│   │   └── GSEA_ridgeplot.jpg
│   │
│   ├── 04_Networks/
│   │   └── STRING_global_network.jpg
│   │
│   └── 05_TF/
│       └── TF_activity.jpg
│
├── README.md
└── .gitignore
```
