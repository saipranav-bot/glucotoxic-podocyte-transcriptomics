## Repository Structure

```text
glucotoxic-podocyte-transcriptomics/
│
├── R/
│   └── glucotoxic_podocyte_analysis.R
│
├── data/
│   └── csv/
│       ├── ALL_SIGNIFICANT_DEGs_M_vs_H.csv
│       ├── ALL_SIGNIFICANT_DEGs_N_vs_H.csv
│       ├── Core_Shared_DEGs.csv
│       ├── DESeq2_M_vs_H_results.csv
│       ├── DESeq2_Medium_vs_High_shrunken.csv
│       ├── DESeq2_N_vs_H_results.csv
│       ├── DESeq2_Normal_vs_High_shrunken.csv
│       ├── GO_BP_top20.csv
│       ├── KEGG_enrichment_results.csv
│       ├── KEGG_top20.csv
│       ├── PCA_coordinates.csv
│       ├── SIGNIFICANT_DEGs_Medium_vs_High.csv
│       └── SIGNIFICANT_DEGs_Normal_vs_High.csv
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
│   │   └── Volcano_plot_FINAL.jpg
│   │
│   ├── 03_Enrichment/
│   │   ├── KEGG_dotplot_FINAL.jpg
│   │   ├── KEGG_network_FINAL.jpg
│   │   ├── GO_dotplot_FINAL.jpg
│   │   ├── GO_network_FINAL.jpg
│   │   └── GSEA_ridgeplot.jpg
│   │
│   ├── 04_Networks/
│   │   ├── STRING_global_network.jpg
│   │   ├── NOMO3_network.jpg
│   │   └── ER_MEMBRANE_PPI.jpg
│   │
│   └── 05_TF/
│       └── TF_activity.jpg
│
├── LINKEDIN_POST.md
├── README.md
└── .gitignore
```
