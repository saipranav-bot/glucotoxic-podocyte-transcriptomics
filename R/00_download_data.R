# ==============================================================================
# DOWNLOAD PUBLIC DATA AND BUILD INPUT MATRICES — GSE307956
# ==============================================================================
# Rebuilds every input strictly from the public GEO record (GSE307956,
# "Glucotoxic podocytes — core table", GEO series supplementary file
# GSE307956_core_table_gene.txt.gz). No committed data is used as a source.
#
# Outputs (written to data/):
#   GSE307956_raw_count_matrix.csv    Entrez gene id x 9 samples (read counts)
#   GSE307956_gene_annotation.csv     Entrez gene id -> official gene symbol
#
# The raw 120 MB source file is cached under data/cache/ and is git-ignored.
# ==============================================================================

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307956/suppl/",
  "GSE307956_core_table_gene.txt.gz"
)

dir.create("data/cache", showWarnings = FALSE, recursive = TRUE)
srce <- "data/cache/GSE307956_core_table_gene.txt.gz"

if (!file.exists(srce)) {
  cat("Downloading public GSE307956 core table (~120 MB)...\n")
  options(timeout = 900)
  download.file(url, srce, mode = "wb", quiet = TRUE)
}

tab <- read.table(gzfile(srce), header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "", comment.char = "")

rc <- c("read_count_H1", "read_count_H2", "read_count_H3",
        "read_count_M1", "read_count_M2", "read_count_M3",
        "read_count_N1", "read_count_N2", "read_count_N3")

sub <- tab[!is.na(tab$gene_id) & tab$gene_id != "", c("gene_id", "gene_symbol", rc)]
sub$gene_id <- as.character(sub$gene_id)

# Drop spec-supplied symbol strings that are not real symbols (missing values
# and collapsed GO terms), then keep one row per Entrez id.
is_junk <- is.na(sub$gene_symbol) |
  grepl("^GO(BP|MF|CC)_", sub$gene_symbol) |
  grepl("///", sub$gene_symbol, fixed = TRUE)
sub <- sub[!is_junk, ]
sub <- sub[!duplicated(sub$gene_id), ]

mat <- sub[, rc]
rownames(mat) <- sub$gene_id
mat <- as.matrix(mat)
storage.mode(mat) <- "integer"
stopifnot(ncol(mat) == 9, all(colSums(is.na(mat)) == 0))

annot <- data.frame(gene_symbol = sub$gene_symbol, row.names = sub$gene_id)

write.csv(mat, "data/GSE307956_raw_count_matrix.csv")
write.csv(annot, "data/GSE307956_gene_annotation.csv")

cat("Wrote data/GSE307956_raw_count_matrix.csv (",
    nrow(mat), "genes x ", ncol(mat), " samples)\n", sep = "")
cat("Wrote data/GSE307956_gene_annotation.csv (", nrow(annot), " ids)\n", sep = "")
cat("Provenance: ", url, "\n", sep = "")