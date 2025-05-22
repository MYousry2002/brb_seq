#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(biomaRt)
library(ggrepel)

# === Setup output dirs ===
results_dir <- "../results/deseq2"
volcano_dir <- "../results/volcano_plots"
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)

# === List DE result files ===
files <- list.files(results_dir, pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# === Get all unique gene IDs for mapping ===
all_genes <- unique(unlist(lapply(files, function(f) fread(f)$gene)))

# === Map Ensembl IDs → HGNC symbols ===
ensembl  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map  <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = all_genes,
  mart       = ensembl
)
setDT(gene_map)
setnames(gene_map, c("ensembl_gene_id", "hgnc_symbol"), c("gene", "symbol"))

# === Process each DE file and make volcano plot ===
for (file in files) {
  # strip off directory + .tsv
  base_name <- tools::file_path_sans_ext(basename(file))

  # read DE results
  res <- fread(file)

  # preserve original order
  res[, idx__ := .I]

  # merge in gene symbols
  res <- merge(res, gene_map, by = "gene", all.x = TRUE, sort = FALSE)
  res[, symbol := fifelse(is.na(symbol) | symbol == "", gene, symbol)]
  setorder(res, idx__); res[, idx__ := NULL]

  # overwrite TSV with new symbol column (optional)
  orig_cols <- setdiff(names(res), "symbol")
  setcolorder(res, c(orig_cols, "symbol"))
  fwrite(res, file, sep = "\t")

  # define status for volcano
  res[, status := fifelse(padj < 0.05 & log2FoldChange >  1, "up",
                   fifelse(padj < 0.05 & log2FoldChange < -1, "down", "ns"))]

  # pick top 40 by smallest padj for labeling
  top_genes <- res[!is.na(padj)][order(padj)][1:40]

  # build volcano
  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = status)) +
       geom_point(alpha = 0.5, size = 1.5) +
       scale_color_manual(values = c(up = "red", down = "blue", ns = "grey")) +
       geom_text_repel(
         data        = top_genes,
         aes(label    = symbol),
         size        = 2.5,
         max.overlaps= 50
       ) +
       theme_classic() +
       labs(
         title = base_name,
         x     = "log2 Fold Change",
         y     = "-log10(p-value)"
       ) +
       theme(legend.position = "none")

  # save without the .tsv in the name
  out_png <- file.path(volcano_dir, paste0(base_name, ".png"))
  ggsave(out_png, plot = p, width = 6, height = 5, units = "in", dpi = 150)
  message("Saved volcano plot: ", out_png)
}