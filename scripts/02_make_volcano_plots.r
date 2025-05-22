#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(biomaRt)
library(ggrepel)

# === Setup output dirs ===
results_dir <- "../results/deseq2"
volcano_dir <- "../results/volcano_plots"
dir.create(volcano_dir, showWarnings = FALSE)

# === List DE result files ===
files <- list.files(results_dir, pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# === Get all unique gene IDs ===
all_genes <- unique(unlist(lapply(files, function(f) fread(f)$gene)))

# === Map Ensembl IDs to gene symbols ===
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = ensembl
)
setDT(gene_map)
setnames(gene_map, c("ensembl_gene_id", "hgnc_symbol"), c("gene", "symbol"))

# === Process each file ===
for (file in files) {
  res <- fread(file)

  # Preserve original order by adding an index
  res[, idx__ := .I]

  # Add symbol column without reordering
  res <- merge(res, gene_map, by = "gene", all.x = TRUE, sort = FALSE)
  res[, symbol := ifelse(is.na(symbol) | symbol == "", gene, symbol)]

  # Reorder to original by index and remove helper column
  setorder(res, idx__)
  res[, idx__ := NULL]

  # Write updated TSV with symbol at the end
  original_cols <- setdiff(names(res), c("symbol"))
  setcolorder(res, c(original_cols, "symbol"))
  fwrite(res, file, sep = "\t")

  # Volcano plot
  res[, status := fifelse(padj < 0.05 & log2FoldChange > 1, "up",
                   fifelse(padj < 0.05 & log2FoldChange < -1, "down", "ns"))]

  top_genes <- res[!is.na(padj)][order(padj)][1:40]

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = status)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
    geom_text_repel(data = top_genes, aes(label = symbol), size = 2.5, max.overlaps = 50) +
    theme_classic() +
    labs(
      title = basename(file),
      x = "log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    theme(legend.position = "none")

  ggsave(file.path(volcano_dir, paste0(basename(file), ".png")), plot = p, width = 6, height = 5)
}