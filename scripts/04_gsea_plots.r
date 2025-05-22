#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(stringr)

# === Paths and output dir ===
gsea_dir    <- "../results/gsea/tables"
dotplot_dir <- file.path(dirname(gsea_dir), "plots")
dir.create(dotplot_dir, showWarnings = FALSE, recursive = TRUE)

# === Gather GSEA files ===
gsea_files <- list.files(gsea_dir, pattern = "\\.tsv$", full.names = TRUE)

# === Dot-plots (top 40 up & top 40 down; 16″×20″) ===
for (file in gsea_files) {
  # strip both the "GSEA_DE_" prefix AND the ".tsv" extension:
  contrast <- sub("^GSEA_DE_", "",
              tools::file_path_sans_ext(basename(file)))

  gsea <- fread(file)
  gsea <- gsea[!is.na(NES) & !is.na(padj) & is.finite(NES) & is.finite(padj)]

  # pick top 40 up and top 40 down
  top_up   <- gsea[padj < 0.05 & NES >  0][order(-NES)][1:40]
  top_down <- gsea[padj < 0.05 & NES <  0][order( NES)][1:40]
  tp       <- rbind(top_up, top_down)
  tp       <- tp[!is.na(pathway) & is.finite(NES) & is.finite(padj)]

  if (nrow(tp) < 2) {
    message("Not enough enriched pathways in ", basename(file),
            " (found ", nrow(tp), "), skipping dot plot.")
    next
  }

  p <- ggplot(tp, aes(
         x     = reorder(pathway, NES),
         y     = NES,
         color = NES > 0,
         size  = -log10(padj)
       )) +
       geom_point(alpha = 0.8) +
       scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
       coord_flip() +
       scale_color_manual(values = c("TRUE" = "firebrick",
                                     "FALSE"= "steelblue")) +
       scale_size_continuous(range = c(2, 6)) +
       theme_classic() +
       labs(
         title = paste("GSEA –", contrast),
         x     = "Pathway",
         y     = "Normalized Enrichment Score (NES)",
         color = "Upregulated",
         size  = "-log10(adj. p-value)"
       )

  out_file <- file.path(dotplot_dir, paste0("gsea_dotplot_", contrast, ".png"))
  tryCatch({
    ggsave(out_file, plot = p, width = 16, height = 20, units = "in", dpi = 150)
    message("Saved dotplot for ", contrast)
  }, error = function(e) {
    message("Failed to save dotplot for ", contrast, ": ", e$message)
  })
}
