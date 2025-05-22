#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(pheatmap)
library(stringr)

# === Paths and output dirs ===
gsea_dir    <- "../results/gsea"
dotplot_dir <- file.path(gsea_dir, "plots", "dotplots")
heatmap_dir <- file.path(gsea_dir, "plots")

dir.create(dotplot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# === Gather GSEA files ===
gsea_files <- list.files(gsea_dir, pattern = "\\.tsv$", full.names = TRUE)
nes_matrix <- list()

# === Dot‐plots (point only, wider & taller) ===
for (file in gsea_files) {
  contrast <- str_replace(basename(file), "GSEA_DE_|\\.tsv", "")
  gsea     <- fread(file)
  gsea     <- gsea[!is.na(NES) & !is.na(padj) & is.finite(NES) & is.finite(padj)]

  # record NES for heatmap
  nes_matrix[[contrast]] <- setNames(gsea$NES, gsea$pathway)

  # select top 30 up and top 30 down
  top_up   <- gsea[padj < 0.05 & NES >  0][order(-NES)][1:30]
  top_down <- gsea[padj < 0.05 & NES <  0][order( NES)][1:30]
  tp       <- rbind(top_up, top_down)
  tp       <- tp[!is.na(pathway) & is.finite(NES) & is.finite(padj)]

  if (nrow(tp) < 2) {
    message("Not enough enriched pathways in ", basename(file),
            " (found ", nrow(tp), "), skipping dot plot.")
    next
  }

  # build point‐only plot
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

  tryCatch(
    {
      # make it larger: 16" wide x 20" tall
      ggsave(out_file, plot = p, width = 16, height = 20)
      message("✔ Saved dotplot for ", contrast)
    },
    error = function(e) {
      message("✗ Failed to save dotplot for ", contrast, ": ", e$message)
    }
  )
}

# === Build NES matrix for heatmap ===
all_pathways  <- unique(unlist(lapply(nes_matrix, names)))
all_contrasts <- names(nes_matrix)

nes_df <- data.frame(pathway = all_pathways, stringsAsFactors = FALSE)
for (c in all_contrasts) {
  vec        <- nes_matrix[[c]]
  nes_df[[c]] <- vec[match(all_pathways, names(vec))]
}
rownames(nes_df) <- nes_df$pathway
nes_mat         <- as.matrix(nes_df[, -1, drop = FALSE])
nes_mat[is.na(nes_mat)] <- 0

# filter to |NES| > 1 in any condition
keep    <- rowSums(abs(nes_mat) > 1) > 0
nes_mat <- nes_mat[keep, , drop = FALSE]

# === Heatmap ===
if (nrow(nes_mat) > 0) {
  png(file.path(heatmap_dir, "gsea_heatmap.png"),
      width = 1200, height = 1000, res = 150)
  pheatmap(nes_mat,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color        = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
           border_color = NA,
           fontsize_row = 7,
           fontsize_col = 9,
           main         = "GSEA – NES Across Conditions")
  dev.off()
  message("✔ Saved GSEA heatmap")
} else {
  message("No pathways with |NES| > 1 across any contrast; skipping heatmap.")
}