library(data.table)
library(ggplot2)
library(pheatmap)
library(stringr)

# === Load all GSEA result files ===
gsea_dir <- "../results/gsea"
gsea_files <- list.files(gsea_dir, pattern = "\\.tsv$", full.names = TRUE)

# === Output directories ===
dotplot_dir <- "../results/gsea/plots/dotplots"
dir.create(dotplot_dir, showWarnings = FALSE)

# Store NES values for heatmap
nes_matrix <- list()

# === Generate Dot Plot per Comparison ===
for (file in gsea_files) {
  contrast <- str_replace(basename(file), "GSEA_DE_|\\.tsv", "")
  gsea <- fread(file)

  # Clean
  gsea <- gsea[!is.na(NES) & !is.na(padj)]

  # Top 10 up and downregulated pathways (padj < 0.05)
  top_up <- gsea[padj < 0.05 & NES > 0][order(-NES)][1:10]
  top_down <- gsea[padj < 0.05 & NES < 0][order(NES)][1:10]
  top_pathways <- rbind(top_up, top_down)

  # Store for heatmap
  nes_matrix[[contrast]] <- setNames(gsea$NES, gsea$pathway)

  # Dot plot with padj as size
  p <- ggplot(top_pathways, aes(x = reorder(pathway, NES),
                                y = NES,
                                color = NES > 0,
                                size = -log10(padj))) +
    geom_point(alpha = 0.8) +
    coord_flip() +
    scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
    scale_size_continuous(range = c(2, 6)) +
    theme_minimal() +
    labs(title = paste("GSEA –", contrast),
         x = "Pathway",
         y = "Normalized Enrichment Score (NES)",
         color = "Upregulated",
         size = "-log10(adj. p-value)")

  ggsave(file.path(dotplot_dir, paste0("gsea_dotplot_", contrast, ".png")),
         plot = p, width = 8, height = 6)
}


# === Build NES matrix for heatmap ===
all_pathways <- unique(unlist(lapply(nes_matrix, names)))
all_contrasts <- names(nes_matrix)

nes_df <- data.frame(pathway = all_pathways)
for (contrast in all_contrasts) {
  nes_vector <- nes_matrix[[contrast]]
  nes_df[[contrast]] <- nes_vector[match(all_pathways, names(nes_vector))]
}

# Convert to matrix
rownames(nes_df) <- nes_df$pathway
nes_mat <- as.matrix(nes_df[, -1])
nes_mat[is.na(nes_mat)] <- 0  # Fill missing NES with 0

# filter pathways with at least one significant enrichment
keep <- rowSums(abs(nes_mat) > 1) > 0
nes_mat <- nes_mat[keep, ]

# === Save heatmap to file ===
heatmap_dir <- "../results/gsea/plots"
dir.create(heatmap_dir, showWarnings = FALSE)

# Save as PNG
png(file.path(heatmap_dir, "gsea_heatmap.png"), width = 1200, height = 1000, res = 150)
pheatmap(nes_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
         border_color = NA,
         main = "GSEA – NES Across Conditions")
dev.off()
