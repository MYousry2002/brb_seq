library(data.table)
library(ggplot2)

# List all DE result files
files <- list.files("../results/deseq2", pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# Output folder
volcano_dir <- "../results/volcano_plots"
dir.create(volcano_dir, showWarnings = FALSE)

for (file in files) {
  res <- fread(file)
  res$gene <- rownames(res)
  
  res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "yes", "no")

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    theme_minimal() +
    labs(title = basename(file), x = "log2 Fold Change", y = "-log10(p-value)") +
    theme(legend.position = "none")

  ggsave(filename = file.path(volcano_dir, paste0(basename(file), ".png")), plot = p, width = 6, height = 5)
}