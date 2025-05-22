library(fgsea)
library(data.table)
library(msigdbr)
library(org.Hs.eg.db)
library(tidyverse)

# === Load gene sets ===
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
msig_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
msig_c7 <- msigdbr(species = "Homo sapiens", category = "C7")

# Combine all gene sets
msig_all <- rbind(msig_h, msig_c2, msig_c7)
pathways <- split(msig_all$gene_symbol, msig_all$gs_name)

# === DE result files ===
de_files <- list.files("../results/deseq2", pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# Output directory
gsea_dir <- "../results/gsea"
dir.create(gsea_dir, showWarnings = FALSE)

# === Run GSEA on each file ===
for (file in de_files) {
  message("Processing ", basename(file))

  res <- fread(file)
  res <- res[!is.na(padj) & !is.na(log2FoldChange)]

  # Ensure gene symbols are present
  if (!"symbol" %in% colnames(res)) {
    stop("No 'symbol' column found in ", file)
  }

  # Remove duplicated genes (keep most significant)
  res <- res[!duplicated(symbol)]

  # Rank genes
  ranked <- res$log2FoldChange
  names(ranked) <- res$symbol
  ranked <- sort(ranked, decreasing = TRUE)

  # Filter to genes in MSigDB
  ranked <- ranked[names(ranked) %in% unique(unlist(pathways))]

  # Run fgsea
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranked,
    nperm = 10000,
    minSize = 15,
    maxSize = 500
  ) %>% as.data.table() %>%
    arrange(padj)

  # Save result
  out_file <- file.path(gsea_dir, paste0("GSEA_", tools::file_path_sans_ext(basename(file)), ".tsv"))
  fwrite(fgsea_res, out_file, sep = "\t")
}

message("All GSEA analyses complete.")