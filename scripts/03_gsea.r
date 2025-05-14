library(fgsea)
library(data.table)
library(msigdbr)
library(org.Hs.eg.db)
library(tidyverse)

# === Load Hallmark gene sets (MSigDB Category H) ===
msigdb <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msigdb$gene_symbol, msigdb$gs_name)

# === DE result files ===
de_files <- list.files("../results/deseq2", pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# Output directory
gsea_dir <- "../results/gsea_msigdb"
dir.create(gsea_dir, showWarnings = FALSE)

# === Run GSEA on each file ===
for (file in de_files) {
  message("Processing ", basename(file))

  res <- fread(file)
  res <- res[!is.na(padj) & !is.na(log2FoldChange)]
  res <- res[!duplicated(res$V1)]  # In case of duplicated genes
  res$symbol <- res$V1

  # Create ranked gene list
  ranked <- res$log2FoldChange
  names(ranked) <- res$symbol
  ranked <- sort(ranked, decreasing = TRUE)

  # Run fgsea
  fgsea_res <- fgsea(pathways = pathways,
                     stats = ranked,
                     nperm = 10000,
                     minSize = 15,
                     maxSize = 500)

  # Save result
  out_file <- file.path(gsea_dir, paste0("GSEA_", tools::file_path_sans_ext(basename(file)), ".tsv"))
  fwrite(fgsea_res, out_file, sep = "\t")
}

message("All GSEA analyses complete.")