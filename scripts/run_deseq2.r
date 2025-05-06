# === Load libraries ===
library(DESeq2)
library(data.table)

# === Load data ===
counts <- fread("umi.counts.txt")
metadata <- fread("sample_metadata.tsv")

# Set rownames (gene IDs)
genes <- counts[[1]]
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- genes

# Ensure sample order matches
metadata <- metadata[metadata$Barcode %in% colnames(count_matrix)]
count_matrix <- count_matrix[, metadata$Barcode]

# === Create DESeq2 object ===
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Condition
)

dds <- DESeq(dds)

# === Differential Expression Analysis ===
# Control mapping
control_map <- list(
  rapa = "dmso",
  tuni = "dmso",
  mg = "dmso",
  oligo = "dmso",
  hs = "37"
)

# Output directory
results_dir <- "../results/deseq2"
dir.create(results_dir, showWarnings = FALSE)

# Run each comparison: drug-dose vs control
all_conditions <- unique(metadata$Condition)

for (cond in all_conditions) {
  for (prefix in names(control_map)) {
    if (startsWith(cond, prefix)) {
      control <- control_map[[prefix]]
      if (control %in% all_conditions) {
        message("Running DE: ", cond, " vs ", control)
        res <- results(dds, contrast = c("Condition", cond, control))
        res_df <- as.data.frame(res[order(res$pvalue), ])
        out_file <- file.path(results_dir, paste0("DE_", cond, "_vs_", control, ".tsv"))
        fwrite(res_df, out_file, sep = "\t")
      }
    }
  }
}

# === Electroporation effect ===
if ("dmso" %in% all_conditions && "ne" %in% all_conditions) {
  message("Running DE: dmso vs ne")
  res <- results(dds, contrast = c("Condition", "dmso", "ne"))
  res_df <- as.data.frame(res[order(res$pvalue), ])
  fwrite(res_df, file.path(results_dir, "DE_dmso_vs_ne.tsv"), sep = "\t")
}

if ("37" %in% all_conditions && "ne" %in% all_conditions) {
  message("Running DE: 37 vs ne")
  res <- results(dds, contrast = c("Condition", "37", "ne"))
  res_df <- as.data.frame(res[order(res$pvalue), ])
  fwrite(res_df, file.path(results_dir, "DE_37_vs_ne.tsv"), sep = "\t")
}

message("All DE comparisons completed.")