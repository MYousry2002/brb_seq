# === Load libraries ===
library(DESeq2)
library(data.table)
library(ggplot2)
library(ggrepel)

# === Load data ===
counts <- fread("../results/matrix/umi.counts.txt")
metadata <- fread("../data/sample_metadata.tsv")

# === Set gene IDs as rownames and extract count matrix ===
genes <- counts[[1]]
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- genes

# === Keep only barcodes that are in metadata ===
valid_barcodes <- colnames(count_matrix)[colnames(count_matrix) %in% metadata$Barcode]
count_matrix <- count_matrix[, valid_barcodes, drop = FALSE]

# === Subset and reorder metadata to match count matrix columns ===
metadata <- metadata[metadata$Barcode %in% valid_barcodes, ]
metadata <- metadata[match(valid_barcodes, metadata$Barcode), ]
rownames(metadata) <- metadata$Barcode

# === Create DESeq2 object ===
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Condition
)

# === Run DESeq2 normalization and modeling ===
dds <- DESeq(dds, fitType="local")

# === PCA Plot (for QC) ===
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Add SampleName to PCA data using Barcode match
pca_data$SampleName <- metadata$SampleName[match(pca_data$name, metadata$Barcode)]

# Save PCA plot
dir.create("../results/validation", showWarnings = FALSE)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = SampleName)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of VST-transformed counts") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../results/validation/pca_plot.png", width = 6, height = 5)


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
        res_df$gene <- rownames(res_df)  # <- Add gene names
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
  res_df$gene <- rownames(res_df)
  fwrite(res_df, file.path(results_dir, "DE_dmso_vs_ne.tsv"), sep = "\t")
}

if ("37" %in% all_conditions && "ne" %in% all_conditions) {
  message("Running DE: 37 vs ne")
  res <- results(dds, contrast = c("Condition", "37", "ne"))
  res_df <- as.data.frame(res[order(res$pvalue), ])
  res_df$gene <- rownames(res_df)
  fwrite(res_df, file.path(results_dir, "DE_37_vs_ne.tsv"), sep = "\t")
}

message("All DE comparisons completed.")