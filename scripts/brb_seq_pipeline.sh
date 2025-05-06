#!/bin/bash
set -euo pipefail

# ==================== CONFIGURATION ====================
THREADS=8
FASTQ_DIR="../data/fastq"
MERGED_DIR="../workdir/merged"
QC_DIR="../results/qc"
GENOME_DIR="../genomes/human/STAR_index"
GENOME_FA="../genomes/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="../genomes/human/Homo_sapiens.GRCh38.108.gtf"
BARCODE_WHITELIST="../data/barcodes_96_V5A_star.txt"
BAM_DIR="../workdir/aligned"
COUNT_MATRIX_DIR="../results/matrix"
LOG_DIR="../logs"

mkdir -p "$MERGED_DIR" "$QC_DIR" "$GENOME_DIR" "$BAM_DIR" "$COUNT_MATRIX_DIR" "$LOG_DIR"

# ==================== STEP 1: Merge FASTQ Files ====================
echo "[Step 1] Merging FASTQ files"
cat $FASTQ_DIR/*_R1.fastq.gz > "$MERGED_DIR/library_R1.fastq.gz"
cat $FASTQ_DIR/*_R2.fastq.gz > "$MERGED_DIR/library_R2.fastq.gz"

# ==================== STEP 2: Run FastQC ====================
echo "[Step 2] Running FastQC"
fastqc --outdir "$QC_DIR" "$MERGED_DIR/library_R1.fastq.gz"
fastqc --outdir "$QC_DIR" "$MERGED_DIR/library_R2.fastq.gz"

# ==================== STEP 3: Generate STAR Genome Index (Only once) ====================
if [ ! -f "$GENOME_DIR/SA" ]; then
  echo "[Step 3] Generating STAR genome index"
  STAR --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF" \
    --runThreadN "$THREADS"
fi

# ==================== STEP 4: Run STARsolo Alignment ====================
echo "[Step 4] Running STARsolo alignment and count generation"
STAR \
  --runMode alignReads \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$MERGED_DIR/library_R2.fastq.gz" "$MERGED_DIR/library_R1.fastq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$BAM_DIR/" \
  --runThreadN "$THREADS" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMmapqUnique 60 \
  --outBAMsortingThreadN "$THREADS" \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 14 \
  --soloUMIstart 15 --soloUMIlen 14 \
  --soloCBwhitelist "$BARCODE_WHITELIST" \
  --soloStrand Forward \
  --soloCellFilter None \
  --soloUMIdedup NoDedup 1MM_All \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene \
  --quantMode GeneCounts \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --outFilterMultimapNmax 1 > "$LOG_DIR/star_solo.log"

# ==================== STEP 5: Convert MTX to Count Matrix ====================
echo "[Step 5] Converting .mtx to count matrix using R"
Rscript - <<EOF
library(data.table)
library(Matrix)
matrix_dir <- "$BAM_DIR/Solo.out/Gene/raw/"
mat <- as.data.frame(as.matrix(readMM(file.path(matrix_dir, "matrix.mtx"))))
features <- fread(file.path(matrix_dir, "features.tsv"), header = FALSE)
barcodes <- fread(file.path(matrix_dir, "barcodes.tsv"), header = FALSE)
colnames(mat) <- barcodes\$V1
rownames(mat) <- features\$V1
fwrite(mat, file = "$COUNT_MATRIX_DIR/umi.counts.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
EOF

# ==================== DONE ====================
echo "Pipeline complete. UMI count matrix saved at: $COUNT_MATRIX_DIR/umi.counts.txt"
echo "FastQC reports are available in $QC_DIR"
echo "STAR logs in $LOG_DIR"
echo "Aligned BAM and STARsolo outputs in $BAM_DIR"
