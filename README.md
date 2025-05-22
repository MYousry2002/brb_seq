# BRB-seq Pipeline

This repository contains scripts and resources for processing BRB-seq (Bulk RNA Barcoding and Sequencing) data to perform:

- **Alignment** and count matrix generation  
- **Differential gene expression (DGE)** analysis using DESeq2  
- **Gene Set Enrichment Analysis (GSEA)** using MSigDB hallmark and curated gene sets  
- **Visualization** including PCA plots, volcano plots, and dotplots


## 📁 Directory Structure

```
BRB_SEQ/
├── data/
│   ├── fastq/                 # Raw FASTQ files (R1, R2 for each lane)
│   ├── barcodes_96_V5A.tsv    # whitelist 
│   ├── sample_metadata.tsv    # metadata
│   └── wells_samples.tsv      # plate layout
├── logs/                      # STARsolo and pipeline logs
│   ├── brbseq_pipeline.err
│   ├── brbseq_pipeline.out
│   └── star_solo.log
├── results/
│   ├── matrix/                # Output count matrix (UMI counts)
│   ├── deseq2/                # Differential expression result tables
│   ├── volcano_plots/         # Volcano plots from DE results
│   ├── gsea/                  # GSEA results + plots
│   ├── pca/                   # PCA plots
│   ├── qc/                    # (optional) QC stats or reports
│   └── validation/            # Additional summary plots (PCA, etc.)
├── scripts/                   # All main scripts in order of execution
│   ├── 00_brb_seq_pipeline.sh         # Alignment & quantification using STARsolo
│   ├── 01_run_deseq2.r                # DESeq2-based differential analysis
│   ├── 02_make_volcano_plots.r        # Volcano plots of DEGs
│   ├── 03_gsea.r                      # GSEA with MSigDB (Hallmark, C2, C7)
│   ├── 04_gsea_plots.r               # Dotplot + heatmap of GSEA output
│   ├── generate_sample_metadata.py   # Parses plate layout into metadata
│   └── run_pipeline_job.sh           # Job submission wrapper for SCC
├── workdir/
│   ├── aligned/              # STARsolo output
│   └── merged/               # Merged outputs if applicable
├── env.yml                   # Conda environment definition
├── env.sh                    # Environment/module loading
└── .gitignore
```

---

## Workflow

0. **Run alignment & quantification**  
   Submit job via:
   ```
   qsub scripts/run_pipeline_job.sh
   ```

1. **Run DESeq2 analysis**
   ```
   Rscript scripts/01_run_deseq2.r
   ```

2. **Create volcano plots**
   ```
   Rscript scripts/02_make_volcano_plots.r
   ```

3. **Run GSEA**
   ```
   Rscript scripts/03_gsea.r
   ```

4. **Plot GSEA dotplots**
   ```
   Rscript scripts/04_gsea_plots.r
   ```

---

## Conda Environment

Install the pipeline dependencies via:

```bash
conda env create -f env.yml
conda activate brb_seq
```

---

## Author
Mohamed Yousry ElSadec

PhD Student in Bioinformatics, Boston University

Juan Fuxman Bass Lab

---

## License

MIT License
