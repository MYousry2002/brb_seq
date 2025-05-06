import pandas as pd

# Load files
wells_samples = pd.read_csv("../data/wells_samples.tsv", sep="\t")
barcodes = pd.read_csv("../data/barcodes_96_V5A.tsv", sep="\t")

# Merge on 'Name' column (well ID)
metadata = pd.merge(barcodes, wells_samples, on="Name", how="inner")

# Rename columns for clarity
metadata = metadata.rename(columns={"B1": "Barcode", "Sample": "Condition"})

# Reorder and export
metadata[["Barcode", "Condition"]].to_csv("../data/sample_metadata.tsv", sep="\t", index=False)