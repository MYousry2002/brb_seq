import pandas as pd

# Load files
wells_samples = pd.read_csv("../data/wells_samples.tsv", sep="\t")
barcodes = pd.read_csv("../data/barcodes_96_V5A.tsv", sep="\t")

# Merge on 'Name' column (well ID)
metadata = pd.merge(barcodes, wells_samples, on="Name", how="inner")

# Rename columns
metadata = metadata.rename(columns={"B1": "Barcode", "Sample": "Condition"})

# Remove samples with undefined conditions (underscore)
metadata = metadata[metadata["Condition"] != "_"]

# Add replicate number per condition
metadata["Replicate"] = metadata.groupby("Condition").cumcount() + 1

# Create SampleName: e.g., rapa1_1, rapa1_2, ...
metadata["SampleName"] = metadata["Condition"] + "_" + metadata["Replicate"].astype(str)

# Reorder and export
metadata[["Barcode", "Condition", "SampleName"]].to_csv("../data/sample_metadata.tsv", sep="\t", index=False)