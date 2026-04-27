# BM CoCnT Arrow File Creation

This directory contains R scripts for generating ArchR Arrow files from T cell depleted bone marrow CoCnT BED files.

## Files

- `run_create_arrow_file_mooc.R`: Creates Arrow files from mono-occupancy BED files.
- `run_create_arrow_file_cooc.R`: Creates Arrow files from co-occupancy BED files after merging paired BED inputs within each library.

## Input Data

Both scripts expect the following relative inputs to exist.

- Sample sheets:
  - `../data/sample_sheet_PR001798.tsv`
  - `../data/sample_sheet_PR001799.tsv`
  - `../data/sample_sheet_PR001855.tsv`
  - `../data/sample_sheet_PR001856.tsv`
- BED directories:
  - `../data/Janssens_Lab_Data/BMMC_Beds/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
  - `../data/Janssens_Lab_Data/BMMC_Beds/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
  - `../data/Janssens_Lab_Data/BMMC_Beds/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
  - `../data/Janssens_Lab_Data/BMMC_Beds/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`

You can download the input from Zenodo by following the instructions in the main README file.

## `run_create_arrow_file_mooc.R`

This script creates Arrow files directly from per-sample BED files across the four batches.

### Output

- Arrow files written to the current working directory, typically as `sample_name.arrow`

## `run_create_arrow_file_cooc.R`

This script creates Arrow files for co-occupancy data. 

## How to run

From this directory:

```r
Rscript run_create_arrow_file_mooc.R
Rscript run_create_arrow_file_cooc.R
```

