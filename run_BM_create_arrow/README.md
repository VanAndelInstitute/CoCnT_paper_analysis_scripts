# BM CoCnT Arrow File Creation

This directory contains two R scripts for generating ArchR Arrow files from bone marrow CoCnT BED files.

## Files

- `run_create_arrow_file_mooc.R`: Creates Arrow files from mono-occupancy BED files.
- `run_create_arrow_file_cooc.R`: Creates Arrow files from co-occupancy BED files after merging paired BED inputs within each library.

## Requirements

- R packages:
  - `ArchR`
  - `data.table`
  - `magrittr`
  - `readr` (`run_create_arrow_file_cooc.R` only)
  - `stringr` (`run_create_arrow_file_cooc.R` only)
- Reference genome configured in ArchR: `hg38`
- `bgzip` available on the command line for `run_create_arrow_file_cooc.R`

## Input Data

Both scripts expect the following relative inputs to exist:

- Sample sheets:
  - `../../data/sample_sheet_PR001798.tsv`
  - `../../data/sample_sheet_PR001799.tsv`
  - `../../data/sample_sheet_PR001855.tsv`
  - `../../data/sample_sheet_PR001856.tsv`
- BED directories:
  - `../../data/Janssens_Lab_Data/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
  - `../../data/Janssens_Lab_Data/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
  - `../../data/Janssens_Lab_Data/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
  - `../../data/Janssens_Lab_Data/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`

## `run_create_arrow_file_mooc.R`

This script creates Arrow files directly from per-sample BED files across the four batches.

### Behavior

- Builds BED paths by combining each `sample_name` with its batch-specific BED directory.
- Combines all sample sheets into one table.
- Generates one Arrow file per sample in the current working directory.

### ArchR settings

- Genome: `hg38`
- Threads: `20`
- `minTSS = 0`
- `minFrags = 100`
- `maxFrags = Inf`
- `nChunk = 1`
- `addTileMat = FALSE`
- `addGeneScoreMat = FALSE`
- `excludeChr = c("chrM")`
- `force = TRUE`

### Output

- Arrow files written to the current working directory, typically as `sample_name.arrow`

## `run_create_arrow_file_cooc.R`

This script creates Arrow files for co-occupancy data by first merging paired BED files belonging to the same `library_id`.

### ArchR settings

- Genome: `hg38`
- Threads: `20`
- `minTSS = 0`
- `minFrags = 0`
- `maxFrags = Inf`
- `nChunk = 1`
- `addTileMat = FALSE`
- `addGeneScoreMat = FALSE`
- `excludeChr = c("chrM")`
- `force = TRUE`

### Output

- Merged intermediate BED files in `./tmp/`
- Compressed merged BED files in `./tmp/*.bed.gz`
- Arrow files written to the current working directory

## How to run

From this directory:

```r
Rscript run_create_arrow_file_mooc.R
Rscript run_create_arrow_file_cooc.R
```

## Notes

- The mono-occupancy workflow filters out cells with fewer than 100 fragments.
- The co-occupancy workflow applies no fragment filtering (`minFrags = 0`).
- Neither script adds tile matrices or gene score matrices during Arrow creation.
- Both scripts exclude mitochondrial reads by removing `chrM`.
