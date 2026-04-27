
# BM Souporcell Pipeline

This directory contains the scripts used to prepare sciCoCnT BAM files from four bone marrow samples and run [Souporcell](https://github.com/wheaton5/souporcell) for doublet detection and donor assignment.

Samples included in this run:

- `PR001798`
- `PR001799`
- `PR001855`
- `PR001856`

## Files

- `run_prepare_bam_PR001798.sh`
- `run_prepare_bam_PR001799.sh`
- `run_prepare_bam_PR001855.sh`
- `run_prepare_bam_PR001856.sh`
  Preprocess per-sample BAM files: remove duplicates, add read-group style tags, merge, sort, and index.
- `run_prepare_barcode.R`
  Exports the cell barcode whitelist used by Souporcell.
- `run_souporcell_PR001798_PR001799_PR001855_PR001856.sh`
  Merges the per-sample BAMs and runs Souporcell with `k=4`.

## Input data

Per-sample BAM inputs are expected under, which is not included in this repo:

- `../data/Janssens_Lab_Data/BMMC_Bams/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BAM`
- `../data/Janssens_Lab_Data/BMMC_Bams/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BAM`
- `../data/Janssens_Lab_Data/BMMC_Bams/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BAM`
- `../data/Janssens_Lab_Data/BMMC_Bams/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BAM`

Additional required inputs:

- Sample sheets:
  `../data/sample_sheet_PR001798.tsv`,
  `../data/sample_sheet_PR001799.tsv`,
  `../data/sample_sheet_PR001855.tsv`,
  `../data/sample_sheet_PR001856.tsv`
- Arrow files from `../run_BM_coCnT_create_arrow/`, which can be downloaded from Zenodo (TODO: add link)
- Reference FASTA: `../data/ref_genome/hg38_gencode.fa`, which is not included in this repo.
- Tagging helper: `../lib/bin/bin_add_tags.py`

## Output

The only relevant output of Souporcell is attached, and can be found in `./tmp/PR001798_PR001799_PR001855_PR001856_souporcell/clusters.tsv`. 

