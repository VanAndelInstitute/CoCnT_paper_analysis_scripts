# run_BM_bed_profile

This folder contains the bone marrow CoCnT fragment QC and pseudo-bulk BED/BigWig generation workflow.

## Overview

The workflow is split into two R scripts:

- `run_qc.R`: loads fragment BED files, builds fragment-per-cell summaries, merges them with sample metadata, and produces QC plots.
- `run_pseudo_bulk_bed.R`: filters fragments to annotated cells, writes pseudo-bulk BED files by cluster and histone mark combination, and converts those BED files into bedGraph and BigWig tracks.

## Inputs

Common metadata inputs:

- `../data/sample_sheet_PR001798.tsv`
- `../data/sample_sheet_PR001799.tsv`
- `../data/sample_sheet_PR001855.tsv`
- `../data/sample_sheet_PR001856.tsv`

QC-specific raw inputs:

- BED files under `../data/Janssens_Lab_Data/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
- BED files under `../data/Janssens_Lab_Data/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
- BED files under `../data/Janssens_Lab_Data/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
- BED files under `../data/Janssens_Lab_Data/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`

Pseudo-bulk-specific inputs:

- `../run_BM_archR_analysis/tmp/table_cell_cluster_annotation_final_round2.tsv`
- `../run_BM_archR_analysis/tmp/table_cellColData_single_all.tsv`
- `./tmp/all_fragments.rds` or `./tmp/table_good_cells_fragments.tsv`
- `../../data/ref_genome/hg38.sizes`

## Outputs

Files written by the workflow include:

- `./tmp/all_fragments.rds`
- `../tmp/cell_fragments.tsv`
- `./tmp/table_good_cells_fragments.tsv`
- `./tmp/pseudo_bulk_bed/*_pseudo_bulk.bed`
- `./tmp/pseudo_bulk_bed/*.sorted.bed`
- `./tmp/pseudo_bulk_bed/*.bedGraph`
- `./tmp/pseudo_bulk_bed/*.bw`
- `./tmp/pseudo_bulk_bed/*.scale_co.bedGraph`
- `./tmp/pseudo_bulk_bed/*.scale_co.bw`
- `./tmp/pseudo_bulk_bed/*.scale_cell_co.bedGraph`
- `./tmp/pseudo_bulk_bed/*.scale_cell_co.bw`

`run_qc.R` also renders QC plots to the active graphics device via `httpgd`.

## Script Notes

### `run_qc.R`

This script:

- combines sample-sheet metadata across four batches
- reads raw fragment BED files and caches them as `all_fragments.rds`
- computes fragment counts per cell
- merges counts with sample metadata
- generates bar plots and fragment-per-cell histograms across library type, donor, antibody target, and antibody species

The script includes cached sections marked with `#+ eval=F` and `#+ eval=T`, so it can be run either from raw BED inputs or from saved intermediate tables.

### `run_pseudo_bulk_bed.R`

This script:

- loads cluster annotations and filtered fragment tables
- creates pseudo-bulk BED files for each cluster and antibody target
- creates pseudo-bulk BED files for antibody combinations in `CoCnt1` and `CoCnt2`
- converts BED files to sorted BED, bedGraph, and BigWig tracks using UCSC tools
- generates additional normalized coverage tracks using coverage-based and cell-count-based scaling

The script expects `bedtools` and `bedGraphToBigWig` to be available through the `ucsc` conda environment.

## Typical Order

1. Run `run_qc.R` to inspect fragment distributions and generate cached fragment summaries.
2. Prepare or confirm cluster annotations in `../run_BM_archR_analysis/tmp/`.
3. Run `run_pseudo_bulk_bed.R` to generate pseudo-bulk BED and BigWig tracks.

## Working Assumptions

- The scripts are written as analysis scripts rather than reusable functions.
- Relative paths assume execution from within `run_BM_bed_profile/`.
- The `tmp/` directory stores intermediate and derived files for this workflow.
