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

- BED files under `../data/Janssens_Lab_Data/BMMC_Beds/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
- BED files under `../data/Janssens_Lab_Data/BMMC_Beds/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
- BED files under `../data/Janssens_Lab_Data/BMMC_Beds/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`
- BED files under `../data/Janssens_Lab_Data/BMMC_Beds/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/`

Pseudo-bulk-specific inputs:

- `../run_BM_archR_analysis/tmp/table_cell_cluster_annotation_final_round2.tsv`
    This file should be prepared in advance by running the ArchR analysis workflow in `../run_BM_archR_analysis/run_117_joint_umpa.R`.
- `./tmp/all_fragments.rds` 
- `../data/ref_genome/hg38.sizes`, please prepare this file in advance by downloading the reference genome FASTA and running `samtools faidx` on it.

## Outputs

Files written by the workflow include:

- `./tmp/all_fragments.rds`
- `./tmp/cell_fragments.tsv`
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

## Typical Order

1. Run `run_qc.R`
3. Run `run_pseudo_bulk_bed.R` to generate pseudo-bulk BED and BigWig tracks.

