# CoCnT Analysis Scripts

This repository collects the analysis scripts used to process and analyze the PBMC CoCnT dataset and T cell depleted bone marrow CoCnT dataset. 

The code is organized by analysis workflow directory: each workflow directory contains scripts for one part of preprocessing, donor deconvolution, ArchR analysis, or pseudo-bulk track generation.

Detailed instructions for each stage are documented in the `README.md` file within the corresponding subdirectory.

## What's not included in this repo:

About sequencing data demultiplexing and mapping, please refering to our Nextflow workflow: https://github.com/vari-bbc/sciCT_pipeline

It needs `BAM` files to use SouporCells to identify the cell doublets. 
But the `BAM` files are very big for shipping with this repo or in Zenodo. 
So the souporcell step will only provide the script and the final output. 
If you are using scCoCnT and want to run SouporCells, please check out our Nextflow workflow here: https://github.com/vari-bbc/souper-star

## Repository Layout

- [`run_PBMC_coCnT_analysis/`](run_PBMC_coCnT_analysis/)
    Analysis scripts for the PBMC CoCnT dataset.
- [`run_BM_create_arrow/`](run_BM_create_arrow/)
    Generate ArchR Arrow files from mono-occupancy and co-occupancy BED inputs.
- [`run_BM_souporcell/`](run_BM_souporcell/)
    Run Souporcell for donor deconvolution and doublet identification.
- [`run_BM_archR_analysis/`](run_BM_archR_analysis/)
    Perform the core ArchR-based analysis.
- [`run_BM_bed_profile/`](run_BM_bed_profile/)
    Generate fragment QC summaries and pseudo-bulk BED/BigWig tracks.
- [`run_BM_motif`](run_BM_motif/)
    Perform motif enrichment analysis and TF footprinting.
- [`lib`](lib)
    Helper functions sourced by multiple analysis scripts.
- [`data/`](./data)
    The fold containing the raw input, it is not included in this repo.

## Prepare the `data` folder

Please download the `data` folder from Zenodo, which contains the raw Bed files and the sample sheet for the analysis.

Also for the reference scRNASeq data, please download from here: https://figshare.com/articles/dataset/Expression_of_97_surface_markers_and_RNA_transcriptome_wide_in_13165_cells_from_a_healthy_young_bone_marrow_donor/13397987?file=41038073

And place the file under the `data` folder.

Also for the reference genome files, please download the `hg38` reference genome files and the genome size files, and place them under the `data/ref_genome` folder.

They should named as `hg38.sizes`, `hg38_gencode.fa`.

