# PBMC ArchR Analysis Run Order

This folder contains the analysis scripts for PBMC samples.

## Data

Please follow the README.md in the root folder to download the data from zenodo.

## Run order

1. `archR_CoCnT_PBMCs_H3K4me3_Doublet_Removal_Final.Rmd`
2. `archR_CoCnT_PBMCs_H3K27me3_Doublet_Removal_Final.Rmd`
3. `archR_CoCnT_PBMCs_H3K4me3_H3K27me3_Export_Final.Rmd`
4. `CoCnT_PBMC_AnnData_Integration_Final.ipynb`
5. `archR_CoCnT_PBMCs_H3K4me3_H3K27me3_Cluster_analysis_Final.Rmd`

## Dependency notes

- The two `Doublet_Removal` Rmd files are the main starting points.
- They are independent of each other, so they can be run in either order.
- Both must finish before `archR_CoCnT_PBMCs_H3K4me3_H3K27me3_Export_Final.Rmd`, because the export notebook loads the final doublet-removed ArchR projects from both marks.
- `CoCnT_PBMC_AnnData_Integration_Final.ipynb` is downstream of the export step, and generate the WNN cluster assignment file which is needed for next step.

## Other files

- `CoCnT_PBMCs_PR001695_PR001743_Read_Counts_FRiPs_Final.ipynb` is separate QC/statistics work for read counts and FRiP. It is not required for the main ArchR -> export -> integration pipeline. And the raw data files are not included in the zenodo upload, please contact the author if you want to run this notebook.

