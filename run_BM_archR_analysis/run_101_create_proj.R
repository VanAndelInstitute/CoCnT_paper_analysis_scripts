########################################################################3
# Script: run_n_create_proj.R
# Purpose: Create ArchR project for all samples, add metadata, and subset to CoCnT samples with H3K27me3, H3K4me1, H3K4me2, H3K4me3
# INPUT: 
#   - Arrow files for all samples
#   * ../run_BM_create_arrow/
#   - Sample sheets for all samples
#   * ../data/sample_sheet_PR001799.tsv
#   * ../data/sample_sheet_PR001798.tsv
#   * ../data/sample_sheet_PR001855.tsv
#   * ../data/sample_sheet_PR001856.tsv
#   - Souporcell results for all samples
#   * ../run_BM_souporcell/tmp/PR001798_PR001799_PR001855_PR001856_souporcell_output/clusters.tsv
# OUTPUT:
#  - ArchR project for all samples with metadata added
#  * ./tmp/ArchRProject_4
#  * ./tmp/cell_metadata_4.tsv
#  - Subsetted ArchR projects for CoCnT samples with H3K27me3, H3K4me1, H3K4me2, H3K4me3
#  * ./tmp/H3K27me3_coCnT_all
#  * ./tmp/H3K4me1_coCnT
#  * ./tmp/H3K4me2_coCnT
#  * ./tmp/H3K4me3_coCnT

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

## Load ArchR
library(ArchR)
addArchRGenome("hg38")

## Input
## Load sample sheet
sample_sheet <- rbind(
    fread("../data/sample_sheet_PR001799.tsv"),
    fread("../data/sample_sheet_PR001798.tsv"),
    fread("../data/sample_sheet_PR001855.tsv"),
    fread("../data/sample_sheet_PR001856.tsv")
)
## Load souporcell results
d_soupor_cell = fread("../run_BM_souporcell/tmp/PR001798_PR001799_PR001855_PR001856_souporcell_output/clusters.tsv")

## Check if any sample is missing
arrow_file <- dir("../run_BM_create_arrow/", pattern = "arrow$", full = T)
sample_name_x <- str_replace(basename(arrow_file), ".arrow", "")
sample_sheet[!(sample_name %in% sample_name_x)]

#+ eval=FALSE
## create an ArchRProject
proj <- ArchRProject(
    ArrowFiles = arrow_file,
    outputDirectory = "./tmp/ArchRProject_4",
    copyArrows = TRUE
)

## Add metadata
old_meta <- proj@cellColData
old_meta$cell_id_ar = rownames(old_meta)
sample_sheet$library_type %>% table
new_meta <- merge(old_meta, sample_sheet, by.x = "Sample", by.y = "sample_name")
new_meta$cell_id <- sub(".*#", "", new_meta$cell_id_ar)
new_meta = merge(new_meta, d_soupor_cell, by.x = "cell_id", by.y = "barcode", all.x = TRUE)
rownames(new_meta) <- new_meta$cell_id_ar

proj@cellColData <- new_meta
dim(new_meta)

saveArchRProject(proj, outputDirectory = "./tmp/ArchRProject_4", load = FALSE, overwrite = TRUE)

write_tsv(data.frame(new_meta), "./tmp/cell_metadata_4.tsv")

#+ eval=TRUE
proj = loadArchRProject("./tmp/ArchRProject_4", force = TRUE)
## Create subset

d_meta <- data.table(data.frame(proj@cellColData))
d_meta$cell_id_ar <- rownames(proj@cellColData)

d_meta[, .N, by = .(antibody_target, library_type)]

d_meta[ antibody_target == "H3K4me1" & library_type == "CoCnT" & nFrags > 250, cell_id_ar] %>% length()
d_meta[ antibody_target == "H3K4me2" & library_type == "CoCnT" & nFrags > 250, cell_id_ar] %>% length()
d_meta[ antibody_target == "H3K4me3" & library_type == "CoCnT" & nFrags > 250, cell_id_ar] %>% length()
d_meta[ antibody_target == "H3K27me3" & library_type == "CoCnT" & nFrags > 250, cell_id_ar] %>% length()

d_meta[, .(N = .N), by = .(batch, library_id)]


## Select the H3K27me3
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_target == "H3K27me3" & library_type == "CoCnT" & !is.na(status), cell_id_ar],
    outputDirectory = "./tmp/H3K27me3_coCnT_all",
    dropCells = TRUE,
    force = TRUE
)

## Select the H3K4me3
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_target == "H3K4me3" & library_type == "CoCnT" & !is.na(status), cell_id_ar],
    outputDirectory = "./tmp/H3K4me3_coCnT",
    dropCells = TRUE,
    force = TRUE
)

## Select the H3K4me2
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_target == "H3K4me2" & library_type == "CoCnT" & !is.na(status), cell_id_ar],
    outputDirectory = "./tmp/H3K4me2_coCnT",
    dropCells = TRUE,
    force = TRUE
)

## Select the H3K4me1
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_target == "H3K4me1" & library_type == "CoCnT" & !is.na(status), cell_id_ar],
    outputDirectory = "./tmp/H3K4me1_coCnT",
    dropCells = TRUE,
    force = TRUE
)
