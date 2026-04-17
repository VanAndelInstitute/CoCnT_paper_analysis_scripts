##################################################
## Script: run_n_create_proj_cooc.R
## Purpose: Create ArchR project for co-occupancy reads
## INPUT:
## - Arrow files for all samples with co-occupancy reads
## * ../run_BM_coCnT_create_arrow/
## - Sample sheets for all samples
## * ../data/sample_sheet_PR001799.tsv
## * ../data/sample_sheet_PR001798.tsv
## * ../data/sample_sheet_PR001855.tsv
## * ../data/sample_sheet_PR001856.tsv
## - Cell cluster annotation for all cells
## * ./tmp/table_cell_cluster_annotation_final_round2.tsv
## OUTPUT:
## - ArchR project for co-occupancy reads with metadata added
## * ./tmp/ArchRProject_cooc_final
## * ./tmp/cell_metadata_cooc_4.tsv
## - Subsetted ArchR projects for co-occupancy reads with H3K27me3, H3K4me1, H3K4me2, H3K4me3
## * ./tmp/H3K27me3_H3K4me1_coCnT_all
## * ./tmp/H3K27me3_H3K4me2_coCnT_all
## * ./tmp/H3K27me3_H3K4me3_coCnT_all
## 
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

d_ann = fread("./tmp/table_cell_cluster_annotation_final_round2.tsv")

## Check if any sample is missing
arrow_file <- dir("../run_BM_create_arrow/", pattern = "arrow$", full = T)
sample_name_x <- str_replace(basename(arrow_file), ".arrow", "")
sample_sheet[(sample_name %in% sample_name_x)]

#+ eval=FALSE
## create an ArchRProject
proj <- ArchRProject(
    ArrowFiles = arrow_file,
    outputDirectory = "./tmp/ArchRProject_cooc",
    copyArrows = TRUE
)

## Add metadata
old_meta <- proj@cellColData
old_meta$cell_id_ar = rownames(old_meta)
old_meta = as.data.table(old_meta)
old_meta[, cell_id := sub(".*#", "", cell_id_ar)]

new_meta <- merge(old_meta, d_ann, by = "cell_id")
new_meta = merge(new_meta, sample_sheet, by.x = "Sample", by.y = "sample_name")

## subset the project to only cells with annotation
proj <- subsetArchRProject(
    ArchRProj = proj,
    cells = new_meta$cell_id_ar,
    outputDirectory = "./tmp/ArchRProject_cooc_final",
    dropCells = TRUE,
    force = TRUE
)
proj

new_meta = DataFrame(new_meta)

rownames(new_meta) <- new_meta$cell_id_ar

proj@cellColData <- new_meta

saveArchRProject(proj, outputDirectory = "./tmp/ArchRProject_cooc_final", load = FALSE, overwrite = TRUE)

write_tsv(data.frame(new_meta), "./tmp/cell_metadata_cooc_4.tsv")

#+ eval=TRUE
proj = loadArchRProject("./tmp/ArchRProject_cooc_final", force = TRUE)
## Create subset

d_meta <- data.table(data.frame(proj@cellColData))
d_meta$cell_id_ar <- rownames(proj@cellColData)

d_meta[, antibody_target := paste0(antibody_target)]
d_meta[, co_antibody_target := paste0(co_antibody_target)]
d_meta[, antibody_combination := paste(sort(c(antibody_target, co_antibody_target)), collapse = "_"), by = cell_id_ar]
d_meta$antibody_combination %>% table()


d_meta[antibody_combination == "H3K27me3_H3K4me1"] %>% nrow()
d_meta[antibody_combination == "H3K27me3_H3K4me2"] %>% nrow()
d_meta[antibody_combination == "H3K27me3_H3K4me3"] %>% nrow()

## select the H3K27me3_H3K4me1 only
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_combination == "H3K27me3_H3K4me1", cell_id_ar],
    outputDirectory = "./tmp/H3K27me3_H3K4me1_coCnT_all",
    dropCells = TRUE,
    force = TRUE
)

## select the H3K27me3_H3K4me2 only
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_combination == "H3K27me3_H3K4me2", cell_id_ar],
    outputDirectory = "./tmp/H3K27me3_H3K4me2_coCnT_all",
    dropCells = TRUE,
    force = TRUE
)

## select the H3K27me3_H3K4me3 only
subsetArchRProject(
    ArchRProj = proj,
    cells = d_meta[antibody_combination == "H3K27me3_H3K4me3", cell_id_ar],
    outputDirectory = "./tmp/H3K27me3_H3K4me3_coCnT_all",
    dropCells = TRUE,
    force = TRUE
)

