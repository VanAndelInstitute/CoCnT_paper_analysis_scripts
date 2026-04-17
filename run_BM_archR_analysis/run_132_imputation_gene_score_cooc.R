##################################################
## Script: run_n_imputation_gene_score_cooc.R
## Purpose: Impute gene scores across co-occopancy datasets
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final3/
## * ./tmp/H3K27me3_H3K4me1_coCnT_all/
## * ./tmp/H3K27me3_H3K4me2_coCnT_all/
## * ./tmp/H3K27me3_H3K4me3_coCnT_all/
## OUTPUT:
## - ArchR project directories (subsetted):
## * ./tmp/H3K4me1_coCnT_cooc_final3/
## * ./tmp/H3K4me2_coCnT_cooc_final3/
## * ./tmp/H3K4me3_coCnT_cooc_final3/
## - RDS files (gene score matrices):
## * ./tmp/matrix_original_gene_score_h3k4me1_cooc.rds
## * ./tmp/matrix_original_gene_score_h3k4me2_cooc.rds
## * ./tmp/matrix_original_gene_score_h3k4me3_cooc.rds
## * ./tmp/matrix_imputed_gene_score_h3k4me1_cooc.rds
## * ./tmp/matrix_imputed_gene_score_h3k4me2_cooc.rds
## * ./tmp/matrix_imputed_gene_score_h3k4me3_cooc.rds
##

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 11)
color = ArchRPalettes$stallion

source("../../lib/R/lib_imputation_gene_score.R")

set.seed(123)

## Load project
proj_h3k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")

proj_h3k4me1_cooc = loadArchRProject("./tmp/H3K27me3_H3K4me1_coCnT_all/")
proj_h3k4me2_cooc = loadArchRProject("./tmp/H3K27me3_H3K4me2_coCnT_all/")
proj_h3k4me3_cooc = loadArchRProject("./tmp/H3K27me3_H3K4me3_coCnT_all/")

## Get the metadata
d_meta_k27 = as.data.table(proj_h3k27@cellColData)

d_meta_h3k4me1_cooc = as.data.table(proj_h3k4me1_cooc@cellColData)
d_meta_h3k4me2_cooc = as.data.table(proj_h3k4me2_cooc@cellColData)
d_meta_h3k4me3_cooc = as.data.table(proj_h3k4me3_cooc@cellColData)

## Merge the metadata
map_k27tok4me1_cooc = merge(
    d_meta_k27[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    d_meta_h3k4me1_cooc[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    by = "cell_id", all = F
)

map_k27tok4me2_cooc = merge(
    d_meta_k27[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    d_meta_h3k4me2_cooc[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    by = "cell_id", all = F
)

map_k27tok4me3_cooc = merge(
    d_meta_k27[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    d_meta_h3k4me3_cooc[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    by = "cell_id", all = F
)


## Subset the datasets
subsetArchRProject(
    proj_h3k4me1_cooc,
    cells = map_k27tok4me1_cooc$cell_id_ar.y,
    outputDirectory = "./tmp/H3K4me1_coCnT_cooc_final3/",
    dropCells = T,
    force = T
)

subsetArchRProject(
    proj_h3k4me2_cooc,
    cells = map_k27tok4me2_cooc$cell_id_ar.y,
    outputDirectory = "./tmp/H3K4me2_coCnT_cooc_final3/",
    dropCells = T,
    force = T
)

## BUG: Need to subset based on the intersected cells only
m = proj_h3k4me3_cooc@reducedDims$Harmony
subsetArchRProject(
    proj_h3k4me3_cooc,
    cells = intersect(map_k27tok4me3_cooc$cell_id_ar.y, rownames(m@listData$matDR)),
    outputDirectory = "./tmp/H3K4me3_coCnT_cooc_final3/",
    dropCells = T,
    force = T
)

proj_h3k4me1_cooc = loadArchRProject("./tmp/H3K4me1_coCnT_cooc_final3/")
proj_h3k4me2_cooc = loadArchRProject("./tmp/H3K4me2_coCnT_cooc_final3/")
proj_h3k4me3_cooc = loadArchRProject("./tmp/H3K4me3_coCnT_cooc_final3/")

## save original gene score matrix
m = getMatrixFromProject(proj_h3k4me1_cooc, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k4me1_orig = assay(m)
rownames(m_h3k4me1_orig) = gene_name
write_rds(m_h3k4me1_orig, "./tmp/matrix_original_gene_score_h3k4me1_cooc.rds")

m = getMatrixFromProject(proj_h3k4me2_cooc, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k4me2_orig = assay(m)
rownames(m_h3k4me2_orig) = gene_name
write_rds(m_h3k4me2_orig, "./tmp/matrix_original_gene_score_h3k4me2_cooc.rds")

m = getMatrixFromProject(proj_h3k4me3_cooc, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k4me3_orig = assay(m)
rownames(m_h3k4me3_orig) = gene_name
write_rds(m_h3k4me3_orig, "./tmp/matrix_original_gene_score_h3k4me3_cooc.rds")



m_h3k4me1_cooc = get_imputed_gene_score(proj_h3k27, proj_h3k4me1_cooc)
m_h3k4me2_cooc = get_imputed_gene_score(proj_h3k27, proj_h3k4me2_cooc)
m_h3k4me3_cooc = get_imputed_gene_score(proj_h3k27, proj_h3k4me3_cooc)

write_rds(m_h3k4me1_cooc, "./tmp/matrix_imputed_gene_score_h3k4me1_cooc.rds")
write_rds(m_h3k4me2_cooc, "./tmp/matrix_imputed_gene_score_h3k4me2_cooc.rds")
write_rds(m_h3k4me3_cooc, "./tmp/matrix_imputed_gene_score_h3k4me3_cooc.rds")


