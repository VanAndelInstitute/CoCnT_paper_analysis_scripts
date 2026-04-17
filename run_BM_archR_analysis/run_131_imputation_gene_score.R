################################################################################
# Script: run_120_imputation_gene_score.R
# Purpose: Impute gene scores across multiple histone mark datasets (H3K27me3, 
#          H3K4me1, H3K4me2, H3K4me3) for cell lineage analysis
#
# INPUT:
#   - ArchR projects (directories):
#     * ./tmp/H3K27me3_coCnT_final2/
#     * ./tmp/H3K4me1_coCnT_final2/
#     * ./tmp/H3K4me2_coCnT_final2/
#     * ./tmp/H3K4me3_coCnT_final2/
#   - External library: ../../lib/R/lib_imputation_gene_score.R
#
# OUTPUT:
#   - ArchR project directories (subsetted):
#     * ./tmp/H3K4me1_coCnT_final3/
#     * ./tmp/H3K4me2_coCnT_final3/
#     * ./tmp/H3K4me3_coCnT_final3/
#     * ./tmp/H3K27me3_coCnT_final3/
#   - RDS files (gene score matrices):
#     * ./tmp/matrix_original_gene_score_h3k27me3.rds
#     * ./tmp/matrix_original_gene_score_h3k4me1.rds
#     * ./tmp/matrix_original_gene_score_h3k4me2.rds
#     * ./tmp/matrix_original_gene_score_h3k4me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me1_from_h3k27me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me2_from_h3k27me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me3_from_h3k27me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k27me3.rds
################################################################################

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

# all lineage 
all_cell_clusters = c(
    "C10",

    "C16",
    "C17",
    "C3",
    "C5",
    "C6",
    "C4",


    "C12",
    "C13",
    "C14",
    
    "C18",
    "C8",
    "C7",

    "C15",
    "C19",
    "C11",

    "C1",
    "C9",
    "C2"
    )


## Load project
# INPUT: ArchR projects from final2 round
proj_h3k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final2/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final2/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final2/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final2/")

## Get the metadata
d_meta_k27 = as.data.table(proj_h3k27@cellColData)
d_meta_h3k4me1 = as.data.table(proj_h3k4me1@cellColData)
d_meta_h3k4me2 = as.data.table(proj_h3k4me2@cellColData)
d_meta_h3k4me3 = as.data.table(proj_h3k4me3@cellColData)

## Merge the metadata
# Create mappings between H3K27me3 and other histone marks for matching cells
map_k27tok4me1 = merge(
    d_meta_k27[, .(cell_id = cell_id, cell_id_ar = cell_id_ar, Clusters_plus)],
    d_meta_h3k4me1[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    by = "cell_id", all = F
)[Clusters_plus %in% all_cell_clusters]

map_k27tok4me2 = merge(
    d_meta_k27[, .(cell_id = cell_id, cell_id_ar = cell_id_ar, Clusters_plus)],
    d_meta_h3k4me2[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    by = "cell_id", all = F
)[Clusters_plus %in% all_cell_clusters]

map_k27tok4me3 = merge(
    d_meta_k27[, .(cell_id = cell_id, cell_id_ar = cell_id_ar, Clusters_plus)],
    d_meta_h3k4me3[, .(cell_id = cell_id, cell_id_ar = cell_id_ar)],
    by = "cell_id", all = F
)[Clusters_plus %in% all_cell_clusters]

## Subset the datasets
# OUTPUT: Subsetted ArchR projects (final3) containing only matched cells from all lineage clusters
subsetArchRProject(
    proj_h3k4me1,
    cells = map_k27tok4me1$cell_id_ar.y,
    outputDirectory = "./tmp/H3K4me1_coCnT_final3/",
    dropCells = T,
    force = T
)

subsetArchRProject(
    proj_h3k4me2,
    cells = map_k27tok4me2$cell_id_ar.y,
    outputDirectory = "./tmp/H3K4me2_coCnT_final3/",
    dropCells = T,
    force = T
)

subsetArchRProject(
    proj_h3k4me3,
    cells = map_k27tok4me3$cell_id_ar.y,
    outputDirectory = "./tmp/H3K4me3_coCnT_final3/",
    dropCells = T,
    force = T
)

subsetArchRProject(
    proj_h3k27,
    cells = unique(c(
	    map_k27tok4me1$cell_id_ar.x,
	    map_k27tok4me2$cell_id_ar.x,
	    map_k27tok4me3$cell_id_ar.x
	    )),
    outputDirectory = "./tmp/H3K27me3_coCnT_final3/",
    dropCells = T,
    force = T
)

# Reload subsetted projects
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final3/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final3/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final3/")
proj_h3k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")

httpgd::hgd(port = 4322)

## save original gene score matrix
# OUTPUT: Original (non-imputed) gene score matrices from each histone mark dataset
m = getMatrixFromProject(proj_h3k27, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k27_orig = assay(m)
rownames(m_h3k27_orig) = gene_name
write_rds(m_h3k27_orig, "./tmp/matrix_original_gene_score_h3k27me3.rds")

m = getMatrixFromProject(proj_h3k4me1, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k4me1_orig = assay(m)
rownames(m_h3k4me1_orig) = gene_name
write_rds(m_h3k4me1_orig, "./tmp/matrix_original_gene_score_h3k4me1.rds")

m = getMatrixFromProject(proj_h3k4me2, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k4me2_orig = assay(m)
rownames(m_h3k4me2_orig) = gene_name
write_rds(m_h3k4me2_orig, "./tmp/matrix_original_gene_score_h3k4me2.rds")

m = getMatrixFromProject(proj_h3k4me3, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k4me3_orig = assay(m)
rownames(m_h3k4me3_orig) = gene_name
write_rds(m_h3k4me3_orig, "./tmp/matrix_original_gene_score_h3k4me3.rds")

## save imputed gene score matrix
# OUTPUT: Imputed gene score matrices using Harmony-based imputation weights
# H3K4me1, H3K4me2, H3K4me3 imputed from H3K27me3 reference; H3K27me3 self-imputed

proj_h3k27 = addImputeWeights(proj_h3k27, reducedDims = "Harmony")
m = getMatrixFromProject(proj_h3k27, useMatrix = "GeneScoreMatrix")
gene_name = m@elementMetadata$name
m_h3k27 = imputeMatrix(mat = assay(m), getImputeWeights(proj_h3k27))
rownames(m_h3k27) = gene_name

m_h3k4me1 = get_imputed_gene_score(proj_h3k27, proj_h3k4me1)
m_h3k4me2 = get_imputed_gene_score(proj_h3k27, proj_h3k4me2)
m_h3k4me3 = get_imputed_gene_score(proj_h3k27, proj_h3k4me3)

write_rds(m_h3k4me1, "./tmp/matrix_imputed_gene_score_h3k4me1_from_h3k27me3.rds")
write_rds(m_h3k4me2, "./tmp/matrix_imputed_gene_score_h3k4me2_from_h3k27me3.rds")
write_rds(m_h3k4me3, "./tmp/matrix_imputed_gene_score_h3k4me3_from_h3k27me3.rds")
write_rds(m_h3k27, "./tmp/matrix_imputed_gene_score_h3k27me3.rds")
