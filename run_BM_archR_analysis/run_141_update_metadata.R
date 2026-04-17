# =============================================================================
# Script: run_141_update_metadata.R
#
# Input:
#   - ArchR project directory: ./tmp/H3K27me3_coCnT_final3/
#     (expected to contain cell-level metadata and UMAP embeddings for H3K27me3)
#
# Output:
#   - UMAP coordinate table: ./tmp/table_umap_h3k4me27_final3.tsv
#
# Description:
#   - Loads the H3K27me3 ArchR project.
#   - Extracts cell metadata and UMAP (Harmony) embeddings.
#   - Combines metadata and embeddings for downstream usage.
#   - Saves the UMAP coordinates as a TSV file.
# =============================================================================

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(ArchR)
color = ArchRPalettes$stallion

proj = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")

d_meta = as.data.table(proj@cellColData)
d_umap = as.data.table(getEmbedding(proj, embedding = "UMAP_Harmony"))
names(d_umap) = c("UMAP_1", "UMAP_2")

d_meta = cbind(d_meta, d_umap)

write_tsv(d_meta, "./tmp/table_metadata_h3k27me3_final3.tsv")

## save the metadata for H3K4me1-3
proj_k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final3/")
proj_k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final3/")
proj_k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final3/")
proj_k4me1_cooc = loadArchRProject("./tmp/H3K4me1_coCnT_cooc_final3/")
proj_k4me2_cooc = loadArchRProject("./tmp/H3K4me2_coCnT_cooc_final3/")
proj_k4me3_cooc = loadArchRProject("./tmp/H3K4me3_coCnT_cooc_final3/")

d_meta_k4me1 = as.data.table(proj_k4me1@cellColData)
d_meta_k4me2 = as.data.table(proj_k4me2@cellColData)
d_meta_k4me3 = as.data.table(proj_k4me3@cellColData)
d_meta_k4me1_cooc = as.data.table(proj_k4me1_cooc@cellColData)
d_meta_k4me2_cooc = as.data.table(proj_k4me2_cooc@cellColData)
d_meta_k4me3_cooc = as.data.table(proj_k4me3_cooc@cellColData)

write_tsv(d_meta_k4me1, "./tmp/table_metadata_h3k4me1_final3.tsv")
write_tsv(d_meta_k4me2, "./tmp/table_metadata_h3k4me2_final3.tsv")
write_tsv(d_meta_k4me3, "./tmp/table_metadata_h3k4me3_final3.tsv")
write_tsv(d_meta_k4me1_cooc, "./tmp/table_metadata_h3k4me1_final3_with_cooc.tsv")
write_tsv(d_meta_k4me2_cooc, "./tmp/table_metadata_h3k4me2_final3_with_cooc.tsv")
write_tsv(d_meta_k4me3_cooc, "./tmp/table_metadata_h3k4me3_final3_with_cooc.tsv")
    

