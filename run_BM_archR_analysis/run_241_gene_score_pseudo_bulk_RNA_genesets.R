library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d = fread("./tmp/table_triana_2021_rna_pseudo_bulk_by_ct.tsv")

celltype_order <- c(
    # Stem & multipotent progenitors
    "HSCs & MPPs",
    "Lymphomyeloid prog",
    "Erythro-myeloid progenitors",

    # B cell lineage
    "Pre-pro-B cells",
    "Pro-B cells",
    "Pre-B cells",
    "Small pre-B cell",
    "Immature B cells",
    "Mature naive B cells",
    "Nonswitched memory B cells",
    "Class switched memory B cells",
    "Plasma cells",

    # Myeloid lineage
    "Early promyelocytes",
    "Late promyelocytes",
    "Myelocytes",
    "Eosinophil-basophil-mast cell progenitors",
    "Classical Monocytes",
    "Non-classical monocytes",

    # Erythroid / Megakaryocytic lineage
    "Early erythroid progenitor",
    "Late erythroid progenitor",
    "Aberrant erythroid",
    "Megakaryocyte progenitors",

    # Dendritic lineage
    "Conventional dendritic cell 1",
    "Conventional dendritic cell 2",
    "Plasmacytoid dendritic cell progenitors",
    "Plasmacytoid dendritic cells",

    # T cell lineage
    "CD4+ naive T cells",
    "CD4+ memory T cells",
    "CD69+PD-1+ memory CD4+ T cells",
    "CD8+ naive T cells",
    "CD8+ central memory T cells",
    "CD8+ effector memory T cells",
    "CD8+CD103+ tissue resident memory T cells",
    "GammaDelta T cells",
    "NK T cells",

    # NK lineage
    "NK cell progenitors",
    "CD56brightCD16- NK cells",
    "CD56dimCD16+ NK cells"
)


d_prc2 = fread("../../data/Gene_Lists/gene_list_PRC2.csv")
d_control = fread("../../data/Gene_Lists/gene_list_control.csv")
d_h3k4_demethylase = fread("../../data/Gene_Lists/gene_list_H3K4_h3k27_demethylase.csv")
d_h3k4_methyltransferase = fread("../../data/Gene_Lists/gene_list_H3K4_methyltransferase.csv")

d_B_cell = data.table(gene = c("PAX5", "EBF1", "TCF3", "RAG1", "RAG2", "PAX5", "DNTT", "ID1", "ID2", "ID3", "ID4"))

draw_heatmap <- function(gene_list, out_pdf) {
  d_sub = d[gene %in% gene_list$gene]
  mat = as.matrix(d_sub[, -1, with = FALSE])
  rownames(mat) = d_sub$gene
  celltype_order_filt = celltype_order[celltype_order %in% colnames(mat)]
  mat = mat[, celltype_order_filt, drop = FALSE]
  pdf(out_pdf, width = 6, height = 8)
  pheatmap::pheatmap(
    mat,
    scale = "row",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    cluster_cols = FALSE,
    fontsize_col = 8,
    color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  )
  dev.off()
}

draw_heatmap(d_prc2, "./figures/241_heatmap_PR2C_genes_pseudo_bulk_RNA_by_cell_type.pdf")
draw_heatmap(d_h3k4_demethylase, "./figures/241_heatmap_H3K4_demethylase_genes_pseudo_bulk_RNA_by_cell_type.pdf")
draw_heatmap(d_h3k4_methyltransferase, "./figures/241_heatmap_H3K4_methyltransferase_genes_pseudo_bulk_RNA_by_cell_type.pdf")
draw_heatmap(d_control, "./figures/241_heatmap_control_genes_pseudo_bulk_RNA_by_cell_type.pdf")

celltype_order <- c(
    # Stem & multipotent progenitors
    "HSCs & MPPs",
    "Lymphomyeloid prog",

    # B cell lineage
    "Pre-pro-B cells",
    "Pro-B cells",
    "Pre-B cells",
    "Small pre-B cell",
    "Immature B cells",
    "Mature naive B cells",
    "Nonswitched memory B cells",
    "Class switched memory B cells",
    "Plasma cells"
)

draw_heatmap(d_B_cell, "./figures/241_heatmap_B_cell_genes_pseudo_bulk_RNA_by_cell_type.pdf")
