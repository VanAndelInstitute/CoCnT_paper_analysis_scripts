##################################################
## Script: run_n_rna_map2.R
## Purpose: Generate heatmaps to visualize the mean scaled expression of marker genes across clusters for H3K4me3
## INPUT:
## - Seurat object containing RNA-seq data with cluster annotations:
## * ../data/triana_2021/WTA_projected.rds
## - Table of marker genes for each cluster and direction (up/down) for each dataset:
## * ./tmp/table_marker_genes_all_clusters_both_directions_round2.tsv
## OUTPUT:
## - Heatmap showing the mean scaled expression of marker genes across clusters:
## * ./figures/116_figure_marker_gene_heatmap_ordered_Reds_H3K4me1.pdf
## * ./figures/116_figure_marker_gene_heatmap_ordered_Reds_H3K4me2.pdf
## * ./figures/116_figure_marker_gene_heatmap_ordered_Reds_H3K4me3.pdf
##

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(Seurat)
library(ArchR)

color = ArchRPalettes$stallion
names(color) = paste0("C", 1:length(color))

so = read_rds("../../data/triana_2021/WTA_projected.bk.rds")
Idents(so) = so@meta.data$ct

x = so@meta.data$ct %>% table
selected = names(x[x > 10])

DimPlot(so, label = F, group.by = "ct")

## subset so
so_sub = subset(so, idents = selected)

celltype_order <- c(
    # Stem & multipotent progenitors
    "HSCs & MPPs",

    # B cell lineage
    "Lymphomyeloid prog",
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
    "Classical Monocytes",
    "Non-classical monocytes",
    "Conventional dendritic cell 1",
    "Conventional dendritic cell 2",

    # Erythroid / Megakaryocytic lineage
    "Erythro-myeloid progenitors",
    "Early erythroid progenitor",
    "Late erythroid progenitor",
    # "Aberrant erythroid",
    "Megakaryocyte progenitors",

    # Dendritic lineage
    # T cell lineage
    "CD4+ naive T cells",
    "CD4+ memory T cells",
    # "CD69+PD-1+ memory CD4+ T cells",
    "CD8+ naive T cells",
    "CD8+ central memory T cells",
    "CD8+ effector memory T cells",
    "CD8+CD103+ tissue resident memory T cells",
    "GammaDelta T cells",
    # "NK T cells",

    # NK lineage
    "NK cell progenitors",
    "CD56brightCD16- NK cells",
    "CD56dimCD16+ NK cells",

    "Eosinophil-basophil-mast cell progenitors",
    "Plasmacytoid dendritic cell progenitors",
    "Plasmacytoid dendritic cells"
)

ordered_cluster = c(
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


## Prepare the marker gene list
d_marker_all_genes = fread("./tmp/table_marker_genes_all_clusters_both_directions_round2.tsv")
d_marker_all_genes$Cluster %>% table

d_marker_gene = d_marker_all_genes[Direction == "Up" & marker == "K4me3" & Log2FC > 1.25, .(Cluster, gene_name = name, Log2FC, FDR)]
## choose the top 50 genes per cluster with smallest FDR
d_marker_gene = d_marker_gene[, head(.SD, 50), by = .(Cluster)]

d_marker_gene[Cluster == "C3"]$gene_name

d_marker_gene$Cluster %<>% factor(., levels = names(color))
doup_gene = d_marker_gene[, .N, by = .(gene_name)][N > 1, gene_name]
d_marker_gene = d_marker_gene[!gene_name %in% doup_gene]
d_marker_gene = d_marker_gene[order(Cluster)]

## get genes of interest
genes_of_interest = d_marker_gene$gene_name

## get data matrix
expr_mat = GetAssayData(so_sub, assay = "RNA", slot = "data") %>% as.matrix
genes_in_data = rownames(expr_mat)
expr_mat_sub = expr_mat[intersect(genes_of_interest, genes_in_data), ]
colnames(expr_mat_sub) = Idents(so_sub)
# colnames(expr_mat_sub) = so@meta.data$ct
table(Idents(so_sub), so_sub@meta.data$ct)

d_plot = expr_mat_sub %>% melt %>% data.table
names(d_plot)[1:2] = c("gene_name", "cluster")
d_plot = d_plot[, .(mean_expr = mean(value)), by = .(gene_name, cluster)]

## scale by gene
d_plot[, mean_expr := scale(mean_expr, center = F), by = gene_name]
d_plot = merge(d_plot, d_marker_gene, by = "gene_name")
d_plot$gene_name = factor(d_plot$gene_name, levels = d_marker_gene$gene_name)
d_plot$cluster %<>% factor(levels = rev(celltype_order))

d_plot_summary = d_plot[, .(index = mean(mean_expr)), by = .(Cluster, cluster)]

## Only average heatmap
m_plot = dcast(d_plot_summary, cluster ~ Cluster, value.var = "index")
rowname = m_plot$cluster
m_plot = as.matrix(m_plot[, -1])
rownames(m_plot) = rowname
library(pheatmap)

pdf("./figures/116_figure_marker_gene_heatmap_ordered_Reds_H3K4me3.pdf", width = 8, height = 8)
pheatmap(
    m_plot[nrow(m_plot):1, ordered_cluster],
	# rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
	colorRampPalette(c("white", RColorBrewer::brewer.pal(9, "Reds")))(100),
	cluster_rows = F,
	cluster_cols = F,
	fontsize_row = 8,
	fontsize_col = 8,
	main = "Mean Scaled Expression of Marker Genes Across Clusters"
)
dev.off()

## H3K4me1
d_marker_gene = d_marker_all_genes[Direction == "Up" & marker == "K4me1" & Log2FC > 1.25, .(Cluster, gene_name = name, Log2FC, FDR)]
## choose the top 50 genes per cluster
d_marker_gene = d_marker_gene[, head(.SD, 50), by = .(Cluster)]

d_marker_gene$Cluster %<>% factor(., levels = names(color))
doup_gene = d_marker_gene[, .N, by = .(gene_name)][N > 1, gene_name]
d_marker_gene = d_marker_gene[!gene_name %in% doup_gene]
d_marker_gene = d_marker_gene[order(Cluster)]

## get genes of interest
genes_of_interest = d_marker_gene$gene_name

## get data matrix
expr_mat = GetAssayData(so_sub, assay = "RNA", slot = "data") %>% as.matrix
dim(expr_mat)
genes_in_data = rownames(expr_mat)
expr_mat_sub = expr_mat[intersect(genes_of_interest, genes_in_data), ]
colnames(expr_mat_sub) = Idents(so_sub)
# colnames(expr_mat_sub) = so@meta.data$ct
table(Idents(so_sub), so_sub@meta.data$ct)

d_plot = expr_mat_sub %>% melt %>% data.table
names(d_plot)[1:2] = c("gene_name", "cluster")
d_plot = d_plot[, .(mean_expr = mean(value)), by = .(gene_name, cluster)]

## scale by gene
d_plot[, mean_expr := scale(mean_expr, center = F), by = gene_name]
d_plot = merge(d_plot, d_marker_gene, by = "gene_name")
d_plot$gene_name = factor(d_plot$gene_name, levels = d_marker_gene$gene_name)
# d_plot$Cluster = factor(d_plot$Cluster, levels = 0:10)
# d_plot$cluster %<>% factor(levels = cluster_name[res$order])
d_plot$cluster %<>% factor(levels = rev(celltype_order))

d_plot_summary = d_plot[, .(index = mean(mean_expr)), by = .(Cluster, cluster)]

## Only average heatmap
m_plot = dcast(d_plot_summary, cluster ~ Cluster, value.var = "index")
rowname = m_plot$cluster
m_plot = as.matrix(m_plot[, -1])
rownames(m_plot) = rowname
library(pheatmap)
pdf("./figures/116_figure_marker_gene_heatmap_ordered_Reds_H3K4me1.pdf", width = 8, height = 8)
pheatmap(
    m_plot[nrow(m_plot):1, ordered_cluster[ordered_cluster %in% colnames(m_plot)]],
	# rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
	colorRampPalette(c("white", RColorBrewer::brewer.pal(9, "Reds")))(100),
	cluster_rows = F,
	cluster_cols = F,
	fontsize_row = 8,
	fontsize_col = 8,
	main = "Mean Scaled Expression of Marker Genes Across Clusters"
)
dev.off()



## H3K4me2
d_marker_gene = d_marker_all_genes[Direction == "Up" & marker == "K4me2" & Log2FC > 1.25, .(Cluster, gene_name = name, Log2FC, FDR)]
## choose the top 50 genes per cluster
d_marker_gene = d_marker_gene[, head(.SD, 50), by = .(Cluster)]

d_marker_gene$Cluster %<>% factor(., levels = names(color))
doup_gene = d_marker_gene[, .N, by = .(gene_name)][N > 1, gene_name]
d_marker_gene = d_marker_gene[!gene_name %in% doup_gene]
d_marker_gene = d_marker_gene[order(Cluster)]

## get genes of interest
genes_of_interest = d_marker_gene$gene_name

## get data matrix
expr_mat = GetAssayData(so_sub, assay = "RNA", slot = "data") %>% as.matrix
dim(expr_mat)
genes_in_data = rownames(expr_mat)
expr_mat_sub = expr_mat[intersect(genes_of_interest, genes_in_data), ]
colnames(expr_mat_sub) = Idents(so_sub)
# colnames(expr_mat_sub) = so@meta.data$ct
table(Idents(so_sub), so_sub@meta.data$ct)

d_plot = expr_mat_sub %>% melt %>% data.table
names(d_plot)[1:2] = c("gene_name", "cluster")
d_plot = d_plot[, .(mean_expr = mean(value)), by = .(gene_name, cluster)]

## scale by gene
d_plot[, mean_expr := scale(mean_expr, center = F), by = gene_name]
d_plot = merge(d_plot, d_marker_gene, by = "gene_name")
d_plot$gene_name = factor(d_plot$gene_name, levels = d_marker_gene$gene_name)
# d_plot$Cluster = factor(d_plot$Cluster, levels = 0:10)
# d_plot$cluster %<>% factor(levels = cluster_name[res$order])
d_plot$cluster %<>% factor(levels = rev(celltype_order))

d_plot_summary = d_plot[, .(index = mean(mean_expr)), by = .(Cluster, cluster)]

## Only average heatmap
m_plot = dcast(d_plot_summary, cluster ~ Cluster, value.var = "index")
rowname = m_plot$cluster
m_plot = as.matrix(m_plot[, -1])
rownames(m_plot) = rowname
library(pheatmap)


pdf("./figures/116_figure_marker_gene_heatmap_ordered_Reds_H3K4me2.pdf", width = 8, height = 8)
pheatmap(
    m_plot[nrow(m_plot):1, ordered_cluster[ordered_cluster %in% colnames(m_plot)]],
	# rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
	colorRampPalette(c("white", RColorBrewer::brewer.pal(9, "Reds")))(100),
	cluster_rows = F,
	cluster_cols = F,
	fontsize_row = 8,
	fontsize_col = 8,
	main = "Mean Scaled Expression of Marker Genes Across Clusters"
)
dev.off()



