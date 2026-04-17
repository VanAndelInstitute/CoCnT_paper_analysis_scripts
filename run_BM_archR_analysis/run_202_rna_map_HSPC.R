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

httpgd::hgd(port = 4323)
so = read_rds("../../data/triana_2021/WTA_projected.bk.rds")
Idents(so) = so@meta.data$ct

x = so@meta.data$ct %>% table
selected = names(x[x > 10])

## subset so
so_sub = subset(so, idents = selected)

d_marker_all_genes = fread("./tmp/table_HSPC_marker_genes_all_markers.tsv")

## Label transfering
d_marker_gene = d_marker_all_genes[Direction == "Up" & marker == "K4me3" & Log2FC > 1.25, .(Cluster, gene_name = name, Log2FC, FDR)]

d_marker_gene[, .N, Cluster]
# d_marker_gene = d_marker_all_genes[Direction == "Down" & marker == "K27" & Log2FC < -1.25, .(Cluster, gene_name = name, Log2FC, FDR)]
## choose the top 100 genes per cluster
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
    "Eosinophil-basophil-mast cell progenitors",
    "Classical Monocytes",
    "Non-classical monocytes",

    # Erythroid / Megakaryocytic lineage
    "Erythro-myeloid progenitors",
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

## scale by gene
d_plot[, mean_expr := scale(mean_expr, center = F), by = gene_name]
d_plot = merge(d_plot, d_marker_gene, by = "gene_name")
d_plot$gene_name = factor(d_plot$gene_name, levels = d_marker_gene$gene_name)
d_plot$cluster %<>% factor(levels = rev(celltype_order))

d_plot_summary = d_plot[, .(index = mean(mean_expr)), by = .(Cluster, cluster)]

## Only average
ggplot(d_plot_summary, aes(y = cluster, x = Cluster, size = index, fill = Cluster)) +
    geom_point(shape = 21) +
    scale_size_continuous(range = c(0, 5)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cluster", y = "Reference", fill = "Mean Scaled Expression") +
    scale_fill_manual(values = color) +
    ggtitle("Mean Scaled Expression of Marker Genes Across Clusters")

## Only average heatmap
m_plot = dcast(d_plot_summary, cluster ~ Cluster, value.var = "index")
rowname = m_plot$cluster
m_plot = as.matrix(m_plot[, -1])
rownames(m_plot) = rowname
library(pheatmap)
p = pheatmap(
    m_plot[nrow(m_plot):1, ],
	colorRampPalette(c("white", RColorBrewer::brewer.pal(9, "Reds")))(100),
	cluster_rows = F,
	cluster_cols = F,
	fontsize_row = 8,
	fontsize_col = 8,
	main = "Mean Scaled Expression of Marker Genes Across Clusters"
)

p

