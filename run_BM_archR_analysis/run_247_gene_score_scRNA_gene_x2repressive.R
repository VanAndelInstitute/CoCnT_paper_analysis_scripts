library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(Seurat)
library(slingshot)

# httpgd::hgd(port = 4323)

cpg_gene = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")

celltype_order_l <- list(
    b = c(
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
	),
    e = c(
	"HSCs & MPPs",
	# Erythroid / Megakaryocytic lineage
	"Erythro-myeloid progenitors",
	"Early erythroid progenitor",
	"Late erythroid progenitor"
	# "Aberrant erythroid"
	# "Megakaryocyte progenitors"
	),
    m = c(
	"HSCs & MPPs",
	# Myeloid lineage
	"Early promyelocytes",
	"Late promyelocytes",
	"Myelocytes",
	# "Eosinophil-basophil-mast cell progenitors",
	"Classical Monocytes"
	# "Non-classical monocytes"
    ),
    t_nk = c(
	"HSCs & MPPs",
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
)


## Draw the umap for each lineage {{{

so <- read_rds("../../data/triana_2021/WTA_projected.bk.rds")
Idents(so) <- so@meta.data$ct

pdf("./figures/247_dimplot_scRNASeq_umap.pdf", width = 8, height = 5)
DimPlot(so, group.by = "ct", label = T) + theme(legend.position = "none")
so

so_sub <- subset(so, idents = celltype_order_l[["b"]])
DimPlot(so_sub, group.by = "ct", label = T)

so_sub <- subset(so, idents = celltype_order_l[["e"]])
DimPlot(so_sub, group.by = "ct", label = T)

so_sub <- subset(so, idents = celltype_order_l[["m"]])
DimPlot(so_sub, group.by = "ct", label = T)

dev.off()
## }}}

## Pseudobulk gene expression per cluster {{{

## check the dataset available in the Seurat object

options(future.globals.maxSize = 10000 * 1024^2)

d = so[["RNA"]]@data
# d <- so[["RNA"]]@scale.data
# so_x = so
# so_x = SCTransform(so_x, verbose = F, return.only.var.genes = F)

dim(d)

pdf("./figures/247_x2repressive_gene_expression.pdf", width = 6, height = 4)
for (i in 1:3) {
    print(i)

    x_group = c("B_consensus", "E_consensus", "M_consensus")[i]
    x_celltype_order = celltype_order_l[[i]]

    so_sub = subset(so, idents = x_celltype_order)

    pb <- AggregateExpression(
	so_sub,
	assays = "RNA",
	group.by = c("ct"),
	return.seurat = TRUE,
	normalization.method = "LogNormalize"
    )
    mat <- GetAssayData(pb, assay = "RNA", layer = "data")

    table(so@meta.data$ct)

    ## Un->repressive, and Bivalent->repressive 
    d_stat_trans <- fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")
    gene_un_2repressive <- d_stat_trans[trans %in% c("Un->Repressive") & group == x_group, gene] %>% intersect(rownames(mat))
    gene_bi_2repressive <- d_stat_trans[trans %in% c("Bivalent->Repressive") & group == x_group, gene] %>% intersect(rownames(mat))

    length(gene_un_2repressive)
    length(gene_bi_2repressive)

    x <- so_sub@meta.data$ct

    # (exp(mat) - 1) / (exp(mat[, "HSCs & MPPs"]) - 1) -> mat

    y <- mat[gene_un_2repressive, ]
    d_plot_un2a = melt(y)
    d_plot_un2a$cat = "un2a"

    y <- mat[gene_bi_2repressive, ]
    d_plot_bi2a = melt(y)
    d_plot_bi2a$cat = "bi2a"

    # y = mat
    # d_plot_all = melt(y)
    # d_plot_all$cat = "all"

    d_plot <- rbind(d_plot_un2a, d_plot_bi2a)
    # d_plot <- rbind(d_plot, d_plot_all)
    d_plot <- data.table(d_plot)
    names(d_plot) <- c("gene", "ct", "y", "cat")


    celltype_order = celltype_order_l[[i]]
    d_plot_b <- d_plot[ct %in% celltype_order]
    d_plot_b$ct <- factor(d_plot_b$ct, levels = celltype_order)
    # d_plot_b[, is_cpg := ifelse(gene %in% cpg_gene, "cpg", "non-cpg")]
    d_plot_b$lineage = names(celltype_order_l)[i]

    library(plyr)
    p <- ggplot(d_plot_b, aes(x = ct, y = y, fill = cat)) +
	geom_boxplot(outliers = F) +
	theme_classic() +
	scale_fill_manual(values = c("un2a" = "#1872b5", "bi2a" = "#e18727"), name = "Trajectory") +
	# facet_wrap(~cat, scales = "free_y") +
	xlab(paste0(names(celltype_order_l)[i])) +
	ylab("scRNASeq gene expression") +
	scale_y_log10() +
	stat_compare_means(aes(group = cat), method = "wilcox.test", label = "p.signif") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	ggtitle(paste0(names(celltype_order_l)[i], " lineage with ", length(gene_un_2repressive), " un2a genes and ", length(gene_bi_2repressive), " bi2a genes"))
    print(p)
}

dev.off()



d_b_l = list()

pdf("./figures/247_plot_x2repressive_fold_change > 2.pdf", width = 6, height = 4)
for (i in 1:3) {
    print(i)

    x_group = c("B_consensus", "E_consensus", "M_consensus")[i]
    x_celltype_order = celltype_order_l[[i]]

    so_sub = subset(so, idents = x_celltype_order)

    pb <- AggregateExpression(
	so_sub,
	assays = "RNA",
	group.by = c("ct"),
	return.seurat = TRUE,
	normalization.method = "LogNormalize"
    )
    mat <- GetAssayData(pb, assay = "RNA", layer = "data")

    table(so@meta.data$ct)

    ## Un->repressive, and Bivalent->repressive 
    d_stat_trans <- fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")
    gene_un_2repressive <- d_stat_trans[trans %in% c("Un->Repressive") & group == x_group, gene] %>% intersect(rownames(mat))
    gene_bi_2repressive <- d_stat_trans[trans %in% c("Bivalent->Repressive") & group == x_group, gene] %>% intersect(rownames(mat))

    length(gene_un_2repressive)
    length(gene_bi_2repressive)

    x <- so_sub@meta.data$ct

    (exp(mat) - 1) / (exp(mat[, "HSCs & MPPs"]) - 1) -> mat

    y <- mat[gene_un_2repressive, ]
    d_plot_un2a = melt(y)
    d_plot_un2a$cat = "un2a"

    y <- mat[gene_bi_2repressive, ]
    d_plot_bi2a = melt(y)
    d_plot_bi2a$cat = "bi2a"

    # y = mat
    # d_plot_all = melt(y)
    # d_plot_all$cat = "all"

    d_plot <- rbind(d_plot_un2a, d_plot_bi2a)
    # d_plot <- rbind(d_plot, d_plot_all)
    d_plot <- data.table(d_plot)
    names(d_plot) <- c("gene", "ct", "y", "cat")

    mat["PAX5", ]

    celltype_order = celltype_order_l[[i]]
    d_plot_b <- d_plot[ct %in% celltype_order]
    d_plot_b$ct <- factor(d_plot_b$ct, levels = celltype_order)
    # d_plot_b[, is_cpg := ifelse(gene %in% cpg_gene, "cpg", "non-cpg")]
    d_plot_b$lineage = names(celltype_order_l)[i]

    d_b_l[[i]] <- d_plot_b[y > 2]

    ## The fold change profile of genes 
    p = ggplot(d_plot_b, aes(x = ct, y = y, fill = cat)) + 
	geom_boxplot(outliers = F) +
	theme_classic() +
	scale_fill_manual(values = c("un2a" = "#1872b5", " bi2a" = "#e18727"), name = "Trajectory") +
	xlab(paste0(names(celltype_order_l)[i])) +
	ylab("Fold change of gene expression (compared to HSCs & MPP )") +
	scale_y_log10() +
	stat_compare_means(aes(group = cat), method = "wilcox.test", label = "p.signif") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)


    ## How many genes has fold change > 2 in each cell type?
    d_plot_c = d_plot_b[, .(num_gene = .N, num_gene_fd2 = sum(y > 2, na.rm = T)), by = .(ct, cat)]
    d_plot_c[, prop_fd2 := num_gene_fd2 / num_gene]
    p = ggplot(d_plot_c[ct != "HSCs & MPPs"], aes(x = ct, y = prop_fd2, fill = cat)) +
	geom_bar(stat = "identity", position = "dodge") +
	theme_classic() +
	scale_fill_manual(values = c("un2a" = "#1872b5", "bi2a" = "#e18727"), name = "Trajectory") +
	xlab(paste0(names(celltype_order_l)[i])) +
	ylab("Proportion of genes with fold > 2") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

}

dev.off()

d_b_l

## }}}

## the pseudotime for three lineages {{{

so <- read_rds("../../data/triana_2021/WTA_projected.bk.rds")
so@meta.data
Idents(so) <- so@meta.data$ct


d = so@meta.data
# mat = so@assays$RNA@scale.data
mat = so@assays$RNA@data

lineage_v = c("B.cells", "Erythroid", "Myelocytes")

pdf("./figures/247_scatter_plot_x2repressive_RNA_pseudotime.pdf", width = 6, height = 4)

for (i in 1:3) {
    print(i)

    psudotime_v = d[, lineage_v[i]]
    x_group = c("B_consensus", "M_consensus", "E_consensus")[i]
    x_celltype_order = celltype_order_l[[i]]
    d_activated_genes = d_b_l[[i]][y > 2]

    # d_stat_trans <- fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")
    # gene_un_2repressive <- d_stat_trans[trans %in% c("Un->repressive", "Un->repressive->Un") & group == x_group, gene] %>% intersect(rownames(mat))
    # gene_bi_2repressive <- d_stat_trans[trans %in% c("Bivalent->repressive", "Bivalent->repressive->Bivalent") & group == x_group, gene] %>% intersect(rownames(mat))


    gene_un_2repressive = d_activated_genes[cat == "un2a", gene] %>% intersect(rownames(mat)) %>% unique
    gene_bi_2repressive = d_activated_genes[cat == "bi2a", gene] %>% intersect(rownames(mat)) %>% unique
    # gene_all = d_activated_genes[cat == "all", gene] %>% intersect(rownames(mat)) %>% unique

    y_u2a = mat[gene_un_2repressive, ] %>% apply(2, mean)
    y_b2a = mat[gene_bi_2repressive, ] %>% apply(2, mean)
    # y_all = mat[gene_all, ] %>% apply(2, mean)


    d_plot = data.table(
	ct = d$ct,
	psudotime = psudotime_v,
	Un2repressive = y_u2a,
	Bi2repressive = y_b2a
	# All = y_all
    ) 
    d_plot = d_plot[ct %in% x_celltype_order]
    d_plot$psudotime %<>% rank(ties = "random")
    d_plot$psudotime = d_plot$psudotime / max(d_plot$psudotime) * 100
    d_plot = melt(d_plot, id.vars = c("psudotime", "ct"), variable.name = "cat", value.name = "y")

    p = ggplot(d_plot, aes(x = psudotime, y = y, color = cat)) +
	geom_point(alpha = 0.2, size = 0.1) +
	geom_smooth(method = "loess", se = F, span = 0.7) +
	theme_classic() +
	xlab(paste0(names(celltype_order_l)[i], " pseudotime")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

    d_plot = data.table(
	ct = d$ct,
	psudotime = psudotime_v,
	Un2repressive = scale(y_u2a),
	Bi2repressive = scale(y_b2a)
	# All = scale(y_all)
    ) 
    d_plot = d_plot[ct %in% x_celltype_order]
    d_plot$psudotime %<>% rank(ties = "random")
    d_plot$psudotime = d_plot$psudotime / max(d_plot$psudotime) * 100
    d_plot = melt(d_plot, id.vars = c("psudotime", "ct"), variable.name = "cat", value.name = "y")

    p = ggplot(d_plot, aes(x = psudotime, y = y, color = cat)) +
	geom_point(alpha = 0.2, size = 0.1) +
	geom_smooth(method = "loess", se = F, span = 0.7) +
	scale_color_manual(values = c("Bi2repressive.V1" = "#e18726", "Un2repressive.V1" = "#1972b5")) +
	theme_classic() +
	xlab(paste0(names(celltype_order_l)[i], " pseudotime")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
}
dev.off()



## }}}

## the pseudotime for specific genes {{{

so <- read_rds("../../data/triana_2021/WTA_projected.bk.rds")
so@meta.data
Idents(so) <- so@meta.data$ct


d = so@meta.data
# mat = so@assays$RNA@scale.data
so
mat = so@assays$RNA@data

lineage_v = c("B.cells", "Erythroid", "Myelocytes")

gene_l = list(
    b = c("PAX5", "EBF1", "RAG1", "DNTT"),
    e = c("GATA1", "HBA1", "HBA2", "KLF1"),
    m = c("CEBPA", "MPO", "S100A8", "LYZ") 
)
lineage_v = c("B.cells", "Erythroid", "Myelocytes")

pdf("./figures/247_gene_scRNA_pseudotime.pdf")
for (i in 1:3) {
    gene_v = gene_l[[i]]

    psudotime_v = d[, lineage_v[i]]
    x_celltype_order = celltype_order_l[[i]]

    so@meta.data
    Idents(so) <- so@meta.data$ct
    d = so@meta.data
    # mat = so@assays$RNA@scale.data
    mat = so@assays$RNA@data

    d_plot = data.table(
	ct = d$ct,
	pseudotime = psudotime_v
	) %>% cbind(as.matrix(t(mat[gene_v, ])))

    d_plot = d_plot[ct %in% x_celltype_order]
    d_plot$pseudotime %<>% rank
    max(d_plot$pseudotime)
    d_plot$pseudotime = d_plot$pseudotime / (max(d_plot$pseudotime, na.rm = T)) * 100

    d_plot = melt(d_plot, id.vars = c("ct", "pseudotime"))
    # d_plot[, value := scale(value), by = variable]

    p = ggplot(d_plot) + aes(x = pseudotime, y = value, color = variable) +
	geom_smooth(se = F, method = "loess", span = 0.15) +
	theme_classic()
    print(p)
}
dev.off()

## }}}

## Bivalent genes has higher expression than univalent genes and background {{{

d_plot = rbindlist(d_b_l)

ggplot(d_plot) + aes(x = ct, y = y, fill = cat) +
    geom_boxplot(outliers = F)


## }}}

## Back {{{
so_sub@meta.data$ct %<>% factor(levels = celltype_order)

pdf("./figures/247_violin_plot_b_cell_markers.pdf", width = 12, height = 4)
VlnPlot(
    so_sub,
    features = c("DNTT", "RAG1", "PAX5", "TCF3", "EBF1", "ID2"),
    group.by = "ct",
    combine = F
)
dev.off()


gene_v <- c("DNTT", "RAG1", "RAG2", "PAX5", "EBF1")
x <- so@meta.data$ct

pdf("./figures/247_scatter_plot_b_cell_markers.pdf", width = 6, height = 4)
for (gene in gene_v) {
    print(gene)
    d_plot <- data.table(
        ct = x,
        y = d[gene, ],
        gene = gene
    ) %>% na.omit()

    d_plot <- d_plot[ct %in% celltype_order]
    d_plot$ct <- factor(d_plot$ct, levels = celltype_order)

    p <- ggplot(d_plot, aes(x = ct, y = y)) +
        geom_boxplot() +
        theme_classic() +
        xlab("B cell annotation") +
        ylab(paste0(gene, " (scaled)")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
}
dev.off()

## }}}

