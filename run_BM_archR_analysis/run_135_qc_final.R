################################################################################
# Script: run_135_qc_final.R
# Purpose: Summarize final QC metrics across single-mark and co-occurrence ArchR
#          projects and generate fragment-count QC figures by cell type, marker,
#          donor, and antibody species
#
# INPUT:
#   - ArchR project directories:
#     * ./tmp/H3K27me3_coCnT_final3/
#     * ./tmp/H3K4me1_coCnT_final3/
#     * ./tmp/H3K4me2_coCnT_final3/
#     * ./tmp/H3K4me3_coCnT_final3/
#     * ./tmp/H3K4me1_coCnT_cooc_final3/
#     * ./tmp/H3K4me2_coCnT_cooc_final3/
#     * ./tmp/H3K4me3_coCnT_cooc_final3/
#
# OUTPUT:
#   - QC summary table:
#     * ./tmp/table_qc_final_all.tsv.gz
#   - PDF figures:
#     * ./figures/135_figure_nFrags_by_cell_type.pdf
#     * ./figures/135_figure_nFrags_on_umap_by_marker.pdf
#     * ./figures/135_figure_qc_umap_by_donnor.pdf
#     * ./figures/135_figure_qc_antibody_species_by_cluster.pdf
#     * ./figures/135_figure_qc_fragments_by_cell_type_by_donnor.pdf
#     * ./figures/135_figure_qc_fragments_by_cell_type_by_donnor_by_antibody_species.pdf
################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggpubr)
library(ggsci)

library(ArchR)

addArchRGenome("hg38")
addArchRThreads(threads = 8)

ordered_ct <- c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3",
    "NK cells1",
    "NK cells2",
    "T cells",
    "Plasmacytoid dendritic cells",
    "Non-classical Monocytes",
    "EoBaso progenitors"
)

proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final3/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final3/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final3/")
proj_h3k4me1_cooc = loadArchRProject("./tmp/H3K4me1_coCnT_cooc_final3/")
proj_h3k4me2_cooc = loadArchRProject("./tmp/H3K4me2_coCnT_cooc_final3/")
proj_h3K4me3_cooc = loadArchRProject("./tmp/H3K4me3_coCnT_cooc_final3/")

d_k27 = as.data.table(proj_k27@cellColData)
d_h3k4me1 = as.data.table(proj_h3k4me1@cellColData)
d_h3k4me2 = as.data.table(proj_h3k4me2@cellColData)
d_h3k4me3 = as.data.table(proj_h3k4me3@cellColData)
d_h3k4me1_cooc = as.data.table(proj_h3k4me1_cooc@cellColData)
d_h3k4me2_cooc = as.data.table(proj_h3k4me2_cooc@cellColData)
d_h3k4me3_cooc = as.data.table(proj_h3K4me3_cooc@cellColData)


d_k27$marker = "H3K27me3"
d_h3k4me1$marker = "H3K4me1"
d_h3k4me2$marker = "H3K4me2"
d_h3k4me3$marker = "H3K4me3"
d_h3k4me1_cooc$marker = "H3K27me3_H3K4me1"
d_h3k4me2_cooc$marker = "H3K27me3_H3K4me2"
d_h3k4me3_cooc$marker = "H3K27me3_H3K4me3"

d_all = rbind(
    d_k27[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)],
    d_h3k4me1[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)],
    d_h3k4me2[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)],
    d_h3k4me3[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)],
    d_h3k4me1_cooc[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)],
    d_h3k4me2_cooc[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)],
    d_h3k4me3_cooc[, .(cell_id, nFrags, library_id, antibody_id, ct2, marker)]
)

write_tsv(d_all, "./tmp/table_qc_final_all.tsv.gz")

d_all = fread("./tmp/table_qc_final_all.tsv.gz")

d_k27$ct2 %>% table

## get UMAP
umap_k27 = as.data.table( getEmbedding( ArchRProj = proj_k27, embedding = "UMAP_Harmony"))
names(umap_k27) = c("UMAP1", "UMAP2")
umap_k27$cell_id = d_k27$cell_id

## nFrags per cell type {{{
ordered_ct <- c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3",
    "NK cells1",
    "NK cells2",
    "T cells",
    "Plasmacytoid dendritic cells",
    "Non-classical Monocytes",
    "EoBaso progenitors"
)

ordered_ct <- c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3",
    "Non-classical Monocytes",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3",
    "NK cells1",
    "NK cells2",
    "T cells",
    "Plasmacytoid dendritic cells",
    "EoBaso progenitors"
)

annot_ct2 = c(
    "C1" = "NK cells1",
    "C2" = "T cells",
    "C3" = "Mature B cells1",
    "C4" = "Plasma cells",
    "C5" = "Mature B cells2",
    "C6" = "Memory B cells",
    "C7" = "Erythroid progenitors3",
    "C8" = "Erythroid progenitors2",
    "C9" = "NK cells2",
    "C10" = "HSPC",
    "C11" = "Plasmacytoid dendritic cells",
    "C12" = "Myelocyte/Classical Monocytes1",
    "C13" = "Myelocyte/Classical Monocytes2",
    "C14" = "Myelocyte/Classical Monocytes3",
    "C15" = "Non-classical Monocytes",
    "C16" = "Pre-Pro-B cells",
    "C17" = "Immature B cells",
    "C18" = "Erythroid progenitors1",
    "C19" = "EoBaso progenitors"
)



color <- c(
    C1  = "#8C8178",
    C2  = "#A3623A",
    C3  = "#B4009E",
    C4  = "#7D0033",
    C5  = "#89288F",
    C6  = "#B6B0FF",
    C7  = "#499195",
    C8  = "#208A42",
    C9  = "#4A5B5D",
    C10 = "#D51F26",
    C11 = "#FFBBED",
    C12 = "#8B3E1D",
    C13 = "#F47D2B",
    C14 = "#E09E5A",
    C15 = "#BDA687",
    C16 = "#6C2B85",
    C17 = "#8250C4",
    C18 = "#14532D",
    C19 = "#8AD4EB"
)


color_ct2 = color
names(color_ct2) = annot_ct2


d_all$ct2 %<>% factor(., levels = ordered_ct)

d_all$ct2 %>% unique

title_x = "Number of fragments per cell by cell type"
g3 = ggplot(d_all[!is.na(ct2)], aes(x = ct2, y = nFrags, fill = ct2)) +
    geom_boxplot(outliers = F, position = position_dodge(width = 0.8)) +
    facet_wrap(~marker, ncol = 4, scales = "free_y") +
    stat_compare_means(ref.group = "HSPC", hide.ns=T, label = "p.signif") +
    scale_fill_manual(values = color_ct2) +
    scale_y_log10() +
    theme_classic() +
    theme(
	axis.text.x = element_text(angle = 90, hjust = 1)
	) +
    labs(y = "Number of fragments (log10 scale)", x = "Cell Type", title = title_x)
g3
ggsave("./figures/135_figure_nFrags_by_cell_type.pdf", g3, width = 28)



## }}}

## nFrags in UMAP {{{
marker_v = unique(d_all$marker)

color = c(
    H3K27me3 = '#3B8FC4', 
    H3K4me1 = '#FFD000',
    H3K4me2 = '#F39C12', 
    H3K4me3 = '#E74C3C',
    H3K4me1_cooc = '#56B870', 
    H3K4me2_cooc = '#D16BA5', 
    H3K4me3_cooc = '#9B59B6'
)

pdf("./figures/135_figure_nFrags_on_umap_by_marker.pdf", width = 8, height = 8)
for (i in 1:length(color)) {
    print(i)
    marker_x = marker_v[i]
    color_x = color[i]

    d_plot = d_all[marker == marker_x]
    d_plot = merge(d_plot, umap_k27, by = "cell_id")
    p = ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = log10(nFrags))) +
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_gradient(low = "white", high = color_x) +
	theme_classic() +
	labs(
	    title = paste0(marker_x, ": nFrags on UMAP (H3K27me3 UMAP)"),
	    x = "UMAP1",
	    y = "UMAP2",
	    color = "nFrags (log10)"
	)
    print(p)
}
dev.off()
## }}}


## H3K27me3 UMAP by donnor
p = ggplot(merge(d_k27[assignment %in% c("0", "1", "2", "3")], umap_k27, by = "cell_id"), aes(x = UMAP1, y = UMAP2, color = assignment)) +
    geom_point(size = 0.1, alpha = 0.5) +
    theme_classic() +
    facet_wrap(~assignment) +
    labs(
	title = "H3K27me3 UMAP by library",
	x = "UMAP1",
	y = "UMAP2",
	color = "Library ID"
    ) +
    scale_color_nejm() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))
ggsave("./figures/135_figure_qc_umap_by_donnor.pdf", p, width = 8, height = 8)


## H3K27me3 stack bar plot for antibody species by cluster
d_k27$Clusters_plus = factor(d_k27$Clusters_plus, levels = paste0("C", 1:length(unique(d_k27$Clusters_plus))))
p = ggplot(d_k27, aes(x = Clusters_plus, fill = as.factor(antibody_species))) +
	geom_bar(position = "fill") +
	scale_fill_nejm() +
	theme_classic() +
	labs(
	title = "H3K27me3 antibody species by cluster",
	x = "Cluster",
	y = "Proportion",
	fill = "Antibody Species"
	) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	guides(fill = guide_legend(override.aes = list(size=2, alpha = 1)))

ggsave("./figures/135_figure_qc_antibody_species_by_cluster.pdf", p, width = 6, height = 4)

## Read per cell per donor per cell type boxplot
d_plot = d_k27[assignment %in% c("0", "1", "2", "3")]
d_plot = d_plot[, .(nFrags = median(nFrags)), by = .(assignment, ct2)]
d_plot$ct2 %<>% factor(., levels = ordered_ct)
p = ggplot(d_plot, aes(x = ct2, y = nFrags)) +
    geom_boxplot(outliers = F) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, aes(color = assignment)) +
    ## y-axis begin at 0
    expand_limits(y = 0) +
    scale_color_nejm() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("./figures/135_figure_qc_fragments_by_cell_type_by_donnor.pdf", p, width = 8, height = 4)

## Read per cell per donor and cell type boxplot with facet by antibody species
d_plot = d_k27[assignment %in% c("0", "1", "2", "3")]
d_plot = d_plot[, .(nFrags = median(nFrags)), by = .(assignment, ct2, antibody_species)]
d_plot$ct2 %<>% factor(., levels = ordered_ct)
p = ggplot(d_plot, aes(x = ct2, y = nFrags)) +
	geom_boxplot(outliers = F) +
	geom_jitter(width = 0.2, size = 0.5, alpha = 1, aes(color = assignment)) +
	# geom_line(aes(group = assignment), alpha = 0.5, size = 0.5) +
	## y-axis begin at 0
	expand_limits(y = 0) +
	scale_color_nejm() +
	facet_grid(antibody_species ~ ., scales = "free_y") +
	stat_compare_means(ref.group = "HSPC", hide.ns=T, label = "p.signif", method = "t.test", paired = T, label.y = max(d_plot$nFrags) * 1.1) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./figures/135_figure_qc_fragments_by_cell_type_by_donnor_by_antibody_species.pdf", p, width = 5, height = 4)

# ## H3K27me3 Ct2 by donnor percentage stack bar plot
# d_ct2 = d_k27[, .N, by = .(assignment, ct2)]
# d_ct2[, total := sum(N), by = .(assignment)]
# d_ct2[, prop := N / total]
# d_ct2$ct2 %<>% factor(., levels = ordered_ct)
# p = ggplot(d_ct2[assignment %in% c("0", "1", "2", "3")], aes(x = ct2, y = prop, fill = assignment)) +
# 	geom_bar(stat = "identity", position = "fill") +
# 	scale_fill_nejm() +
# 	theme_classic() +
# 	labs(
# 	title = "H3K27me3 Ct2 by library",
# 	x = "Library ID",
# 	y = "Proportion",
# 	fill = "Cell Type"
# 	) +
# 	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 	guides(fill = guide_legend(override.aes = list(size=2, alpha = 1)))
# p

## Show nFrags on UMAP
# d_plot = merge(d_all, umap_k27, by = "cell_id", all.x = TRUE)
# p_list = lapply(
# 	unique(d_plot$marker),
# 	function(x) {
# 	d_sub = d_plot[marker == x]
# 	p = ggplot(d_sub, aes(x = UMAP1, y = UMAP2, color = log10(nFrags))) +
# 	    geom_point(size = 0.5, alpha = 0.5) +
# 	    scale_color_viridis_c() +
# 	    theme_classic() +
# 	    labs(
# 		title = paste0(x, ": nFrags on UMAP (H3K27me3 UMAP)"),
# 		x = "UMAP1",
# 		y = "UMAP2",
# 		color = "nFrags"
# 	    )
# 	return(p)
# 	}
# )
#
# ggarrange(plotlist = p_list, ncol = 2, nrow = 2)

## show antibody on UMAP
# d_plot
# p_list = lapply(
# 	unique(d_plot$marker),
# 	function(x) {
# 	d_sub = d_plot[marker == x]
# 	p = ggplot(d_sub, aes(x = UMAP1, y = UMAP2, color = antibody_species)) +
# 	    geom_point(size = 0.5, alpha = 0.5) +
# 	    theme_classic() +
# 	    labs(
# 		title = paste0(x, ": antibody on UMAP (H3K27me3 UMAP)"),
# 		x = "UMAP1",
# 		y = "UMAP2",
# 		color = "antibody_id"
# 	    )
# 	return(p)
# 	}
# )
# ggarrange(plotlist = p_list, ncol = 2, nrow = 2)
#
# ## show antibody species with boxplot by cluster
# p_list = lapply(
# 	unique(d_plot$marker),
# 	function(x) {
# 	d_sub = d_plot[marker == x]
# 	p = ggplot(d_sub, aes(x = ct2, fill = as.factor(antibody_species))) +
# 	    geom_bar(position = "fill") +
# 	    scale_fill_nejm() +
# 	    theme_classic() +
# 	    labs(
# 		title = paste0(x, ": antibody species by cluster"),
# 		x = "Cluster",
# 		y = "Proportion",
# 		fill = "antibody species"
# 	    ) +
# 	    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 	    guides(fill = guide_legend(override.aes = list(size=2, alpha = 1)))
# 	return(p)
# 	}
# )
# ggarrange(plotlist = p_list, ncol = 1, nrow = 1)
#
#
#
#
# ## read per cell per batch per antibody per cell type
# title_x = "Number of fragments per cell by library"
# g1 = ggplot(d_all, aes(x = library_id, y = nFrags)) +
#     geom_boxplot(outliers = F, position = position_dodge(width = 0.8)) +
#     facet_wrap(~marker, nrow = 1) +
#     scale_y_log10() +
#     theme_bw() +
#     theme(
# 	axis.text.x = element_text(angle = 90, hjust = 1)
# 	) +
#     labs(y = "Number of fragments (log10 scale)", x = "Library ID", title = title_x)
# g1
#
# title_x = "Number of fragments per cell by antibody"
# g2 = ggplot(d_all, aes(x = antibody_id, y = nFrags)) +
#     geom_boxplot(outliers = F, position = position_dodge(width = 0.8)) +
#     facet_wrap(~marker, nrow = 1) +
#     scale_y_log10() +
#     theme_bw() +
#     theme(
# 	axis.text.x = element_text(angle = 90, hjust = 1)
# 	) +
#     labs(y = "Number of fragments (log10 scale)", x = "Antibody ID", title = title_x)
# g2
#
# ## Number of cells per libray per antibody per marker
# d_count = d_all[, .N, by = .(library_id, antibody_id, marker)]
# title_x = "Number of cells per library"
# g4 = ggplot(d_count, aes(x = library_id, y = N)) +
# 	geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
# 	facet_wrap(~marker, nrow = 1) +
# 	theme_bw() +
# 	theme(
# 	    axis.text.x = element_text(angle = 90, hjust = 1)
# 	) +
# 	labs(y = "Number of cells", x = "Library ID", title = title_x)
# print(g4)
#
# title_x = "Number of cells per antibody"
# g5 = ggplot(d_count, aes(x = antibody_id, y = N)) +
# 	geom_bar(stat = "identity", position = position_dodge(width = 0.8
# )) +
# 	facet_wrap(~marker, nrow = 1) +
# 	theme_bw() +
# 	theme(
# 	    axis.text.x = element_text(angle = 90, hjust = 1)
# 	) +
# 	labs(y = "Number of cells", x = "Antibody ID", title = title_x)
# print(g5)
#
# d_count = d_all[!is.na(ct2), .N, by = .(ct2, marker)]
# title_x = "Number of cells per cell type"
# d_count$ct2 %<>% factor(., levels = ordered_ct)
# g6 = ggplot(d_count[!is.na(ct2)], aes(x = ct2, y = N)) +
# 	geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
# 	facet_wrap(~marker, nrow = 1) +
# 	theme_bw() +
# 	theme(
# 	    axis.text.x = element_text(angle = 90, hjust = 1)
# 	) +
# 	labs(y = "Number of cells", x = "Cell Type", title = title_x)
# print(g6)
