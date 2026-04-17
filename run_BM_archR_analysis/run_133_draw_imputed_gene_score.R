################################################################################
# Script: run_133_draw_imputed_gene_score.R
# Purpose: Draw imputed gene-score visualizations across histone mark datasets
#          using the H3K27me3 UMAP embedding and cell-type annotations
#
# INPUT:
#   - ArchR project directory:
#     * ./tmp/H3K27me3_coCnT_final3/
#   - RDS files (imputed gene-score matrices):
#     * ./tmp/matrix_imputed_gene_score_h3k4me1_from_h3k27me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me2_from_h3k27me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me3_from_h3k27me3.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me1_cooc.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me2_cooc.rds
#     * ./tmp/matrix_imputed_gene_score_h3k4me3_cooc.rds
#     * ./tmp/matrix_imputed_gene_score_h3k27me3.rds
#
# OUTPUT:
#   - RDS file:
#     * ./tmp/list_matrix_imputed_gene_score.rds
#   - PDF figures:
#     * ./figures/133_figure_imputed_umap_gene_score_example.pdf
#     * ./figures/133_figure_imputed_boxplot_gene_score_example.pdf
################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Cairo)

library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 11)

set.seed(123)

color = c(
    H3K27me3 = '#3B8FC4', 
    # H3K4me1 = '#F4D03F', 
    # H3K4me1 = "#8c8023",
    H3K4me1 = '#FFD000',
    H3K4me2 = '#F39C12', 
    H3K4me3 = '#E74C3C',
    H3K4me1_cooc = '#56B870', 
    H3K4me2_cooc = '#D16BA5', 
    H3K4me3_cooc = '#9B59B6'
)



proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")

d_meta_k27 = as.data.table(proj_k27@cellColData)
umap_k27 = as.data.table( getEmbedding( ArchRProj = proj_k27, embedding = "UMAP_Harmony"))
names(umap_k27) = c("UMAP1", "UMAP2")
d_meta = cbind(d_meta_k27, umap_k27)

## Update the Stem cell subcluster
d_meta$ct3 = d_meta$ct2
# d_meta[Clusters_HSPC != "Other_cells", ct3 := Clusters_HSPC]
d_meta

## imputed gene score
m_h3k4me1 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me1_from_h3k27me3.rds")
m_h3k4me2 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me2_from_h3k27me3.rds")
m_h3k4me3 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me3_from_h3k27me3.rds")
m_h3k4me1_cooc = read_rds("./tmp/matrix_imputed_gene_score_h3k4me1_cooc.rds")
m_h3k4me2_cooc = read_rds("./tmp/matrix_imputed_gene_score_h3k4me2_cooc.rds")
m_h3k4me3_cooc = read_rds("./tmp/matrix_imputed_gene_score_h3k4me3_cooc.rds")
m_h3k27 = read_rds("./tmp/matrix_imputed_gene_score_h3k27me3.rds")

colnames(m_h3k4me1) %<>% sub(".*#", "", .)
colnames(m_h3k4me2) %<>% sub(".*#", "", .)
colnames(m_h3k4me3) %<>% sub(".*#", "", .)
colnames(m_h3k27) %<>% sub(".*#", "", .)
colnames(m_h3k4me1_cooc) %<>% sub(".*#", "", .)
colnames(m_h3k4me2_cooc) %<>% sub(".*#", "", .)
colnames(m_h3k4me3_cooc) %<>% sub(".*#", "", .)

m_list = list(
    m_h3k4me1 = m_h3k4me1,
    m_h3k4me2 = m_h3k4me2,
    m_h3k4me3 = m_h3k4me3,
    m_h3k27 = m_h3k27,
    m_h3k4me1_cooc = m_h3k4me1_cooc,
    m_h3k4me2_cooc = m_h3k4me2_cooc,
    m_h3k4me3_cooc = m_h3k4me3_cooc
)

names(m_list) = c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27me3", "H3K4me1_cooc", "H3K4me2_cooc", "H3K4me3_cooc")
write_rds(m_list, "./tmp/list_matrix_imputed_gene_score.rds")

m_list = read_rds("./tmp/list_matrix_imputed_gene_score.rds")

plot_umap_gene_score = function(gene_name = "MECOM", m_list, color_v) {
    lapply(1:length(m_list), function(i) {
	color_i = color_v[names(m_list)[i]]
	m = m_list[[i]]
	gene_score_v = m[gene_name, ]

	## ceil the top by 10 percent by gene
	quants = quantile(gene_score_v, probs = c(0.95), na.rm = T)
	gene_score_v[gene_score_v > quants[1]] = quants[1]

	d_plot = d_meta[cell_id %in% names(gene_score_v)]
	d_plot$gene_score = gene_score_v[match(d_plot$cell_id, names(gene_score_v))]

	ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = gene_score)) +
	    geom_point(size = 0.2) +
	    scale_color_gradient(low = "white", high = color_i) +
	    theme_classic() +
	    theme(
		legend.position = "right",
		axis.title = element_blank()
		) +
	    ggtitle(str_glue("{names(m_list)[i]} - {gene_name}"))
    }) %>% ggarrange(plotlist = ., ncol = 2, nrow = 4)
}

Cairo::CairoPDF("./figures/133_figure_imputed_umap_gene_score_example.pdf", width = 9, height = 14)
gene_list = c("MECOM", "PAX5", "CEBPA", "GATA1", "CD34")
for (gene_name in gene_list) {
	print(gene_name)
	print(plot_umap_gene_score(gene_name, m_list, color) )
}
dev.off()

plot_boxplot_gene_score = function(gene_name, m_list, color_v) {
	lapply(1:length(m_list), function(i) {
	color_i = color_v[names(m_list)[i]]
	m = m_list[[i]]
	gene_score_v = m[gene_name, ]

	d_plot = d_meta[cell_id %in% names(gene_score_v)]
	d_plot$gene_score = gene_score_v[match(d_plot$cell_id, names(gene_score_v))]

	ggplot(d_plot, aes(x = ct3, y = gene_score, fill = ct3)) +
	    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, outliers = F) +
	    scale_fill_manual(values = rep(color_i, length(unique(d_plot$ct3)))) +
	    theme_classic() +
	    theme(
		legend.position = "none",
		axis.title.x = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1)
		) +
	    ggtitle(str_glue("{names(m_list)[i]} - {gene_name}"))
	}) %>% ggarrange(plotlist = ., ncol = 2, nrow = 4)
}

Cairo::CairoPDF("./figures/133_figure_imputed_boxplot_gene_score_example.pdf", width = 9, height = 14)
gene_list = c("MECOM", "PAX5", "CEBPA", "GATA1", "CD34", "SPI1", "RUNX1", "KLF1", "EGR1", "IRF8")
for (gene_name in gene_list) {
	print(gene_name)
	print(plot_boxplot_gene_score(gene_name, m_list, color) )
}
dev.off()


# Cairo::CairoPDF("./tmp/figure_imputed_umap_gene_score_list2.pdf", width = 9, height = 14) {{{
# gene_list1 = c("MECOM", "PRDM16", "HOXB5", "PAX5", "EBF1", "CEBPD", "CEBPA", "CEBPE", "FOSB", "JUNB", "ZFPM1", "GATA1")
# gene_list2 = c( "ASCL4", "ASXL3", "ATOH1", "BMP4", "CPEB2", "ESRRB", "FEZF1", "FEZF2", "FGF4", "FGF10", "FGF13", "FOXA2", "FZD10", "GLI2", "IRX2", "IRX3", "IRX4", "LHX8", "LMO3", "MEOX2", "NKX2-1", "NKX2-8", "NKX6-1", "NKX6-2", "PRDM16", "PAX7", "RXRG", "SOX1", "SOX3", "SOX9", "SOX14", "SOX17", "TBX3", "TRIM49", "TRIM49C", "VGLL2", "VGLL3", "WNT2", "WNT7A")
# for (gene_name in gene_list2) {
#     print(gene_name)
#     print(plot_umap_gene_score(gene_name, m_list, color) )
# }
# dev.off()
#
# Cairo::CairoPDF("./tmp/figure_imputed_umap_gene_score_list1.pdf", width = 9, height = 14)
# gene_list1 = c("MECOM", "PRDM16", "HOXB5", "PAX5", "EBF1", "CEBPD", "CEBPA", "CEBPE", "FOSB", "JUNB", "ZFPM1", "GATA1")
# gene_list2 = c( "ASCL4", "ASXL3", "ATOH1", "BMP4", "CPEB2", "ESRRB", "FEZF1", "FEZF2", "FGF4", "FGF10", "FGF13", "FOXA2", "FZD10", "GLI2", "IRX2", "IRX3", "IRX4", "LHX8", "LMO3", "MEOX2", "NKX2-1", "NKX2-8", "NKX6-1", "NKX6-2", "PRDM16", "PAX7", "RXRG", "SOX1", "SOX3", "SOX9", "SOX14", "SOX17", "TBX3", "TRIM49", "TRIM49C", "VGLL2", "VGLL3", "WNT2", "WNT7A")
# for (gene_name in gene_list1) {
#     print(gene_name)
#     print(plot_umap_gene_score(gene_name, m_list, color) )
# }
# dev.off()
#
# Cairo::CairoPDF("./tmp/figure_imputed_umap_gene_score_heatmap_annotation.pdf", width = 9, height = 14)
# gene_list1 = c("MECOM", "PRDM16", "HOXB5", "PAX5", "EBF1", "CEBPD", "CEBPA", "CEBPE", "FOSB", "JUNB", "ZFPM1", "GATA1")
# gene_list2 = c( "ASCL4", "ASXL3", "ATOH1", "BMP4", "CPEB2", "ESRRB", "FEZF1", "FEZF2", "FGF4", "FGF10", "FGF13", "FOXA2", "FZD10", "GLI2", "IRX2", "IRX3", "IRX4", "LHX8", "LMO3", "MEOX2", "NKX2-1", "NKX2-8", "NKX6-1", "NKX6-2", "PRDM16", "PAX7", "RXRG", "SOX1", "SOX3", "SOX9", "SOX14", "SOX17", "TBX3", "TRIM49", "TRIM49C", "VGLL2", "VGLL3", "WNT2", "WNT7A")
# for (gene_name in gene_list1) {
#     print(gene_name)
#     print(plot_umap_gene_score(gene_name, m_list, color) )
# }
# dev.off()
#
#
#
# Cairo::CairoPDF("./tmp/figure_imputed_umap_gene_score_list3.pdf", width = 9, height = 14)
# gene_list3 = c(
#   "OSR2", "PAX3", "HMX3", "NR2F1", "ZIC2", "SOX21", "SHOX2", "SP8", "IGF2BP3", "ID3", "EGR3", "DACT3", "HBZ", "HBA1", "HBA2", "HBG1", "HBG2", "HBD", "HBB",
#   "ENAH", "WT1", "LOXL4", "VANGL2", "MECOM", "PRDM5", "HOXB5"
# )
# for (gene_name in gene_list3) {
#     print(gene_name)
#     print(plot_umap_gene_score(gene_name, m_list, color) )
# }
# dev.off()
#
# Cairo::CairoPDF("./tmp/figure_imputed_umap_gene_score_list4.pdf", width = 9, height = 14)
# gene_list4 = c("HMX3", "SOX21", "ZIC2", "SHOX2")
# for (gene_name in gene_list4) {
#     print(gene_name)
#     print(plot_umap_gene_score(gene_name, m_list, color) )
# }
# dev.off()
#
#
#
# # m_list = list(
# #     m_h3k4me1 = m_h3k4me1,
# #     m_h3k4me2 = m_h3k4me2,
# #     m_h3k4me3 = m_h3k4me3,
# #     m_h3k27 = m_h3k27_orig
# # )
# # names(m_list) = c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27me3")
# # Cairo::CairoPDF("./tmp/figure_original_umap_gene_score.pdf", width = 12, height = 12)
# # for (gene_name in gene_list) {
# #     print(gene_name)
# #     print(plot_umap_gene_score(gene_name, m_list) )
# # }
# # dev.off()
#
#
# ## Plot gene modules in UMAP
# plot_umap_gene_score_n = function(gene_name, m_list, gene_module, color_v) {
#     lapply(1:length(m_list), function(i) {
# 	color_i = color_v[names(m_list)[i]]
# 	m = m_list[[i]]
# 	gene_score_m = m[gene_name, ]
#
# 	gene_score_v = colMeans(gene_score_m, na.rm = T)
# 	d_plot = d_meta[cell_id %in% names(gene_score_v)]
# 	d_plot$gene_score = gene_score_v[match(d_plot$cell_id, names(gene_score_v))]
#
# 	ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = gene_score)) +
# 	    geom_point(size = 0.5) +
# 	    scale_color_gradient(low = "white", high = color_i) +
# 	    theme_classic() +
# 	    theme(
# 		legend.position = "right",
# 		axis.title = element_blank()
# 		) +
# 	    ggtitle(str_glue("{names(m_list)[i]} - {gene_module}"))
#     }) %>% ggarrange(plotlist = ., ncol = 2, nrow = 4)
# }
#
# d_gene_module = fread("./tmp/table_gene_modules_all_lineage_heatmap_renamed.tsv")
# gene_module_v = unique(d_gene_module$new_module_name)
#
# # httpgd::hgd(port = 4322)
#
# Cairo::CairoPDF("./tmp/figure_imputed_umap_gene_module_score_scale.pdf", width = 9, height = 14)
# for (i in 1:length(gene_module_v)) {
#     gene_module_i = gene_module_v[i]
#     gene_name = d_gene_module[new_module_name == gene_module_i, name]
#     print(gene_module_i)
#     print(plot_umap_gene_score_n(gene_name, m_list_scale, paste1("Module ", gene_module_i), color))
# }
# dev.off() }}}
