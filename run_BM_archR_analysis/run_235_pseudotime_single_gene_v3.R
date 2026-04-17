library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

table(proj_h3k27$ct3)

d_meta_k27 = fread("./tmp/table_metadata_h3k27me3_final3_with_trajectories.tsv")

umap_h3k27 = fread("./tmp/table_umap_h3k27me3_final3.tsv")


## TODO: Use the following code to

## imputed gene score
m_h3k4me1 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me1_from_h3k27me3.rds")
m_h3k4me2 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me2_from_h3k27me3.rds")
m_h3k4me3 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me3_from_h3k27me3.rds")
m_h3k4me1_cooc = read_rds("./tmp/matrix_imputed_gene_score_h3k4me1_cooc.rds")
m_h3k4me2_cooc = read_rds("./tmp/matrix_imputed_gene_score_h3k4me2_cooc.rds")
m_h3k4me3_cooc = read_rds("./tmp/matrix_imputed_gene_score_h3k4me3_cooc.rds")
m_h3k27 = read_rds("./tmp/matrix_imputed_gene_score_h3k27me3.rds")
# m_h3k27_orig = read_rds("./tmp/matrix_original_gene_score_h3k27me3.rds")

colnames(m_h3k4me1) %<>% sub(".*#", "", .)
colnames(m_h3k4me2) %<>% sub(".*#", "", .)
colnames(m_h3k4me3) %<>% sub(".*#", "", .)
colnames(m_h3k27) %<>% sub(".*#", "", .)
# colnames(m_h3k27_orig) %<>% sub(".*#", "", .)
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
names(m_list)

gene_list1 = c("MECOM", "PRDM16", "HOXB5", "PAX5", "EBF1", "CEBPD", "CEBPA", "CEBPE", "FOSB", "JUNB", "ZFPM1", "GATA1")
gene_list2 = c( "ASCL4", "ASXL3", "ATOH1", "BMP4", "CPEB2", "ESRRB", "FEZF1", "FEZF2", "FGF4", "FGF10", "FGF13", "FOXA2", "FZD10", "GLI2", "IRX2", "IRX3", "IRX4", "LHX8", "LMO3", "MEOX2", "NKX2-1", "NKX2-8", "NKX6-1", "NKX6-2", "PRDM16", "PAX7", "RXRG", "SOX1", "SOX3", "SOX9", "SOX14", "SOX17", "TBX3", "TRIM49", "TRIM49C", "VGLL2", "VGLL3", "WNT2", "WNT7A")
gene_list3 = c(
    "OSR2", "PAX3", "HMX3", "NR2F1", "ZIC2", "SOX21", "SHOX2", "SP8", "IGF2BP3", "ID3", "EGR3", "DACT3", "HBZ", "HBA1", "HBA2", "HBG1", "HBG2", "HBD", "HBB",
    "ENAH", "WT1", "LOXL4", "VANGL2", "MECOM", "PRDM5", "HOXB5"
)

gene_list_lab_meeting = c("MECOM", "PAX5", "CEBPA", "GATA1", "SOX9", "IGF2BP3", "IGF2BP1")

color = c(
    H3K27me3 = '#3B8FC4', 
    # H3K4me1 = '#F4D03F', 
    H3K4me1 = "#8c8023",
    H3K4me1 = '#FFD000',
    H3K4me2 = '#F39C12', 
    H3K4me3 = '#E74C3C', 
    H3K4me1_cooc = '#56B870', 
    H3K4me2_cooc = '#D16BA5', 
    H3K4me3_cooc = '#9B59B6'
)

marker_ord = c("H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me1_cooc", "H3K4me2_cooc", "H3K4me3_cooc")

gene_list_l = list(
    # gene_list1,
    # gene_list2,
    # gene_list3
    gene_list_lab_meeting
)

for (i in 1:length(gene_list_l)) {
    print(i)
    gene_list = gene_list_l[[i]]
    pdf_name = str_glue("./tmp/figure_trajectory_single_gene_scores_list_bi{i}.pdf")
    Cairo::CairoPDF(pdf_name, width = 12, height = 4)
    for (gene_name in c(gene_list)) {
	print(gene_name)
	d_gene_score = lapply(1:length(m_list), function(i) {
	    m = m_list[[i]]

	    gene_score_v = log1p(m[gene_name, ])
	    # quants = quantile(gene_score_v, probs = c(0.95), na.rm = T)
	    # gene_score_v[gene_score_v > quants[1]] = quants[1]
	    data.table(
		cell_id = names(gene_score_v),
		gene_score = gene_score_v,
		# gene_score = scale(gene_score_v)[, 1],
		marker = names(m_list)[i]
	    )
	}) %>% rbindlist

	d_plot = merge(d_meta_k27, d_gene_score, by = "cell_id")
	d_plot = melt(d_plot, id.vars = c("cell_id", "Bcell_Trajectory", "Monocyte_Trajectory", "Erythroid_Trajectory", "marker"), measure.vars = "gene_score", value.name = "gene_score")
	d_plot$marker = factor(d_plot$marker, levels = marker_ord)
	# d_plot = d_plot[(grepl("H3K4me2_cooc", marker))]
	# d_plot = d_plot[marker %in% c("H3K4me3_cooc", "H3K27me3", "H3K4me3")]
	d_plot = d_plot[marker %in% c("H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3")]

	g1 = ggplot(d_plot, aes(x = Bcell_Trajectory, y = gene_score, color = marker)) +
	    geom_point(alpha = 0.1) +
	    geom_smooth(method = "loess", se = F) +
	    theme_classic() +
	    ggtitle(str_glue("B cell trajectory - {gene_name} gene score")) +
	    scale_color_manual(values = color) 

	g2 = ggplot(d_plot, aes(x = Monocyte_Trajectory, y = gene_score, color = marker)) +
	    geom_point(alpha = 0.1) +
	    geom_smooth(method = "loess", se = F) +
	    theme_classic() +
	    scale_color_manual(values = color) +
	    ggtitle(str_glue("Monocyte trajectory - {gene_name} gene score")) 

	g3 = ggplot(d_plot, aes(x = Erythroid_Trajectory, y = gene_score, color = marker)) +
	    geom_point(alpha = 0.1) +
	    geom_smooth(method = "loess", se = F) +
	    theme_classic() +
	    scale_color_manual(values = color) +
	    ggtitle(str_glue("Erythroid trajectory - {gene_name} gene score"))

	print(ggarrange(g1, g2, g3, ncol = 3, nrow = 1))
    }
    dev.off()
}


## HACK: Try to make the code better
gene_name = "PAX5"
d_gene_score = lapply(1:length(m_list), function(i) {
    m = m_list[[i]]

    gene_score_v = log1p(m[gene_name, ])
    # quants = quantile(gene_score_v, probs = c(0.95), na.rm = T)
    # gene_score_v[gene_score_v > quants[1]] = quants[1]
    data.table(
	cell_id = names(gene_score_v),
	gene_score = gene_score_v,
	# gene_score = scale(gene_score_v)[, 1],
	marker = names(m_list)[i]
    )
}) %>% rbindlist
d_plot = merge(d_meta_k27, d_gene_score, by = "cell_id")
d_plot = melt(d_plot, id.vars = c("cell_id", "Bcell_Trajectory", "Monocyte_Trajectory", "Erythroid_Trajectory", "marker"), measure.vars = "gene_score", value.name = "gene_score")
d_plot$marker = factor(d_plot$marker, levels = marker_ord)
# d_plot = d_plot[(grepl("H3K4me2_cooc", marker))]
d_plot = d_plot[marker %in% c("H3K4me2_cooc", "H3K27me3", "H3K4me2")]
d_plot = d_plot[!is.na(Bcell_Trajectory)]

d_x =dcast(d_plot[marker %in% c("H3K4me2", "H3K27me3", "H3K4me2_cooc")], cell_id + Bcell_Trajectory ~ marker, value.var = "gene_score") %>% na.omit
d_x$H3K4me2_cooc %<>% log1p
d_x$H3K4me2 %<>% log1p
d_x$H3K27me3 %<>% log1p

d_x$x = d_x$H3K4me2 * d_x$H3K27me3

lm_fit = lm(H3K4me2_cooc ~ x, data = d_x)
summary(lm_fit)
lm_fit = lm(H3K4me2_cooc ~ H3K4me2 + H3K27me3, data = d_x)
summary(lm_fit)

plot(lm_fit)

## show to predicted vs actual
d_x$predicted = predict(lm_fit, newdata = d_x)
g_pred = ggplot(d_x, aes(x = predicted, y = H3K4me2_cooc)) +
	geom_point(alpha = 0.2) +
	geom_smooth(method = "lm", se = F, color = "red") +
	theme_classic() +
	xlab("Predicted PAX5 gene score (H3K4me2_cooc)") +
	ylab("Actual PAX5 gene score (H3K4me2_cooc)") +
	ggtitle("PAX5 gene score prediction using H3K4me2 and H3K27me3")
g_pred

ggplot(d_x, aes(x = Bcell_Trajectory, y = predicted)) +
    geom_point(aes(y = H3K4me2_cooc), color = "blue", alpha = 0.2) +
	geom_point(alpha = 0.2) +
	geom_smooth(method = "loess", se = F, color = "blue") +
	theme_classic() +
	ylab("PAX5 gene score (H3K4me2_cooc)")








## multiple genes with one marker
gene_name_v = c("HMX3", "SOX21", "ZIC2", "SHOX2")
color = c(
    HMX3 =  '#0072B2',
    SOX21 = '#E69F00',
    ZIC2 = '#009E73',
    SHOX2 = '#D55E00'
)

d_gene_score = lapply(1:length(gene_name_v), function(i) {
    m = m_list[["H3K4me2_cooc"]]

    gene_name = gene_name_v[i]
    gene_score_v = m[gene_name, ]
    quants = quantile(gene_score_v, probs = c(0.95), na.rm = T)
    gene_score_v[gene_score_v > quants[1]] = quants[1]
    data.table(
	cell_id = names(gene_score_v),
	gene_score = gene_score_v,
	# gene_score = scale(gene_score_v)[, 1],
	marker = names(m_list)[i],
	gene_name = gene_name
    )
}) %>% rbindlist

d_plot = merge(d_meta_k27, d_gene_score, by = "cell_id")
d_plot = melt(d_plot, id.vars = c("cell_id", "Bcell_Trajectory", "Monocyte_Trajectory", "Erythroid_Trajectory", "gene_name"), measure.vars = "gene_score", value.name = "gene_score")
d_plot$marker = factor(d_plot$marker, levels = marker_ord)

g3 = ggplot(d_plot, aes(x = Erythroid_Trajectory, y = gene_score, color = gene_name)) +
    # geom_point(alpha = 0.1) +
    geom_smooth(method = "loess", se = F) +
    theme_classic() +
    scale_color_manual(values = color) +
    ggtitle(str_glue("Erythroid trajectory - {gene_name} gene score"))

g3
ggsave("./tmp/figure_erythroid_trajectory_multiple_genes_h3k4me2_cooc_v1.pdf", g3, width = 4, height = 4)
ggsave("./tmp/figure_erythroid_trajectory_multiple_genes_h3k4me2_cooc_v2.pdf", g3, width = 6, height = 4)
ggsave("./tmp/figure_erythroid_trajectory_multiple_genes_h3k4me2_cooc_v3.pdf", g3, width = 8, height = 4)
ggsave("./tmp/figure_erythroid_trajectory_multiple_genes_h3k4me2_cooc_v4.pdf", g3, width = 12, height = 4)

## gene module scores in HSPC cluster
m = m_list[["H3K4me2_cooc"]]
cell_v = d_meta_k27[ct2 == "HSPC", cell_id]
m = m[, colnames(m) %in% cell_v]
gene_module_score_v = rowMeans(m)



plot_trajectory_gene_score = function(gene_name = "MECOM", m_list) {
    lapply(1:length(m_list), function(i) {
	m = m_list[[i]]
	gene_score_v = m[gene_name, ]

	## ceil the top by 10 percent by gene
	quants = quantile(gene_score_v, probs = c(0.95), na.rm = T)
	gene_score_v[gene_score_v > quants[1]] = quants[1]

	d_plot = d_meta_k27[cell_id %in% names(gene_score_v)]
	d_plot$gene_score = gene_score_v[match(d_plot$cell_id, names(gene_score_v))]

	ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = gene_score)) +
	    geom_point(size = 0.5) +
	    scale_color_gradient(low = "lightgrey", high = "red") +
	    theme_classic() +
	    theme(
		legend.position = "right",
		axis.title = element_blank()
		) +
	    ggtitle(str_glue("{names(m_list)[i]} - {gene_name}"))
    }) %>% ggarrange(plotlist = ., ncol = 2, nrow = 4)
}

color = c(
    H3K27me3 = '#3B8FC4', 
    # H3K4me1 = '#F4D03F', 
    H3K4me1 = "#8c8023",
    H3K4me1 = '#FFD000',
    H3K4me2 = '#F39C12', 
    H3K4me3 = '#E74C3C', 
    H3K4me1_cooc = '#56B870', 
    H3K4me2_cooc = '#D16BA5', 
    H3K4me3_cooc = '#9B59B6'
)

marker_ord = c("H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me1_cooc", "H3K4me2_cooc", "H3K4me3_cooc")


Cairo::CairoPDF("./tmp/figure_trajectory_gene_module_scores.pdf", width = 16, height = 12)
d_pseudotime = d_meta_k27[, .(cell_id, pseudotime = Bcell_Trajectory)] %>% na.omit
d_plot = merge(d_module_score, d_pseudotime, by = "cell_id")
d_plot$marker = factor(d_plot$marker, levels = marker_ord)
g = ggplot(d_plot, aes(x = pseudotime, y = module_score, color = marker)) +
    geom_smooth(method = "loess", se = F) +
    facet_wrap(~gene_module, scales = "free_y") +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.title = element_blank()
	) +
    ggtitle("B cell trajectory - Gene module scores")
print(g)

d_pseudotime = d_meta_k27[, .(cell_id, pseudotime = Monocyte_Trajectory)] %>% na.omit
d_plot = merge(d_module_score, d_pseudotime, by = "cell_id")
d_plot$marker = factor(d_plot$marker, levels = marker_ord)
g = ggplot(d_plot, aes(x = pseudotime, y = module_score, color = marker)) +
    geom_smooth(method = "loess", se = F) +
    facet_wrap(~gene_module, scales = "free_y") +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.title = element_blank()
	) +
    ggtitle("Monocyte trajectory - Gene module scores")
print(g)

d_pseudotime = d_meta_k27[, .(cell_id, pseudotime = Erythroid_Trajectory)] %>% na.omit
d_plot = merge(d_module_score, d_pseudotime, by = "cell_id")
d_plot$marker = factor(d_plot$marker, levels = marker_ord)d
g = ggplot(d_plot, aes(x = pseudotime, y = module_score, color = marker)) +
    geom_smooth(method = "loess", se = T) +
    facet_wrap(~gene_module, scales = "free_y") +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.title = element_blank()
	) +
    ggtitle("Erythroid trajectory - Gene module scores")
print(g)
dev.off()
