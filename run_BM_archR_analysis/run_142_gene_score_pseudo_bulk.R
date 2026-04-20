# =============================================================================
# Script: run_130_gene_score_pseudo_bulk.R
#
# Input:
#   - Gene score matrices (RDS files) from ./tmp/:
#       * matrix_original_gene_score_h3k27me3.rds
#       * matrix_original_gene_score_h3k4me1.rds
#       * matrix_original_gene_score_h3k4me2.rds
#       * matrix_original_gene_score_h3k4me3.rds
#       * matrix_original_gene_score_h3k4me1_cooc.rds
#       * matrix_original_gene_score_h3k4me2_cooc.rds
#       * matrix_original_gene_score_h3k4me3_cooc.rds
#   - Cell cluster annotation table: table_cell_cluster_annotation_final_round2.tsv
#
# Output:
#   - Aggregated gene marker average scores by cell type:
#       * ./tmp/table_gene_marker_average_score_by_celltype.tsv.gz
# =============================================================================

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

## Load project
m_h3k27 = read_rds("./tmp/matrix_original_gene_score_h3k27me3.rds")
m_h3k4me1 = read_rds("./tmp/matrix_original_gene_score_h3k4me1.rds")
m_h3k4me2 = read_rds("./tmp/matrix_original_gene_score_h3k4me2.rds")
m_h3k4me3 = read_rds("./tmp/matrix_original_gene_score_h3k4me3.rds")
m_h3k4me1_cooc = read_rds("./tmp/matrix_original_gene_score_h3k4me1_cooc.rds")
m_h3k4me2_cooc = read_rds("./tmp/matrix_original_gene_score_h3k4me2_cooc.rds")
m_h3k4me3_cooc = read_rds("./tmp/matrix_original_gene_score_h3k4me3_cooc.rds")

m_list = list(
    H3K27me3 = m_h3k27,
    H3K4me1 = m_h3k4me1,
    H3K4me2 = m_h3k4me2,
    H3K4me3 = m_h3k4me3,
    H3K4me1_cooc = m_h3k4me1_cooc,
    H3K4me2_cooc = m_h3k4me2_cooc,
    H3K4me3_cooc = m_h3k4me3_cooc
)

d_celltype = fread("./tmp/table_cell_cluster_annotation_final_round2.tsv")

d_list = lapply(1:length(m_list), function(i) {
    print(i)
    marker_name = names(m_list)[i]
    m = m_list[[i]] %>% as.matrix
    cell_id_m = colnames(m) %>% sub(".*#", "", .)
    cell_type = d_celltype[match(cell_id_m, cell_id), ct2]
    m = m[, !is.na(cell_type)]
    cell_type = cell_type[!is.na(cell_type)]
    split_idx <- split(seq_along(cell_type), cell_type)

    res <- sapply(split_idx, function(cols) {
	rowMeans(m[, cols, drop = FALSE])
    })

    d_res = res %>% melt
    names(d_res) = c("gene", "cell_type", "gene_score")

    d_res$marker = marker_name

    d_res
}) %>% rbindlist

write_tsv(d_list, "./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")


## {{{
# d_module = fread("./tmp/table_gene_modules_all_lineage_heatmap_renamed.tsv")[, .(name, module = new_module_name)]
# module_v = unique(d_module$module)
# Cairo::CairoPDF("./tmp/figure_gene_module_marker_score_by_celltype_violin_boxplot_var_cell.pdf", width = 12, height = 12)
# for (i in 1:length(module_v)) {
# # for (i in 1:3) {
#     module_name = module_v[i]
#     print(module_name)
#     lapply(1:length(m_list), function(j) {
# 	marker_name = names(m_list)[j]
# 	gene_v = d_module[module == module_name, name]
# 	d_sub = m_list[[marker_name]][gene_v, ] %>% as.matrix %>% melt %>% data.table
# 	names(d_sub) = c("gene", "cell_id", "gene_score")
# 	d_sub$cell_id = sub(".*#", "", d_sub$cell_id)
# 	d_sub = merge(d_sub, d_celltype[, .(cell_id, ct2)], by = "cell_id")
# 	d_sub[, gene_score := scale(gene_score), by = gene]
# 	d_plot = d_sub[, .(gene_score = mean(gene_score)), by = .(ct2, cell_id)]
# 	## remove top 5% gene score to avoid outlier for each cell type
# 	quantile_cutoff = d_plot[, quantile(gene_score, 0.95, na.rm = T), by = ct2]
# 	names(quantile_cutoff) = c("ct2", "quantile_cutoff")
# 	d_plot = merge(d_plot, quantile_cutoff, by = "ct2")
# 	d_plot = d_plot[gene_score <= quantile_cutoff]
# 	d_plot[, marker := marker_name]
# 	d_plot
#     }) %>% rbindlist -> d_plot
#     d_plot = na.omit(d_plot)
#     d_plot$ct2 = factor(d_plot$ct2, levels = ordered_ct)
#     g1 = ggplot(d_plot, aes(x = ct2, y = gene_score)) +
# 	geom_violin() +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) +
# 	facet_grid(marker ~ ., scales = "free_y")
#     g2 = ggplot(d_plot, aes(x = ct2, y = gene_score)) +
# 	geom_boxplot(outlier.size = 0.5) +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) +
# 	facet_grid(marker ~ ., scales = "free_y")
#     print(g1)
#     print(g2)
# }
# dev.off()
#
#
# module_v = unique(d_module$module)
# Cairo::CairoPDF("./tmp/figure_gene_module_marker_score_by_celltype_violin_boxplot_var_cell_no_scale.pdf", width = 12, height = 12)
# for (i in 1:length(module_v)) {
# # for (i in 1:3) {
#     module_name = module_v[i]
#     print(module_name)
#     lapply(1:length(m_list), function(j) {
# 	marker_name = names(m_list)[j]
# 	gene_v = d_module[module == module_name, name]
# 	d_sub = as.matrix(m_list[[marker_name]])[gene_v, ] %>% melt %>% data.table
# 	names(d_sub) = c("gene", "cell_id", "gene_score")
# 	d_sub$cell_id = sub(".*#", "", d_sub$cell_id)
# 	d_sub = merge(d_sub, d_celltype[, .(cell_id, ct2)], by = "cell_id")
# 	# d_sub[, gene_score := scale(gene_score), by = gene]
# 	d_plot = d_sub[, .(gene_score = mean(gene_score)), by = .(ct2, cell_id)]
# 	## remove top 5% gene score to avoid outlier for each cell type
# 	# quantile_cutoff = d_plot[, quantile(gene_score, 0.95, na.rm = T), by = ct2]
# 	# names(quantile_cutoff) = c("ct2", "quantile_cutoff")
# 	# d_plot = merge(d_plot, quantile_cutoff, by = "ct2")
# 	# d_plot = d_plot[gene_score <= quantile_cutoff]
# 	d_plot[, marker := marker_name]
# 	d_plot
#     }) %>% rbindlist -> d_plot
#     d_plot = na.omit(d_plot)
#     d_plot$ct2 = factor(d_plot$ct2, levels = ordered_ct)
#     g1 = ggplot(d_plot, aes(x = ct2, y = gene_score)) +
# 	geom_violin() +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) +
# 	facet_grid(marker ~ ., scales = "free_y")
#     g2 = ggplot(d_plot, aes(x = ct2, y = gene_score)) +
# 	geom_boxplot(outlier.size = 0.5) +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) +
# 	facet_grid(marker ~ ., scales = "free_y")
#     print(g1)
#     print(g2)
# }
# dev.off()
#
#
#
# module_v = unique(d_module$module)
# Cairo::CairoPDF("./tmp/figure_gene_module_marker_score_by_celltype_violin_boxplot_var_gene_no_scale.pdf", width = 12, height = 12)
# for (i in 1:length(module_v)) {
# # for (i in 1:3) {
#     module_name = module_v[i]
#     print(module_name)
#     lapply(1:length(m_list), function(j) {
# 	marker_name = names(m_list)[j]
# 	gene_v = d_module[module == module_name, name]
# 	d_sub = as.matrix(m_list[[marker_name]])[gene_v, ] %>% melt %>% data.table
# 	names(d_sub) = c("gene", "cell_id", "gene_score")
# 	d_sub$cell_id = sub(".*#", "", d_sub$cell_id)
# 	d_sub = merge(d_sub, d_celltype[, .(cell_id, ct2)], by = "cell_id")
# 	# d_sub[, gene_score := scale(gene_score), by = gene]
# 	d_plot = d_sub[, .(gene_score = mean(gene_score)), by = .(ct2, gene)]
# 	## remove top 5% gene score to avoid outlier for each cell type
# 	# quantile_cutoff = d_plot[, quantile(gene_score, 0.95, na.rm = T), by = ct2]
# 	# names(quantile_cutoff) = c("ct2", "quantile_cutoff")
# 	# d_plot = merge(d_plot, quantile_cutoff, by = "ct2")
# 	# d_plot = d_plot[gene_score <= quantile_cutoff]
# 	d_plot[, marker := marker_name]
# 	d_plot
#     }) %>% rbindlist -> d_plot
#     d_plot = na.omit(d_plot)
#     d_plot$ct2 = factor(d_plot$ct2, levels = ordered_ct)
#     g1 = ggplot(d_plot, aes(x = ct2, y = gene_score)) +
# 	geom_violin() +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) +
# 	facet_grid(marker ~ ., scales = "free_y")
#     g2 = ggplot(d_plot, aes(x = ct2, y = gene_score)) +
# 	geom_boxplot(outlier.size = 0.5) +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) +
# 	facet_grid(marker ~ ., scales = "free_y")
#     print(g1)
#     print(g2)
# }
# dev.off()
#
#
# ## Boxplot/violin plot merging some cell types
# d_celltype[, ct3 := revalue(ct2, c(
# 	"Pre-Pro-B cells" = "Immature B cells",
# 	"Mature B cells1" = "Mature B cells",
# 	"Mature B cells2" = "Mature B cells",
# 	"Memory B cells" = "Mature B cells",
# 	"Plasma cells" = "Mature B cells",
# 	"Myelocyte/Classical Monocytes1" = "Monocytes",
# 	"Myelocyte/Classical Monocytes2" = "Monocytes",
# 	"Myelocyte/Classical Monocytes3" = "Monocytes",
# 	"Non-classical Monocytes" = "Monocytes",
# 	"Erythroid progenitors2" = "Erythroid progenitors",
# 	"Erythroid progenitors2" = "Erythroid progenitors",
# 	"Erythroid progenitors3" = "Erythroid progenitors"
# ))]
#
# ordered_ct3 <- c(
#     "HSPC",
#     "Immature B cells",
#     "Mature B cells",
#     "Monocytes",
#     "Erythroid progenitors"
# )
#
#
# ordered_markers <- c(
# 	"H3K27me3",
# 	"H3K4me1",
# 	"H3K4me2",
# 	"H3K4me3",
# 	"H3K4me1_cooc",
# 	"H3K4me2_cooc",
# 	"H3K4me3_cooc"
# )
#
# module_v = unique(d_module$module)
# Cairo::CairoPDF("./tmp/figure_gene_module_marker_score_by_celltype_violin_boxplot_merge.pdf", width = 6, height = 4)
# # for (i in 1:length(module_v)) {
# for (i in c(11, 12, 8, 9)) {
#     module_name = module_v[i]
#     print(module_name)
#     lapply(1:length(m_list), function(j) {
# 	marker_name = names(m_list)[j]
# 	gene_v = d_module[module == module_name, name]
# 	d_sub = m_list[[marker_name]][gene_v, ] %>% as.matrix %>% melt %>% data.table
# 	names(d_sub) = c("gene", "cell_id", "gene_score")
# 	d_sub$cell_id = sub(".*#", "", d_sub$cell_id)
# 	d_sub = merge(d_sub, d_celltype[, .(cell_id, ct3)], by = "cell_id")
# 	d_sub[, gene_score := scale(gene_score), by = gene]
# 	d_plot = d_sub[, .(gene_score = mean(gene_score)), by = .(ct3, cell_id)]
# 	## remove top 5% gene score to avoid outlier for each cell type
# 	quantile_cutoff = d_plot[, quantile(gene_score, 0.95, na.rm = T), by = ct3]
# 	names(quantile_cutoff) = c("ct3", "quantile_cutoff")
# 	d_plot = merge(d_plot, quantile_cutoff, by = "ct3")
# 	d_plot = d_plot[gene_score <= quantile_cutoff]
# 	d_plot[, marker := marker_name]
# 	d_plot
#     }) %>% rbindlist -> d_plot
#     d_plot = na.omit(d_plot)
#     d_plot$ct3 = factor(d_plot$ct3, levels = ordered_ct3)
#     # g1 = ggplot(d_plot, aes(x = ct3, y = gene_score, fill = marker)) +
# 	# geom_violin() +
# 	# theme_classic() +
# 	# scale_fill_nejm() +
# 	# theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	# ggtitle(paste0("module ", module_name))
#     d_plot$marker = factor(d_plot$marker, levels = ordered_markers)
#     g2 = ggplot(d_plot, aes(x = ct3, y = gene_score, fill = marker)) +
# 	geom_boxplot(outlier.size = 0.5) +
# 	theme_classic() +
# 	scale_fill_nejm() +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	ggtitle(paste0("module ", module_name)) 
#     # print(g1)
#     print(g2)
# }
# dev.off()

## }}}
