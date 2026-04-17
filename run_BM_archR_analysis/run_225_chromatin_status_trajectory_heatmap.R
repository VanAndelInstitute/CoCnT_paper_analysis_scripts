# =============================================================================
# Script: run_224_chromatin_status_trajectory_status.R
#
# Input files:
#   - ./tmp/table_gene_chromatin_status_wide_me1.tsv
#   - ./tmp/table_gene_chromatin_status_wide_me2.tsv
#   - ./tmp/table_gene_chromatin_status_wide_me3.tsv
#
# Output file:
#   - table_gene_chromatin_status_trajectory_remove_dup_wide.tsv
#   - table_gene_chromatin_status_trajectory_remove_dup_long.tsv
#   - table_gene_chromatin_status_trajectory_remove_dup_summary.tsv
#   - table_gene_chromatin_status_trajectory_remove_dup_summary_wide.tsv
#   - table_shared_genes_between_markers_in_trajectory.tsv
#   - table_shared_genes_between_markers_in_trajectory_full_trans.tsv
#
# =============================================================================

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(UpSetR)
library(ggsci)

httpgd::hgd(port = 4322)

## Cell type order
ct_v_B = c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
)

ct_v_E = c(
    "HSPC",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3"
)

ct_v_M = c(
    "HSPC",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3"
)


## get step wise trajectory for each gene
process_trajectory <- function(d, ct_v, suffix) {
    d_b = d[, c("gene", ct_v), with = FALSE]
    trans = apply(d_b, 1, function(x) {
	x = x[x != "Other"]
	paste(rle(x[-1])$value, collapse = "->")
	# paste(rle(x[2:3])$value, collapse = "->")
    })
    d_plot = data.table(gene = d_b$gene, trans)

    d_plot[, group := suffix]
    d_plot
}

## read data
d1 = fread("./tmp/table_gene_chromatin_status_wide_me1.tsv")
d2 = fread("./tmp/table_gene_chromatin_status_wide_me2.tsv")
d3 = fread("./tmp/table_gene_chromatin_status_wide_me3.tsv")


res = list(
    process_trajectory(d1, ct_v_B, "B_me1"),
    process_trajectory(d1, ct_v_E, "E_me1"),
    process_trajectory(d1, ct_v_M, "M_me1"),
    process_trajectory(d2, ct_v_B, "B_me2"),
    process_trajectory(d2, ct_v_E, "E_me2"),
    process_trajectory(d2, ct_v_M, "M_me2"),
    process_trajectory(d3, ct_v_B, "B_me3"),
    process_trajectory(d3, ct_v_E, "E_me3"),
    process_trajectory(d3, ct_v_M, "M_me3")
) %>% rbindlist()

## count of shared genes between trajectories
res[, marker := substr(group, 3, 5)]
res[, lineage := substr(group, 1, 1)]


d_wide = dcast(res, gene ~ group, value.var = "trans", fill = "None")

write_tsv(d_wide, "./tmp/table_gene_chromatin_status_trajectory_remove_dup_wide.tsv")
write_tsv(res, "./tmp/table_gene_chromatin_status_trajectory_remove_dup_long.tsv")


# Count events by each marker and lineage
d_summary = data.table(table(res[, c("group", "trans")]))[order(group, -N)]

d_summary_wide = dcast(d_summary, trans ~ group, value.var = "N", fill = 0)
d_summary_wide[order(-B_me2)][1:30]

write_tsv(d_summary, "./tmp/table_gene_chromatin_status_trajectory_remove_dup_summary.tsv")
write_tsv(d_summary_wide, "./tmp/table_gene_chromatin_status_trajectory_remove_dup_summary_wide.tsv")

pdf("./figures/figure_chromatin_status_trajectory_barplot.pdf", width = 10, height = 8)
d_plot = d_summary[trans %in% c("Active->Un->Repressive", "Active->Bivalent->Repressive", "Active->Repressive")][order(-N)]

d_plot[, lineage := substr(group, 1, 1)]
d_plot[, marker := substr(group, 3, 5)]

ggplot(d_plot, aes(x = lineage, y = N, fill = trans)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    scale_fill_nejm(alpha = 0.6) +
    labs(
	title = "The Path to Repressive from Active",
	x = "Lineage",
	y = "Number of genes"
	) +
    theme(
	text = element_text(size = 16),
	legend.position = "top"
	) +
    facet_wrap(~marker, ncol = 3)

## Genes status at least assigned with two markers
d_plot = res[, .N, by = .(gene, trans, lineage)][N > 1]
d_plot = d_plot[trans %in% c("Active->Un->Repressive", "Active->Bivalent->Repressive", "Active->Repressive")][order(-N)]
d_plot = d_plot[, .N, by = .(trans, lineage)]

ggplot(d_plot, aes(x = lineage, y = N, fill = trans)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    labs(
	title = "Genes with consistent chromatin status between markers",
	x = "Lineage",
	y = "Number of genes"
	) +
    theme(
	text = element_text(size = 16),
	legend.position = "top"
    )
dev.off()


d_summary[, rank := rank(-N), by = group]

ggplot(d_summary, aes(x = rank, y = N, color = group)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    geom_hline(yintercept = 250, linetype = "dashed", color = "grey") +
    theme_classic() +
    labs(title = "Gene chromatin status trajectory",
	x = "Ranked trajectory",
	y = "Number of genes (log10)") +
    theme(
	text = element_text(size = 16),
	legend.position = "top"
    )

v_trans_selected = d_summary[N > 200, trans] %>% unique

## heatmap of total genes in each trajecotry
d_plot = res[trans %in% v_trans_selected][, .N, by = .(trans, group)] %>%
    dcast(trans ~ group, value.var = "N", fill = 0)

library(pheatmap)
mat = as.matrix(d_plot[, -1])
# mat = mat[order(rowSums(mat)), ]
rownames(mat) = d_plot$trans
mat = log10(mat + 1)
# mat = mat[, c("B_me1", "M_me1", "E_me1",
    # "B_me2", "M_me2", "E_me2",
    # "B_me3", "M_me3", "E_me3")]

pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    # display_numbers = TRUE,
    fontsize_number = 12,
    fontsize = 14,
    scale = "none",
    gaps_col = c(3,6),
    main = "Gene chromatin status trajectory\n(trajectories with > 200 genes)"
)


marker_v = c("me1", "me2", "me3")

d_g = expand.grid(marker_v, v_trans_selected)

marker_v = c("me1", "me2", "me3")
lineage_v = c("B", "E", "M")

d_g = expand.grid(lineage_v, v_trans_selected)

lapply(1:nrow(d_g), function(k) {
    trans_i = d_g[k, 2]
    lineage_j = d_g[k, 1]

    d_sub = res[lineage == lineage_j]

    d_sub_sub = d_sub[trans == trans_i]

    y = d_sub_sub[, .(marker = paste(sort(marker), collapse = "+")), by = gene]
    y[, lineage := lineage_j]
    y[, trans := trans_i]

    y
}) %>% rbindlist -> d_shared_genes

d_shared_genes

write_tsv(d_shared_genes, "./tmp/table_shared_genes_between_markers_in_trajectory.tsv")


marker_v = c("me1", "me2", "me3")
lineage_v = c("B", "E", "M")

d_g = expand.grid(lineage_v, unique(res$trans))

lapply(1:nrow(d_g), function(k) {
	trans_i = d_g[k, 2]
	lineage_j = d_g[k, 1]

	d_sub = res[lineage == lineage_j]

	d_sub_sub = d_sub[trans == trans_i]

	y = d_sub_sub[, .(marker = paste(sort(marker), collapse = "+")), by = gene]
	y[, lineage := lineage_j]
	y[, trans := trans_i]

	y
}) %>% rbindlist -> d_shared_genes

d_shared_genes

write_tsv(d_shared_genes, "./tmp/table_shared_genes_between_markers_in_trajectory_full_trans.tsv")


## Draw heatmap
library(pheatmap)

Cairo::CairoPDF("./figures/heatmap_shared_genes_between_lineages_in_trajectory.pdf", width = 12, height = 8)

marker_v = c("me1", "me2", "me3")

d_g = expand.grid(marker_v, v_trans_selected)

lapply(1:nrow(d_g), function(k) {
	trans_i = d_g[k, 2]
	marker_j = d_g[k, 1]
	d_sub = res[marker == marker_j]

	d_sub_sub = d_sub[trans == trans_i]

	y = d_sub_sub[, .(status = paste(sort(lineage), collapse = "+")), by = gene]
	y = y[ , .N, by = status][order(-N)]
	y[, marker := marker_j]
	y[, trans := trans_i]

	y
}) %>% rbindlist -> d_shared

row_order = c("B+E+M", "B+E", "B+M", "E+M", "B", "E", "M")

v_n_status = v_trans_selected %>% str_split("->") %>% lapply(length) %>% unlist
trans_order = v_trans_selected[order(rank(v_trans_selected) + 100 * rank(v_n_status))]


plot_shared_genes_heatmap <- function(d_shared, marker_i, trans_order, main_title) {
    d_plot <- dcast(d_shared[marker == marker_i], trans ~ status, value.var = "N", fill = 0)
    mat <- as.matrix(d_plot[, -1])
    rownames(mat) <- d_plot$trans
    # mat <- log10(mat + 1)
    mat <- mat[, c("B+E+M", "B+E", "B+M", "E+M", "B", "E", "M")]
    mat <- mat[trans_order, ]
    pheatmap(
	mat,
	cluster_rows = FALSE,
	cluster_cols = FALSE,
	display_numbers = TRUE,
	number_format = "%.1f",
	fontsize_number = 12,
	fontsize = 14,
	scale = "none",
	gaps_col = c(1, 4),
	# gaps_row = c(4, 13, 24),
	main = main_title
	# col = colorRampPalette(c("white", "Navy"))(100)
    )
}

plot_shared_genes_heatmap(
    d_shared,
    marker = "me1",
    trans_order = trans_order,
    main_title = "Shared genes between lineages in \n chromatin status trajectory (ME1)"
)

plot_shared_genes_heatmap(
    d_shared,
    marker = "me2",
    trans_order = trans_order,
    main_title = "Shared genes between lineages in \n chromatin status trajectory (ME2)"
)

plot_shared_genes_heatmap(
    d_shared,
    marker = "me3",
    trans_order = trans_order,
    main_title = "Shared genes between lineages in \n chromatin status trajectory (ME3)"
)


## Plot the shared genes_heatmap overlap between markers
plot_shared_genes_heatmap_share <- function(d, trans_order, main_title) {
    d_plot <- dcast(d, trans ~ status, value.var = "N", fill = 0)
    mat <- as.matrix(d_plot[, -1])
    rownames(mat) <- d_plot$trans
    # mat <- log10(mat + 1)
    mat <- mat[, c("B+E+M", "B+E", "B+M", "E+M", "B", "E", "M")]
    trans_order = trans_order[trans_order %in% rownames(mat)]
    mat <- mat[trans_order, ]
    pheatmap::pheatmap(
	mat,
	cluster_rows = T,
	cluster_cols = FALSE,
	display_numbers = FALSE,
	number_format = "%.1f",
	fontsize_number = 12,
	fontsize = 14,
	scale = "row",
	gaps_col = c(1, 4),
	# gaps_row = c(4, 13, 24),
	main = main_title,
	col = colorRampPalette(c("white", "red"))(100)
    )
}

## Genes status at least assigned with two markers
d_plot = res[, .N, by = .(gene, trans, lineage)][N > 1][, .(status = paste(sort(lineage), collapse = "+")), by = .(gene, trans)]

d_plot[gene %in% "PROM1"]

d_plot = d_plot[, .N, .(trans, status)][order(-N)]


plot_shared_genes_heatmap_share(
    d_plot,
    trans_order = trans_order,
    main_title = "Shared genes between lineages in \n chromatin status trajectory (ME1, ME2, ME3 least two agree)"
)

dev.off()

# all_genes <- unique(d_plot[grepl("^Biv", trans), gene]) {{{
#
# run_fisher_summary <- function(trans_val, status_val) {
#   x_genes <- d_plot[trans == trans_val & status == status_val, gene]
#   res <- d_plot[gene %in% x_genes][, .N, .(trans, status)][order(-N)]
#   res <- res[grepl("^Biv", trans)]
#   if (nrow(res) == 0) return(res)
#   x1 <- unique(d_plot$gene[d_plot$trans == "Bivalent->Active" & d_plot$status == status_val])
#   ## pvalue
#   res[, p_value := sapply(1:.N, function(i) {
#     x2 <- unique(d_plot$gene[d_plot$trans == res$trans[i] & d_plot$status == res$status[i]])
#     fisher_overlap_test(x1, x2, all_genes)$fisher_result$p.value
#   })]
#   ## Odds ratio
#   res[, odds_ratio := sapply(1:.N, function(i) {
# 	x2 <- unique(d_plot$gene[d_plot$trans == res$trans[i] & d_plot$status == res$status[i]])
# 	fisher_overlap_test(x1, x2, all_genes)$fisher_result$estimate
#   })]
#   res$p.adjust <- p.adjust(res$p_value, method = "BH")
#   res
# }
#
# res_B <- run_fisher_summary("Bivalent->Active", "B")
# res_E <- run_fisher_summary("Bivalent->Active", "E")
# res_M <- run_fisher_summary("Bivalent->Active", "M")
# res_B
# res_E
# res_M
#
# res_B <- run_fisher_summary("Bivalent->Repressive", "B")
# res_E <- run_fisher_summary("Bivalent->Repressive", "E")
# res_M <- run_fisher_summary("Bivalent->Repressive", "M")
# res_B
# res_E
# res_M }}}


# get_status <- function(x, status_v) { {{{
#     if (status_v[1] %in% x) {
# 	return(status_v[1])
#     } else if (status_v[2] %in% x) {
# 	return(status_v[2])
#     } else if (status_v[3] %in% x) {
# 	return(status_v[3])
#     } else {
# 	return(status_v[4])
#     }
# }
#
# # Read data
# step_m <- fread("./tmp/table_chromatin_status_transitions_M_me1.tsv")
# step_e <- fread("./tmp/table_chromatin_status_transitions_E_me1.tsv")
# step_b <- fread("./tmp/table_chromatin_status_transitions_B_me1.tsv")
#
# # Remove 'Other' genes
# gene_other <- c(
#     step_m[next_cat == "Other" | cat == "Other", gene],
#     step_e[next_cat == "Other", gene],
#     step_b[next_cat == "Other", gene]
#     ) %>% unique
#
# step_m <- step_m[!gene %in% gene_other]
# step_e <- step_e[!gene %in% gene_other]
# step_b <- step_b[!gene %in% gene_other]
#
# step_b
#
# # Define gene sets
# gene_hspc_biv_v <- step_b[cell_type == "HSPC" & cat == "Bivalent", gene]
# gene_hspc_un_v  <- step_b[cell_type == "HSPC" & cat == "Un", gene]
# gene_hspc_act_v <- step_b[cell_type == "HSPC" & cat == "Active", gene]
# gene_hspc_rep_v <- step_b[cell_type == "HSPC" & cat == "Repressive", gene]
#
# # Status ranking
# rank_bi_repressive <- c("Repressive", "Active", "Un", "Bivalent")
# rank_bi_active     <- c("Active", "Repressive", "Un", "Bivalent")
#
# files = dir("./tmp/", pattern = "table_upset_data_.*_ME2.tsv", full.names = T)
# name = str_match(files, "table_upset_data_(.*)_ME2.tsv")[,2]
# d = lapply(1:length(files), function(i) {
#     d = fread(files[i])
#     d[, module := name[i]]
#     return(d)
# }) %>% rbindlist
#
# r2a = d[module %in% c("Repressive -> Active") & B == T]
# r2b = d[module %in% c("Repressive -> Bivalent") & B == T]
# r2u = d[module %in% c("Repressive -> Unmarked") & B == T]
#
# l_plot = list(
# 	r2a = r2a$gene,
# 	r2b = r2b$gene,
# 	r2u = r2u$gene
# )
#
# library(eulerr)
#
# fit <- euler(l_plot)
# p <- plot(
#     fit,
#     fills = c("red", "green", "blue"),
#     alpha = 0.5,
#     labels = list(cex = 2),
#     quantities = list(cex = 2),
# )
# p
# }}}


