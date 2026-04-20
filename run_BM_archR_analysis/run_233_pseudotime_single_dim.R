## script: run_233_pseudotime_single_dim.R
## purpose: plot gene score along pseudotime for genes in different chromatin status transition categories, for each marker separately
## input:
## - Gene score matrix (imputed) for each marker
## * ./tmp/list_matrix_imputed_gene_score.rds
## - Chromatin status trajectory consensus table
## * ./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv
## - Metadata with pseudotime for each lineage
## * ./tmp/table_metadata_h3k4me27_final3_HSPC.tsv
## output:
## - Pseudotime plot for each marker and lineage, comparing Bivalent->Active vs Un->Active genes
## * ./figures/233_pseudotime_single_dim.pdf
## * ./figures/233_pseudotime_single_dim_scale.pdf
## 
library(data.table)
library(ggplot2)
library(ggpubr)
# library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggrastr)

d_meta <- fread("./tmp/table_metadata_h3k4me27_final3_HSPC.tsv")
m_list <- readRDS("./tmp/list_matrix_imputed_gene_score.rds")

gene_category = fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")

marker_colors <- c(
    H3K27me3 = "#3B8FC4",
    H3K4me1 = "#FFD000",
    H3K4me2 = "#F39C12",
    H3K4me3 = "#E74C3C",
    H3K4me1_cooc = "#56B870",
    H3K4me2_cooc = "#D16BA5",
    H3K4me3_cooc = "#9B59B6"
)

trajectory_map <- list(
    "B" = list(name = "Bcell_Trajectory"),
    "M" = list(name = "Monocyte_Trajectory"),
    "E" = list(name = "Erythroid_Trajectory")
)

available_markers <- intersect(names(m_list), names(marker_colors))
available_genes <- rownames(m_list[[1]])
if (is.null(available_genes)) stop("m_list[[1]] must have rownames as genes.")

threshold_v = c(
    H3K27me3 = 0.21,
    H3K4me1 = 0.305,
    H3K4me2 = 0.139,
    H3K4me3 = 0.181,
    H3K4me1_cooc = 0447,
    H3K4me2_cooc = 0.314,
    H3K4me3_cooc = 0.455
)

# -------------------------
# Helpers
# -------------------------
normalize_gene_matrix <- function(mat, method = c("z", "minmax")) {
    method <- match.arg(method)

    if (method == "z") {
	# z-score per gene (row)
	mu <- rowMeans(mat, na.rm = TRUE)
	sdv <- apply(mat, 1, sd, na.rm = TRUE)
	sdv[sdv == 0 | is.na(sdv)] <- 1
	mat <- sweep(mat, 1, mu, "-")
	mat <- sweep(mat, 1, sdv, "/")
    }

    if (method == "minmax") {
	rmin <- apply(mat, 1, min, na.rm = TRUE)
	rmax <- apply(mat, 1, max, na.rm = TRUE)
	denom <- rmax - rmin
	denom[denom == 0 | is.na(denom)] <- 1
	mat <- sweep(mat, 1, rmin, "-")
	mat <- sweep(mat, 1, denom, "/")
    }

    mat
}

parse_genes_from_text <- function(txt) {
    # one gene per line, but also tolerate commas/spaces
    g <- unlist(strsplit(txt, "[,\n\r\t ]+"))
    g <- toupper(trimws(g))
    g <- g[nzchar(g)]
    unique(g)
}

make_geneset_score_dt <- function(
    genes, 
    markers, 
    allow_missing = TRUE,
    normalize = FALSE,
    norm_method = c("z", "minmax")
    ) {
    norm_method <- match.arg(norm_method)
    genes <- toupper(trimws(genes))
    genes <- genes[nzchar(genes)]
    genes <- unique(genes)
    if (length(genes) == 0) stop("No genes provided.")

    markers <- intersect(markers, names(m_list))
    if (length(markers) == 0) stop("No valid markers selected.")

    present <- intersect(genes, available_genes)
    missing <- setdiff(genes, available_genes)

    if (!allow_missing && length(missing) > 0) {
	stop(sprintf(
		"Missing genes (not in matrix rownames): %s",
		paste(missing, collapse = ", ")
		))
    }
    if (length(present) == 0) {
	stop("None of the input genes are found in available_genes.")
    }

    # For each marker: avg across genes per cell (colMeans over selected rows)
    dt_list <- lapply(markers, function(mk) {
	m <- m_list[[mk]]
	idx <- intersect(present, rownames(m))
	if (length(idx) == 0) return(NULL)


	sub <- as.matrix(m[idx, , drop = FALSE])

	if (normalize) {
	    sub <- normalize_gene_matrix(sub, method = norm_method)
	}

	# v <- colMeans(sub, na.rm = TRUE)
	v = apply(sub, 2, median, na.rm = TRUE)

	data.table(
	    cell_id = names(v),
	    gene_score = as.numeric(v),
	    marker = mk
	)
    })

    d_out <- rbindlist(dt_list, fill = TRUE)
    if (nrow(d_out) == 0) stop("No data produced for the selected markers/genes.")

    list(
	dt = d_out,
	present = present,
	missing = missing
    )
}

make_gene_score_dt <- function(gene, markers) {
    gene <- toupper(trimws(gene))
    if (!nzchar(gene)) stop("Gene is empty.")

    if (!(gene %in% available_genes)) {
	stop(sprintf("Gene '%s' not found in m_list rownames.", gene))
    }

    markers <- intersect(markers, names(m_list))
    if (length(markers) == 0) stop("No valid markers selected.")

    lapply(markers, function(mk) {
	m <- m_list[[mk]]
	if (!(gene %in% rownames(m))) return(NULL)
	v <- m[gene, ]
	data.table(
	    cell_id = names(v),
	    gene_score = as.numeric(v),
	    marker = mk
	)
    }) %>% rbindlist(fill = TRUE)
}

make_pseudotime_dt <- function(lineage_label) {
    trj_col <- trajectory_map[[lineage_label]]$name
    if (!(trj_col %in% colnames(d_meta))) {
	stop(sprintf("Trajectory column '%s' not found in ArchR cellColData.", trj_col))
    }
    d_meta[, .(cell_id, pseudotime = get(trj_col))] #%>% na.omit()
}

plot_gene_vs_pseudotime <- function(d_plot, gene, lineage_label, selected_markers) {
    d_plot[, marker := factor(marker, levels = selected_markers)]

    ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = marker)) +
	geom_point(alpha = 0.5, size = 1) +
	geom_smooth(method = "gam", se = FALSE) +
	scale_color_manual(values = marker_colors, drop = FALSE) +
	theme_classic() +
	theme(
	    legend.position = "right",
	    axis.title = element_blank()
	    ) +
	ggtitle(sprintf("%s — %s pseudotime", gene, lineage_label))
}

#  analysis for Bivalent->Active, Un->Active {{{
B2A = gene_category[trans == "Bivalent->Active"]
U2A = gene_category[trans == "Un->Active"]
BAB = gene_category[trans == "Bivalent->Active->Bivalent"]
UAU = gene_category[trans == "Un->Active->Un"]

default_markers <- c("H3K27me3")

marker_v = c("H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me1_cooc", "H3K4me2_cooc", "H3K4me3_cooc")

pdf("./figures/233_pseudotime_single_dim.pdf", width = 5, height = 3)

for (mk in marker_v) {
    default_markers = mk

    for (lineage_x in names(trajectory_map)) {

	group_x = paste0(lineage_x, "_consensus")

	gene_set_l = list(
	    "Bivalent->Active" = c(B2A[group == group_x, gene], BAB[group == group_x, gene]),
	    "Un->Active" = c(U2A[group == group_x, gene], UAU[group == group_x, gene])
	    # "Bivalent->Active->Bivalent" = BAB[group == group_x, gene],
	    # "Un->Active->Un" = UAU[group == group_x, gene]
	)

	lapply(1:length(gene_set_l), function(i) {
	    res = make_geneset_score_dt(
		genes = gene_set_l[[i]],
		markers = default_markers,
		allow_missing = TRUE,
		normalize = F
	    )

	    res$dt[, category := names(gene_set_l)[i]]
	}) %>% rbindlist() -> d_plot

	d_plot = merge(d_plot, make_pseudotime_dt(lineage_x), by = "cell_id", all.x = TRUE)

	p = ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = category)) +
	    rasterise(geom_point(alpha = 0.3, size = 0.1), dpi = 600) +
	    geom_smooth(method = "loess", se = FALSE, span = 0.3) +
	    scale_color_manual(values = c("Bivalent->Active" = "#e18727", "Un->Active" = "#1872b5")) +
	    theme_classic() +
		theme(
		    legend.position = "right",
		    axis.title = element_blank()
		    ) +
		ggtitle(sprintf("%s — %s pseudotime", mk, lineage_x))
	print(p)
    }
}

dev.off()

pdf("./figures/233_pseudotime_single_dim_scale.pdf", width = 5, height = 3)

for (mk in marker_v) {
    default_markers = mk

    for (lineage_x in names(trajectory_map)) {

	group_x = paste0(lineage_x, "_consensus")

	gene_set_l = list(
	    "Bivalent->Active" = c(B2A[group == group_x, gene], BAB[group == group_x, gene]),
	    "Un->Active" = c(U2A[group == group_x, gene], UAU[group == group_x, gene])
	    # "Bivalent->Active->Bivalent" = BAB[group == group_x, gene],
	    # "Un->Active->Un" = UAU[group == group_x, gene]
	)

	lapply(1:length(gene_set_l), function(i) {
	    res = make_geneset_score_dt(
		genes = gene_set_l[[i]],
		markers = default_markers,
		allow_missing = TRUE,
		normalize = T,
		norm_method = "z"
	    )

	    res$dt[, category := names(gene_set_l)[i]]
	}) %>% rbindlist() -> d_plot

	d_plot = merge(d_plot, make_pseudotime_dt(lineage_x), by = "cell_id", all.x = TRUE)

	p = ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = category)) +
	    rasterise(geom_point(alpha = 0.3, size = 0.1), dpi = 600) +
	    geom_smooth(method = "loess", se = FALSE, span = 0.3) +
	    scale_color_manual(values = c("Bivalent->Active" = "#e18727", "Un->Active" = "#1872b5")) +
	    theme_classic() +
		theme(
		    legend.position = "right",
		    axis.title = element_blank()
		    ) +
		ggtitle(sprintf("%s — %s pseudotime", mk, lineage_x))
	print(p)
    }
}

dev.off()

## }}}

############
#  Backup  #
############


# # Bi-markers analysis for Active->Un->Repressive, Active->Bivalent->Repressive, Active->Repressive {{{
#
# pdf("./figures/289_pseudotime_single_bi_markers_dim_consensus_2.pdf", width = 9, height = 2)
#
# gene_consensus = fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")
# gene_consensus[, lineage := substr(group, 1, 1)]
#
#
# out_l = list()
#
# mk_s = "me2"
# mk_full = c(
#     paste0("H3K4", mk_s),
#     paste0("H3K4", mk_s, "_cooc"),
#     "H3K27me3"
# )
#
# ABR = gene_consensus[trans == "Active->Bivalent->Repressive"]
# AUR = gene_consensus[trans == "Active->Un->Repressive"]
# AR = gene_consensus[trans == "Active->Repressive"]
#
# threshold_v_sub = threshold_v[c(paste0("H3K4", mk_s), paste0("H3K4", mk_s, "_cooc"), "H3K27me3")]
#
#
# for (lineage_x in names(trajectory_map)) {
#     gene_set_l = list(
# 	"Active->Bivalent->Repressive" = ABR[lineage == lineage_x, gene],
# 	"Active->Un->Repressive" = AUR[lineage == lineage_x, gene],
# 	"Active->Repressive" = AR[lineage == lineage_x, gene]
#     )
#
#     lapply(1:length(gene_set_l), function(i) {
# 	res = make_geneset_score_dt(
# 	    genes = gene_set_l[[i]],
# 	    markers = mk_full,
# 	    allow_missing = TRUE,
# 	    normalize = FALSE
# 	)
#
# 	res$dt[, category := names(gene_set_l)[i]]
# 	res$dt
#     }) %>% rbindlist() -> d_plot
#
#     d_plot = merge(d_plot, make_pseudotime_dt(lineage_x), by = "cell_id", all.x = TRUE)
#
#     d_plot$mk_s = mk_s
#     d_plot$lineage = lineage_x
#
#     print(mk_s)
#     print(lineage_x)
#
#     out_l = c(out_l, list(d_plot))
#
#     p = ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = category)) +
# 	geom_point(alpha = 0.3, size = 0.1) +
# 	geom_smooth(method = "gam", se = FALSE) +
# 	geom_hline(
# 	    data = data.frame(marker = names(threshold_v_sub),
# 		yintercept = threshold_v_sub),
# 	    aes(yintercept = yintercept),
# 	    linetype = "dashed",
# 	    color = "black"
# 	    ) +
# 	facet_wrap(~ marker, nrow = 1, scales = "free_y") +
# 	scale_color_manual(values = c(
# 		"Active->Bivalent->Repressive" = "#ecb682", 
# 		"Active->Un->Repressive" = "#a4cce5",
# 		"Active->Repressive" = "#d88c80"
# 		)) +
# 	theme_classic() +
# 	theme(
# 	    legend.position = "right",
# 	    axis.title = element_blank()
# 	    ) +
# 	ggtitle(sprintf("%s — %s pseudotime", mk_full, lineage_x))
#     print(p)
# }
#
# out_l %>% rbindlist() -> d_all
#
# write_tsv(d_all, "./tmp/table_a2p_bi_marker_subset_consensus.tsv")
#
# dev.off()
# ## }}}

# #  analysis for Active->Un->Repressive, Active->Bivalent->Repressive, Active->Repressive {{{
#
# pdf("./figures/289_pseudotime_single_dim_2.pdf", width = 9, height = 2)
#
# marker_v_short = c("me1", "me2", "me3")
#
# out_l = list()
#
# for (mk_s in marker_v_short) {
#     mk_full = c(
# 	paste0("H3K4", mk_s),
# 	paste0("H3K4", mk_s, "_cooc"),
# 	"H3K27me3"
# 	)
#
#     ABR = gene_category[trans == "Active->Bivalent->Repressive" & grepl(mk_s, marker)]
#     ABR[marker == "me2"]
#     AUR = gene_category[trans == "Active->Un->Repressive" & grepl(mk_s, marker)]
#     AR = gene_category[trans == "Active->Repressive" & grepl(mk_s, marker)]
#
#
#     for (lineage_x in names(trajectory_map)) {
# 	gene_set_l = list(
# 	    "Active->Bivalent->Repressive" = ABR[lineage == lineage_x, gene],
# 	    "Active->Un->Repressive" = AUR[lineage == lineage_x, gene],
# 	    "Active->Repressive" = AR[lineage == lineage_x, gene]
# 	)
#
# 	lapply(1:length(gene_set_l), function(i) {
# 	    res = make_geneset_score_dt(
# 		genes = gene_set_l[[i]],
# 		markers = mk_full,
# 		allow_missing = TRUE,
# 		normalize = FALSE
# 	    )
#
# 	    res$dt[, category := names(gene_set_l)[i]]
# 	    res$dt
# 	}) %>% rbindlist() -> d_plot
#
# 	d_plot = merge(d_plot, make_pseudotime_dt(lineage_x), by = "cell_id", all.x = TRUE)
#
# 	d_plot$mk_s = mk_s
# 	d_plot$lineage = lineage_x
#
# 	print(mk_s)
# 	print(lineage_x)
#
# 	out_l = c(out_l, list(d_plot))
#
# 	p = ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = category)) +
# 	    geom_point(alpha = 0.3, size = 0.1) +
# 	    geom_smooth(method = "gam", se = FALSE) +
# 	    geom_hline(
# 		data = data.frame(marker = names(threshold_v),
# 		    yintercept = threshold_v),
# 		aes(yintercept = yintercept),
# 		linetype = "dashed",
# 		color = "black"
# 		) +
# 	    facet_wrap(~ marker, nrow = 1, scales = "free_y") +
# 	    scale_color_manual(values = c(
# 		    "Active->Bivalent->Repressive" = "#ecb682", 
# 		    "Active->Un->Repressive" = "#a4cce5",
# 		    "Active->Repressive" = "#d88c80"
# 		    )) +
# 	    theme_classic() +
# 	    theme(
# 		legend.position = "right",
# 		axis.title = element_blank()
# 		) +
# 	    ggtitle(sprintf("%s — %s pseudotime", mk_full, lineage_x))
# 	print(p)
#     }
# }
#
# out_l %>% rbindlist() -> d_all
#
#
# write_tsv(d_all, "./tmp/table_a2p_subset.tsv")
#
# dev.off()
# ## }}}

# #  Two axis figure  {{{
#
# d_all = fread("./tmp/table_a2p_subset.tsv")
#
# ## H3k4me2
#
# scale_factor = 0.139 / 0.21 
#
# Cairo::CairoPDF("./figures/289_pseudotime_single_dim_2_axis.pdf", width = 12, height = 3.5)
#
# d_plot = d_all[lineage == "B" & marker %in% c("H3K4me2", "H3K27me3") & mk_s == "me2"]
#
# ggplot(d_plot[marker == "H3K4me2"], aes(x = pseudotime, y = gene_score, color = marker)) +
# 	geom_point(alpha = 0.3, size = 0.1) +
# 	geom_point(data = d_plot[marker == "H3K27me3"], aes(x = pseudotime, y = gene_score * scale_factor, color = marker), alpha = 0.3, size = 0.1) +
# 	geom_smooth(method = "gam", se = FALSE) +
# 	geom_smooth(data = d_plot[marker == "H3K27me3"], method = "gam", se = FALSE, aes(x = pseudotime, y = gene_score * scale_factor, color = marker)) +
# 	scale_color_manual(values = c("H3K4me2" = "#a97c50", "H3K27me3" = "#006838")) +
# 	facet_wrap(~ category, nrow = 1, scales = "free_y") +
# 	expand_limits(x = 0, y = 0) +
# 	geom_hline(yintercept = 0.139, linetype = "dashed") +
# 	theme_classic() +
# 	theme(
# 	    legend.position = "right",
# 	    axis.title.y.right = element_text(color = "#006838"),
# 	    axis.title.y.left = element_text(color = "#a97c50")
# 	) +
# 	labs(x = "Pseudotime", y = "H3K4me2 gene score", title = "B lineage (me2)") +
# 	scale_y_continuous(
# 	    name = "H3K27me3 gene score",
# 	    sec.axis = sec_axis(~ ./scale_factor, name = "H3K27me3 gene score")
# 	)
#
# d_plot = d_all[lineage == "E" & marker %in% c("H3K4me2", "H3K27me3") & mk_s == "me2"]
#
# ggplot(d_plot[marker == "H3K4me2"], aes(x = pseudotime, y = gene_score, color = marker)) +
# 	geom_point(alpha = 0.3, size = 0.1) +
# 	geom_point(data = d_plot[marker == "H3K27me3"], aes(x = pseudotime, y = gene_score * scale_factor, color = marker), alpha = 0.3, size = 0.1) +
# 	geom_smooth(method = "gam", se = FALSE) +
# 	geom_smooth(data = d_plot[marker == "H3K27me3"], method = "gam", se = FALSE, aes(x = pseudotime, y = gene_score * scale_factor, color = marker)) +
# 	scale_color_manual(values = c("H3K4me2" = "#a97c50", "H3K27me3" = "#006838")) +
# 	facet_wrap(~ category, nrow = 1, scales = "free_y") +
# 	expand_limits(x = 0, y = 0) +
# 	geom_hline(yintercept = 0.139, linetype = "dashed") +
# 	theme_classic() +
# 	theme(
# 	    legend.position = "right",
# 	    axis.title.y.right = element_text(color = "#006838"),
# 	    axis.title.y.left = element_text(color = "#a97c50")
# 	) +
# 	labs(x = "Pseudotime", y = "H3K4me2 gene score", title = "E lineage (me2)") +
# 	scale_y_continuous(
# 	    name = "H3K27me3 gene score",
# 	    sec.axis = sec_axis(~ ./scale_factor, name = "H3K27me3 gene score")
# 	)
#
# dev.off()
#
#
# #######################
# #  CPG island or not  #
# #######################
# d_cg_gene = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")
#
# marker_v_short = c("me1", "me2", "me3")
#
# pdf("./figures/289_pseudotime_single_dim_2_cpg.pdf", width = 12, height = 3.5)
#
# out_l = list()
#
# for (mk_s in marker_v_short) {
#     mk_full = c(
# 	paste0("H3K4", mk_s),
# 	paste0("H3K4", mk_s, "_cooc"),
# 	"H3K27me3"
# 	)
#
#     ABR = gene_category[trans == "Active->Bivalent->Repressive" & grepl(mk_s, marker)]
#     AUR = gene_category[trans == "Active->Un->Repressive" & grepl(mk_s, marker)]
#     AR = gene_category[trans == "Active->Repressive" & grepl(mk_s, marker)]
#
#
#     for (lineage_x in names(trajectory_map)) {
#
# 	all_gene_v = unique(c(
# 	    ABR[lineage == lineage_x, gene],
# 	    AUR[lineage == lineage_x, gene],
# 	    AR[lineage == lineage_x, gene]
# 	))
#
# 	gene_set_l = list(
# 	    "cpg_island" = all_gene_v[all_gene_v %in% d_cpg_gene],
# 	    "non_cpg_island" = all_gene_v[!(all_gene_v %in% d_cpg_gene)]
# 	)
#
# 	lapply(1:length(gene_set_l), function(i) {
# 	    res = make_geneset_score_dt(
# 		genes = gene_set_l[[i]],
# 		markers = mk_full,
# 		allow_missing = TRUE,
# 		normalize = FALSE
# 	    )
#
# 	    res$dt[, category := names(gene_set_l)[i]]
# 	    res$dt
# 	}) %>% rbindlist() -> d_plot
#
# 	d_plot = merge(d_plot, make_pseudotime_dt(lineage_x), by = "cell_id", all.x = TRUE)
#
# 	d_plot$mk_s = mk_s
# 	d_plot$lineage = lineage_x
#
# 	print(mk_s)
# 	print(lineage_x)
#
# 	out_l = c(out_l, list(d_plot))
#
# 	p = ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = category)) +
# 	    geom_point(alpha = 0.3, size = 0.1) +
# 	    geom_smooth(method = "gam", se = FALSE) +
# 	    facet_wrap(~ marker, nrow = 1, scales = "free_y") +
# 	    scale_color_manual(values = c(
# 		    "cpg_island" = "#a97c50",
# 		    "non_cpg_island" = "#006838"
# 		    )) +
# 	    theme_classic() +
# 	    theme(
# 		legend.position = "right",
# 		axis.title = element_blank()
# 		) +
# 	    ggtitle(sprintf("%s — %s pseudotime", mk_full, lineage_x))
# 	print(p)
#     }
# }
#
# out_l %>% rbindlist() -> d_all
#
#
# write_tsv(d_all, "./tmp/table_a2p_cpg_subset.tsv")
#
# dev.off()
#
# ## }}}

#  First dereive  # {{{
###################

# httpgd::hgd(port = 4323)
#
# d_all = fread("./tmp/table_a2p_subset.tsv")
#
# df_subset = d_all[lineage == "E" & mk_s == "me2"] %>% na.omit
#
# library(dplyr)
# library(purrr)
# library(mgcv)
# library(gratia)
# library(ggplot2)
#
# # 1) Nest into 9 subsets
# nested <- df_subset %>%
#     tidyr::nest(data = c(pseudotime, gene_score), .by = c(marker, category))
#
# # 2) Fit 9 GAMs
# nested <- nested %>%
#     mutate(
# 	m = map(data, ~ gam(gene_score ~ s(pseudotime, k = 10), data = .x, method = "REML"))
#     )
#
# # 3) Derivatives for each model
# nested <- nested %>%
#     mutate(
# 	d = map(m, ~ derivatives(.x))
#     )
#
# # 4) Combine derivatives into one dataframe
# D <- nested %>%
#     select(marker, category, d) %>%
#     tidyr::unnest(d)
#
# # 5) Plot derivatives, 3×3 facet grid
# ggplot(D, aes(x = pseudotime, y = .derivative, color = category)) +
#     geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.15) +
#     geom_line(linewidth = 1) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     facet_grid(marker ~. ) +
#     theme_classic() +
#     labs(x = "x", y = "dy/dx", title = "First derivatives (9 separate GAMs)")
#
#
#
# deriv <- derivatives(gam_fit)
#
# draw(deriv)
#
# ggplot(deriv, aes(x = data, y = derivative)) +
#     geom_line(size = 1) +
#     theme_classic() +
#     labs(x = "x", y = "First derivative (dy/dx)")
#
#
# d_all
# }}}

