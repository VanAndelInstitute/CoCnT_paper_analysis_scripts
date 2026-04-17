# =============================================================================
# Script: run_221_gaussian_mixture_model_cutoff.R
#
# Input:
#   - table_cell_cluster_annotation_final.tsv: cell-type/cluster annotations.
#   - table_gene_marker_average_score_by_celltype.tsv.gz: marker gene scores by cell type (compressed).
#
# Output:
#   - A PDF file ('./figures/figure_supplementary_gaussian_mixture_model_cutoff.pdf') visualizing the GMM fit and cutoffs for each marker.
# =============================================================================

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d_ann = fread("./tmp/table_metadata_h3k27me3_final3.tsv")
d_exp = fread("./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")

library(mclust)
library(ggplot2)

marker_v = unique(d_exp$marker)

cell_type_included = c(
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
    "Erythroid progenitors3"
)


pdf("./figures/figure_supplementary_gaussian_mixture_model_cutoff.pdf", width = 7, height = 5)
for (i in 1:length(marker_v)) {
    print(i)
    set.seed(1234)
    marker_x = marker_v[i]
    d_exp_sub = d_exp[marker == marker_x & cell_type %in% cell_type_included]
    gene_score = d_exp_sub$gene_score

    x <- gene_score
    # x = x[x > 0]
    x_log <- log10(x + 0.1)

    fit <- Mclust(x_log, G = 2)

    params <- data.frame(
	component = factor(1:2),
	mean = fit$parameters$mean,
	sd   = sqrt(fit$parameters$variance$sigmasq),
	weight = fit$parameters$pro
    )


    xx <- seq(min(x_log), max(x_log), length.out = 1000)

    dens_df <- do.call(rbind, lapply(1:2, function(k) {
	    data.frame(
		x = xx,
		density = params$weight[k] * dnorm(xx, params$mean[k], params$sd[k]),
		component = paste0("Component ", k)
	    )
    }))

    mixture_df <- aggregate(density ~ x, dens_df, sum)

    f1 <- function(z) params$weight[1] * dnorm(z, params$mean[1], params$sd[1])
    f2 <- function(z) params$weight[2] * dnorm(z, params$mean[2], params$sd[2])

    cutoff <- uniroot(function(z) f1(z) - f2(z),
	interval = range(x_log))$root

    gene_n = nrow(d_exp_sub)
    gene_above_cutoff = sum(x_log > cutoff)

    g = ggplot(data.frame(x = x_log), aes(x)) +

	# Histogram (density-scaled!)
	geom_histogram(aes(y = after_stat(density)),
	    bins = 80,
	    fill = "grey90",
	    color = "white") +

	# Individual Gaussian components
	geom_line(data = dens_df,
	    aes(x, density, color = component),
	    linewidth = 1.2) +

	# Total mixture
	geom_line(data = mixture_df,
	    aes(x, density),
	    color = "black",
	    linewidth = 1.2) +

	scale_color_manual(values = c("#1f77b4", "#d62728")) +

	labs(
	    x = "log10(gene_score + 0.1)",
	    y = "Density",
	    title = "Gaussian mixture model"
	    ) +

	theme_classic() +
	theme(legend.title = element_blank()) + 
	geom_vline(xintercept = cutoff,
	    linetype = "dashed",
	    linewidth = 0.8) +
	labs(title= paste0("Marker: ", marker_x,
		"; Cutoff: ", round(10^cutoff - 0.01, 3)),
	    subtitle = paste0("Total genes * celltype: ", gene_n,
		"; Genes *celltype above cutoff: ", gene_above_cutoff))

	print(g)
}
dev.off()
