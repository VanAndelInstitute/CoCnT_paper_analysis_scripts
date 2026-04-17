library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

httpgd::hgd(port = 4323)

d_list = fread("./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")

## Scatter plot for marker correlation
d_plot = dcast(d_list, gene + cell_type ~ marker, value.var = "gene_score")

# sample(1:nrow(d_plot), 10000, replace = F) -> sampled_idx
cell_type_v = unique(d_plot$cell_type)
sampled_idx = d_plot$cell_type == cell_type_v[4]

library(ggplot2)
library(patchwork)

# Helper function to create correlation hex plots 
make_corr_hex_plot <- function(df, xvar, yvar, bins = 70) {
    cor_val <- cor(df[[xvar]], df[[yvar]], method = "spearman")
    p_val <- cor.test(df[[xvar]], df[[yvar]], method = "spearman")$p.value
    label = paste0("r = ", round(cor_val, 2), "\np = ", signif(p_val, 2))
    ggplot(df, aes_string(x = paste0(xvar, " + 0.001"), y = paste0(yvar, " + 0.001"))) +
	geom_hex(bins = bins) +
	scale_fill_continuous(type = "viridis") +
	scale_x_log10() +
	scale_y_log10() +
	theme_classic() +
	annotate("text", x = Inf, y = Inf, label = label,
	    hjust = 1.1, vjust = 1.5, size = 5) 
}

# Define pairs to plot
corr_pairs <- list(
    c("H3K4me2", "H3K4me1"),
    c("H3K4me3", "H3K4me1"),
    c("H3K27me3", "H3K4me1"),
    c("H3K4me1_cooc", "H3K4me1"),
    c("H3K4me2_cooc", "H3K4me1"),
    c("H3K4me3_cooc", "H3K4me1"),
    c("H3K4me3", "H3K4me2"),
    c("H3K27me3", "H3K4me2"),
    c("H3K4me1_cooc", "H3K4me2"),
    c("H3K4me2_cooc", "H3K4me2"),
    c("H3K4me3_cooc", "H3K4me2"),
    c("H3K27me3", "H3K4me3"),
    c("H3K4me1_cooc", "H3K4me3"),
    c("H3K4me2_cooc", "H3K4me3"),
    c("H3K4me3_cooc", "H3K4me3"),
    c("H3K4me1_cooc", "H3K27me3"),
    c("H3K4me2_cooc", "H3K27me3"),
    c("H3K4me3_cooc", "H3K27me3"),
    c("H3K4me2_cooc", "H3K4me1_cooc"),
    c("H3K4me3_cooc", "H3K4me1_cooc"),
    c("H3K4me3_cooc", "H3K4me2_cooc")
)

## Generate plots for all genes {{{
plots <- lapply(corr_pairs, function(pair) {
    make_corr_hex_plot(d_plot[sampled_idx, ], pair[1], pair[2])
})

pdf("./figures/232_gene_score_correlation_all_genes.pdf", width = 18, height = 12)
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]] | plots[[5]] | plots[[6]]) /
(plot_spacer()| plots[[7]] | plots[[8]] | plots[[9]]| plots[[10]] | plots[[11]]) /
(plot_spacer()| plot_spacer()| plots[[12]] | plots[[13]] | plots[[14]]| plots[[15]]) /
(plot_spacer()| plot_spacer()| plot_spacer()| plots[[16]] | plots[[17]] | plots[[18]]) /
(plot_spacer()| plot_spacer()| plot_spacer()| plot_spacer()| plots[[19]] | plots[[20]]) /
(plot_spacer()| plot_spacer()| plot_spacer()| plot_spacer()| plot_spacer()| plots[[21]])
dev.off()
## }}}

## Gene score correlation across cpg genes {{{
cpg_genes = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")
plots <- lapply(corr_pairs, function(pair) {
    make_corr_hex_plot(d_plot[sampled_idx, ][gene %in% cpg_genes], pair[1], pair[2])
})

pdf("./figures/232_gene_score_correlation_cpg_genes.pdf", width = 18, height = 12)
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]] | plots[[5]] | plots[[6]]) /
(plot_spacer()| plots[[7]] | plots[[8]] | plots[[9]]| plots[[10]] | plots[[11]]) /
(plot_spacer()| plot_spacer()| plots[[12]] | plots[[13]] | plots[[14]]| plots[[15]]) /
(plot_spacer()| plot_spacer()| plot_spacer()| plots[[16]] | plots[[17]] | plots[[18]]) /
(plot_spacer()| plot_spacer()| plot_spacer()| plot_spacer()| plots[[19]] | plots[[20]]) /
(plot_spacer()| plot_spacer()| plot_spacer()| plot_spacer()| plot_spacer()| plots[[21]])
dev.off()
## }}}


## Correlation matrix for all markers {{{
library(ellipse)
d_cor = cor(d_plot[, -(1:2)], method = "spearman")
ord = names(m_list)
d_cor = d_cor[ord, ord]
library(RColorBrewer)
my_colors <- brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)[100:1]

pdf("./figures/232_gene_score_correlation_marker_correlation.pdf", width = 6, height = 6)
plotcorr(d_cor, col = my_colors[cut(d_cor, breaks = 100)], main = "Marker correlation (Spearman)")
## generate the color key
plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
legend_image <- as.raster(matrix(my_colors, ncol = 1))
rasterImage(legend_image, 0.1, 0.1, 0.2, 0.9)
text(x = 0.3, y = seq(0.1, 0.9, length.out = 5), labels = round(seq(-1, 1, length.out = 5), 2), cex = 0.8)
dev.off()
## }}}

## back {{{
# library(GGally)
# ggpairs(
#     d_plot[sampled_idx, ] + 0.1,
#     columns = 3:ncol(d_plot),
#     upper = list(continuous = wrap("cor", size = 3)),
#     lower = list(continuous = wrap("points", size = 0.1))
# ) + scale_x_log10() + scale_y_log10()
## }}}
