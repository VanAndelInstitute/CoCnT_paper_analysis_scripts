library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

httpgd::hgd(port = 4322)

# The data needed
# bulk RNA seq data
d_rna = fread("./tmp/table_triana_2021_rna_pseudo_bulk_by_ct.tsv")
# cpg data
cpg_genes = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")
# bulk HSPC chromatin state
d_state = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")
d_state_long = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype.tsv")
# bulk HSPC histone gene score
d_score = fread("./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")

## gene expression v.s. chromatin state {{{
pdf("./figures/243_chromatin_state_vs_rna_expression.pdf", width = 4, height = 4)
comparisons = list(
	c("Active", "Bivalent"),
	c("Active", "Repressive"),
	c("Active", "Un"),
	c("Bivalent", "Repressive"),
	c("Bivalent", "Un"),
	c("Repressive", "Un")
)

d_plot = merge(
    d_rna[, .(gene, rna_expression = `HSCs & MPPs`)],
    d_state[, .(gene, chromatin_state = `HSPC`)]
    )[chromatin_state %in% c("Active", "Bivalent", "Repressive", "Un")]

ggplot(d_plot, aes(x = chromatin_state, y = rna_expression)) +
	geom_boxplot(outliers = F) +
	theme_classic() +
	labs(x = "Chromatin state", y = "RNA expression (HSPC)") +
	stat_compare_means(method = "wilcox.test", comparisons = comparisons, p.adjust.method = "BH") +
	scale_y_log10() + labs(title = "RNA expression vs. consensus chromatin state (HSPC)")

# H3K4me1
d_plot = merge(
	d_rna[, .(gene, rna_expression = `HSCs & MPPs`)],
	d_state_long[cell_type == "HSPC"][grep("me1", consensus), .(gene, chromatin_state = cat)]
	)[chromatin_state %in% c("Active", "Bivalent", "Repressive", "Un")]
ggplot(d_plot, aes(x = chromatin_state, y = rna_expression)) +
	geom_boxplot(outliers = F) +
	theme_classic() +
	labs(x = "Chromatin state", y = "RNA expression (HSPC)") +
	stat_compare_means(method = "wilcox.test", comparisons = comparisons, p.adjust.method = "BH") +
	scale_y_log10() + labs(title = "RNA expression vs. H3K4me1 chromatin state (HSPC)")

# H3K4me2
d_plot = merge(
	d_rna[, .(gene, rna_expression = `HSCs & MPPs`)],
	d_state_long[cell_type == "HSPC"][grep("me2", consensus), .(gene, chromatin_state = cat)]
	)[chromatin_state %in% c("Active", "Bivalent", "Repressive", "Un")]
ggplot(d_plot, aes(x = chromatin_state, y = rna_expression)) +
	geom_boxplot(outliers = F) +
	theme_classic() +
	labs(x = "Chromatin state", y = "RNA expression (HSPC)") +
	stat_compare_means(method = "wilcox.test", comparisons = comparisons, p.adjust.method = "BH") +
	scale_y_log10() + labs(title = "RNA expression vs. H3K4me2 chromatin state (HSPC)")

# H3K4me3
d_plot = merge(
	d_rna[, .(gene, rna_expression = `HSCs & MPPs`)],
	d_state_long[cell_type == "HSPC"][grep("me3", consensus), .(gene, chromatin_state = cat)]
	)[chromatin_state %in% c("Active", "Bivalent", "Repressive", "Un")]
ggplot(d_plot, aes(x = chromatin_state, y = rna_expression)) +
	geom_boxplot(outliers = F) +
	theme_classic() +
	labs(x = "Chromatin state", y = "RNA expression (HSPC)") +
	stat_compare_means(method = "wilcox.test", comparisons = comparisons, p.adjust.method = "BH") +
	scale_y_log10() + labs(title = "RNA expression vs. H3K4me3 chromatin state (HSPC)")

dev.off()
## }}}

## gene expression v.s. histone gene score all {{{

d_plot = merge(
    d_rna[, .(gene, rna_expression = `HSCs & MPPs`)],
    d_score[cell_type == "HSPC"][, .(gene, marker, gene_score, is_cpg = gene %in% cpg_genes)]
)

pdf("./figures/243_histone_gene_score_vs_rna_expression.pdf", width = 12, height = 10)

ggplot(d_plot, aes(x = gene_score + 0.1, y = rna_expression + 0.1)) +
    geom_hex(bins = 100) +
    theme_classic() +
    facet_wrap(~marker) +
    stat_cor(method = "spearman") +
    scale_fill_viridis_c(trans = "log10", name = "Cell count") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Histone gene score (HSPC)", y = "RNA expression (HSPC)")

ggplot(d_plot[is_cpg == TRUE], aes(x = gene_score + 0.1, y = rna_expression + 0.1)) +
    geom_hex(bins = 100) +
    theme_classic() +
    facet_wrap(~marker) +
    stat_cor(method = "spearman") +
    scale_fill_viridis_c(trans = "log10", name = "Cell count") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Histone gene score (HSPC CpG genes)", y = "RNA expression (HSPC)")

dev.off()

## }}}

## Un->Active, and Bivalent->Active {{{
d_stat_trans = fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")
gene_un_2active = d_stat_trans[trans %in% c("Un->Active") & group == "B_consensus", gene]
gene_bi_2active = d_stat_trans[trans %in% c("Bivalent->Active") & group == "B_consensus", gene]

celltype_order <- c(
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
)

## }}}
