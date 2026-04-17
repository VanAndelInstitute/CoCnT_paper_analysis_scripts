library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

httpgd::hgd(port = 4322)

d_me1 = fread("./tmp/table_gene_chromatin_status_wide_me1.tsv") %>% melt(id.vars = c("gene"), variable.name = "cell_type", value.name = "cat")
d_me2 = fread("./tmp/table_gene_chromatin_status_wide_me2.tsv") %>% melt(id.vars = c("gene"), variable.name = "cell_type", value.name = "cat")
d_me3 = fread("./tmp/table_gene_chromatin_status_wide_me3.tsv") %>% melt(id.vars = c("gene"), variable.name = "cell_type", value.name = "cat")

d = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype.tsv")

d_exp = fread("./tmp/table_triana_2021_rna_pseudo_bulk_by_ct.tsv")
names(d_exp)

d$consensus = factor(d$consensus, levels = c("me1", "me2", "me3", "me1+me2", "me1+me3", "me2+me3", "me1+me2+me3"))

## Shared active is higher expressed than unique H3K4me3 bivalent {{{

## The good quality active genes should show higher expression than active support by one marker

d_plot = d[cat == "Active" & cell_type == "HSPC"]
d_plot = merge(d_plot, d_exp[, .(gene, exp = `HSCs & MPPs`)], by = "gene")

d_plot[, .N, by = consensus][N > 100]

comparison = list(
    c("1", "2"),
    c("1", "3"),
    c("2", "3")
)

pdf("./figures/224_chromatin_status_active_exp_by_consensus_count.pdf", width = 6, height = 4)
ggplot(d_plot, aes(x = factor(consensus_count), y = exp)) +
	geom_boxplot(outliers = F)  +
	stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = comparison) +
	scale_y_log10() + 
	theme_classic()

## two scale plot put the gene number on the left and the expression on the right
d_plot_gene_count = d_plot[, .N, by = consensus_count]
d_plot_gene_count$N = cumsum(d_plot_gene_count$N)

ggplot(d_plot_gene_count) + aes(x = factor(consensus_count), y = N) +
	geom_point() +
	geom_line(group = 1) +
	theme_classic() +
	## axis on the right
	scale_y_continuous(sec.axis = sec_axis(~ ., name = "Cumulative gene count")) +
	ylab("Gene count") +
	ggtitle("Gene count by number of supporting marks for active genes in HSPC")
dev.off()

## }}}

## Shared bivalent is more related to development than unique bivalent {{{

# Get all genes in the Go term "developmental process" (GO:0032502) 
library(AnnotationDbi)
library(org.Hs.eg.db)
go_term_v = c("GO:0007275", "GO:0032502", "GO:0044238", "GO:0010817", "GO:0030097", "GO:0045165", "GO:0048568",
    "GO:0008152", "GO:0050435")

pdf("./figures/224_chromatin_status_vs_percent_in_go_term_by_consensus_count.pdf", width = 6, height = 4)

for (i in 1:length(go_term_v)) {
    ## show key types in org.Hs.eg.db
    go_term = go_term_v[i]
    keytypes(org.Hs.eg.db)
    gene_in_go = AnnotationDbi::select(
	org.Hs.eg.db,
	keys = go_term,
	keytype = "GOALL",
	columns = c("SYMBOL")
	)$SYMBOL
    gene_in_go = unique(na.omit(gene_in_go))
    ## 11516 genes
    gene_in_go %>% length

    ## bivalent genes in each marker
    d_plot = d
    d_plot[, is_in_go_term := ifelse(gene %in% gene_in_go, "In GO term", "Not in GO term")]
    d_plot = d_plot
    # percentage of genes in GO term by consensus count
    d_plot = d_plot[, .N, by = .(consensus_count > 1, is_in_go_term, cell_type, cat)]
    d_plot[, prop := N / sum(N), by = .(consensus_count, cell_type, cat)]
    d_plot = d_plot[is_in_go_term == "In GO term"]

    g = ggplot(d_plot[cat != "Other"], aes(x = factor(consensus_count), y = prop)) +
	geom_boxplot(outliers = F) +
	geom_line(aes(group = cell_type)) +
	geom_point() +
	facet_grid(~ cat) +
	theme_classic() + 
	labs(title = paste0("Proportion of genes in GO term ", go_term, " by number of supporting marks for bivalent genes"), y = "Proportion of genes in GO term")
    print(g)
}
dev.off()

## }}}

## Shared bivalent is more related to conservative between species than unique bivalent {{{

## }}}

## Shared bivalent is more related to the presence of CpG islands than unique bivalent {{{
cpg_genes = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")

# 1, 2, 3 markers. Bivalent genes cpg %

pdf("./figures/224_chromatin_status_vs_percent_cpg_by_consensus_count.pdf", width = 24, height = 4)
# calcualte the percentage of CpG promoter genes by consensus count
d_plot = copy(d)
d_plot[, is_cpg := ifelse(gene %in% cpg_genes, "CpG promoter", "Non-CpG promoter")]
d_plot = d_plot[, .N, by = .(consensus_count, is_cpg, cell_type, cat)]
d_plot = d_plot[, .(prop = N / sum(N), is_cpg), by = .(consensus_count, cell_type, cat)][is_cpg == "CpG promoter"]
ggplot(d_plot[cat != "Other"], aes(x = factor(consensus_count), y = prop)) +
	geom_boxplot(outliers = F) +
	geom_line(aes(group = cell_type)) +
	geom_point() +
	facet_grid(~ cat) +
	theme_classic() +
	labs(title = "Proportion of genes with CpG promoter by number of supporting marks for bivalent genes", y = "Proportion of genes with CpG promoter")

# calcualte the percentage of CpG promoter genes by consensus
d_plot = copy(d)
d_plot[, is_cpg := ifelse(gene %in% cpg_genes, "CpG promoter", "Non-CpG promoter")]
d_plot = d_plot[, .N, by = .(consensus, is_cpg, cell_type, cat)]
d_plot = d_plot[, .(prop = N / sum(N), is_cpg), by = .(consensus, cell_type, cat)][is_cpg == "CpG promoter"]
ggplot(d_plot[cat != "Other"], aes(x = factor(consensus), y = prop)) +
	geom_boxplot(outliers = F) +
	geom_line(aes(group = cell_type)) +
	geom_point() +
	facet_grid(~ cat) +
	theme_classic() +
	labs(title = "Proportion of genes with CpG promoter by number of supporting marks for bivalent genes", y = "Proportion of genes with CpG promoter")

dev.off()



## }}}

## Shared bivalent is more related to the presence of the super enhancers than unique bivalent {{{

## }}}

## Shared bivalent is more related to the BMI1 KO phenotype than unique bivalent {{{
x = fread("../../data/Gene_Lists/E_Biv_wE_Enhancers.csv")
d_bmi_io = fread("./tmp/table_BMI1_KO_gene_log2FC_summary.tsv")[p_g1 < 0.25 & p_g3 < 0.25]

d_plot = merge(d[cell_type == "Erythroid progenitors1", .(gene, cat, consensus_count)], d_bmi_io, by.x = "gene", by.y = "SYMBOL")[cat != "Other"]
d_plot$cat %>% table
# d_plot = merge(d_me2[cell_type == "Erythroid progenitors3", .(gene, cat)], d_bmi_io, by.x = "gene", by.y = "SYMBOL")[cat != "Other"]

pdf("./figures/224_chromatin_status_bmi1_ko_by_cat.pdf", width = 6, height = 4)

ggplot(d_plot) + aes(x = cat, y = (mean_log2FC_g1 + mean_log2FC_g3) / 2, fill = cat) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Active", label.y = 1) +
    facet_wrap(~consensus_count) +
    scale_fill_nejm() +
    theme_classic() + 
    labs(title = "Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby chromatin status in Erythroid progenitors", y = "Mean log2FC BMI1 g1 KO vs Neg2")

ggplot(d_plot) + aes(x = cat, y = mean_log2FC_g1, fill = cat) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Bivalent", label.y = 1) +
    scale_fill_nejm() +
    theme_classic() +
    labs(title = "Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby chromatin status in Erythroid progenitors", y = "Mean log2FC BMI1 g1 KO vs Neg2")

ggplot(d_plot[BMI_g3 > Neg2 & BMI_g1 > Neg2][cat %in% c("Bivalent", "Repressive")]) + aes(x = cat, y = BMI_g3 - Neg2 , fill = cat) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Bivalent", label.y = 1) +
    facet_wrap(~consensus_count) +
    scale_fill_nejm() +
    theme_classic() +
    labs(title = "Gene expression changes upon BMI1 g3 KO vs Neg2 in CD34+ cells\nby chromatin status in Erythroid progenitors", y = "BMI1 g3 KO - Neg2")

ggplot(d_plot[BMI_g3 > Neg2 & BMI_g1 > Neg2][cat %in% c("Bivalent", "Repressive")]) + aes(x = cat, y = BMI_g1 - Neg2 , fill = cat) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Bivalent", label.y = 1) +
    facet_wrap(~consensus_count) +
    scale_fill_nejm() +
    theme_classic() + 
    labs(title = "Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby chromatin status in Erythroid progenitors", y = "BMI1 g1 KO - Neg2")

ggplot(d_plot[BMI_g3 > Neg2 & BMI_g1 > Neg2][cat %in% c("Bivalent", "Repressive")]) + aes(x = cat, y = (BMI_g1 + BMI_g3 - Neg2 * 2), fill = cat) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Bivalent", label.y = 1) +
    facet_wrap(~consensus_count) +
    scale_fill_nejm() +
    theme_classic() + 
    labs(title = "Gene expression changes upon BMI1 g1 and g3 KO vs Neg2 in CD34+ cells\nby chromatin status in Erythroid progenitors", y = "BMI1 g1 and g3 KO - Neg2")
dev.off()

## }}}


############
#  Backup  #
############

## E_Biv_wE_enahncers Positive control {{{
## E enhancer number
d_cat = d[consensus_count > 1 & celltype == "HSPC", .(gene, cat)]
x = fread("../../data/Gene_Lists/E_Enhancers_All_Genes.csv")
# d_bmi_io = fread("./tmp/table_BMI1_KO_gene_log2FC_summary.tsv")
# d_bmi_io = fread("./tmp/limma_EZH2_g1_vs_Neg2.tsv")
# d_bmi_io = fread("../../data/Gene_Lists/sgEZH2_DE.csv")
d_bmi_io = fread("../../data/Gene_Lists/sgBMI1_2_DE_Genes_Up.csv")
# d_bmi_io[, SYMBOL := gene]
# d_bmi_io[, mean_log2FC_g1 := log2FoldChange]

d_plot = merge(x, d_bmi_io, by.x = "gene", by.y = "gene", all.x = T)
d_plot = merge(d_plot, d_cat, by = "gene")
d_plot = d_plot[cat != "Other"]
d_plot$cat %>% table
# d_plot = d_plot[BMI_g1 > Neg2]

# v_cpg = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")

# d_plot[, is_cpg := ifelse(gene %in% v_cpg, "CpG promoter", "Non-CpG promoter")]

## how to get gene body size for all genes
# library(GenomicFeatures)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#
# gr_gene <- genes(txdb)
#
# gene_size_dt <- data.frame(
#   gene_id = names(gr_gene),
#   seqnames = as.character(seqnames(gr_gene)),
#   start = start(gr_gene),
#   end = end(gr_gene),
#   strand = as.character(strand(gr_gene)),
#   gene_body_size = width(gr_gene)
# )
#
# library(AnnotationDbi)
# library(org.Hs.eg.db)
#
# gene_anno <- AnnotationDbi::select(
#   org.Hs.eg.db,
#   keys = gene_size_dt$gene_id,
#   keytype = "ENTREZID",
#   columns = c("SYMBOL", "GENENAME")
# )
#
# gene_size_dt2 <- merge(
#   gene_size_dt,
#   gene_anno,
#   by.x = "gene_id",
#   by.y = "ENTREZID",
#   all.x = TRUE
# )

# head(gene_size_dt2)

# d_plot = merge(d_plot, data.table(gene_size_dt2)[, .(SYMBOL, gene_body_size)], by.x = "gene", by.y = "SYMBOL")


ggplot(d_plot, aes(x = cat, y = Total_E_Enhancer, fill = cat)) +
	geom_boxplot(outliers = F) +
	scale_fill_nejm() +
	# facet_wrap(~is_cpg) +
	theme_classic() +
	stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Active", label.y = 1) +
	labs(title = "Number of E_Enhancers by chromatin status in Erythroid progenitors", y = "Number of E_Enhancers")

d_plot[is.na(log2FoldChange), log2FoldChange := 0]
d_plot = d_plot[log2FoldChange >= 0]

table(d_plot$log2FoldChange == 0, d_plot$cat)

ggplot(d_plot, aes(x = log2FoldChange > 0, y = Total_E_Enhancer, fill = cat)) +
	geom_boxplot(outliers = F) +
	# scale_fill_nejm() +
	# facet_wrap(~is_cpg) +
	theme_classic() +
	stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "Active", label.y = 1) +
	labs(title = "Number of E_Enhancers by chromatin status in Erythroid progenitors", y = "Number of E_Enhancers")
#
#
# ggplot(d_plot, aes(x = Total_E_Enhancer, y = gene_body_size, fill = cat)) +
#     geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.5) +
#     stat_cor(method = "spearman") +
#     theme_classic() +
#     facet_wrap(~cat) +
# 	labs(title = "Gene body size by number of E_Enhancers in Erythroid progenitors", y = "Gene body size")


pdf("./figures/224_chromatin_status_bmi1_ko_by_e_enhancer.pdf", width = 10, height = 6)

d_plot[, y := Total_E_Enhancer]

d_plot[, .N, by = y]
d_plot[y > 4, y := 4]
table(d_plot$y)

d_plot[gene == "GCNT2"]


ggplot(d_plot) + aes(x = factor(y), y = logFC) + #, group = cat) +
    # geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.2) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    # stat_cor(method = "spearman", label.y = 15) +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "0", label.y = 1) +
    theme_classic() +
    facet_wrap(~cat, nrow = 1) +
    ggtitle("Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby number of Total E_Enhancers")

ggplot(d_plot) + aes(x = factor(y), y = BMI_g1 - Neg2) + #, group = cat) +
    # geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.2) +
    geom_boxplot(outliers = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    # stat_cor(method = "spearman", label.y = 15) +
    stat_compare_means(method = "wilcox.test", label = "p.format", ref.group = "0", label.y = 1) +
    theme_classic() +
    facet_wrap(~cat, nrow = 1) +
    ggtitle("Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby number of Total E_Enhancers")



d_plot[, y := Strong_E_Enhancer]

d_plot[, .N, by = y]
table(d_plot$y)

ggplot(d_plot[y < 10]) + aes(x = y, y = mean_log2FC_g1, group = cat) +
    geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.2) +
    stat_cor(method = "spearman", label.y = 15) +
    theme_classic() +
    facet_wrap(~cat, nrow = 1, scales = "free_y") +
    ggtitle("Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby number of Strong E_Enhancers")


d_plot[, y := E_Enhancer]

d_plot[, .N, by = y]
table(d_plot$y)

ggplot(d_plot[y < 10]) + aes(x = y, y = mean_log2FC_g1, group = cat) +
    geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.2) +
    stat_cor(method = "spearman", label.y = 15) +
    theme_classic() +
    facet_wrap(~cat, nrow = 1, scales = "free_y") +
    ggtitle("Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby number of E_Enhancers")

dev.off()


ggplot(d_plot) + aes(x = factor(y), y = (mean_log2FC_g1 + mean_log2FC_g3) / 2, fill = factor(y)) +
    geom_boxplot(outliers = F) +
    stat_compare_means(method = "t.test", label = "p.format", ref.group = "0", label.y = 1) +
    ggtitle("Gene expression changes upon BMI1 g1 KO vs Neg2 in CD34+ cells\nby presence of E_Enhancers") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    scale_fill_nejm() +
    theme_classic()

}}}
