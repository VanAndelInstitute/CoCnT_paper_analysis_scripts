library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d = fread("./tmp/table_triana_2021_rna_pseudo_bulk_by_ct.tsv")
d_transition = fread("./tmp/table_shared_genes_between_markers_in_trajectory.tsv")


celltype_order <- c(
    # Stem & multipotent progenitors
    "HSCs & MPPs",

    # B cell lineage
    # "Lymphomyeloid prog",
    # "Pre-pro-B cells",
    # "Pro-B cells",
    # "Pre-B cells",
    # "Small pre-B cell",
    # "Immature B cells",
    # "Mature naive B cells",
    # "Nonswitched memory B cells",
    # "Class switched memory B cells",
    # "Plasma cells"

    # Erythroid / Megakaryocytic lineage
    "Erythro-myeloid progenitors",
    "Early erythroid progenitor",
    "Late erythroid progenitor"
)

cell_type_order_l = list(
    "E" = c(
	"HSCs & MPPs",
	"Erythro-myeloid progenitors",
	"Early erythroid progenitor",
	"Late erythroid progenitor"
	),
    "B" = c(
	"HSCs & MPPs",
	"Lymphomyeloid prog",
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
)


pdf("./figures/237_active2repressive_scRNASeq.pdf", p, width = 7, height = 5)

for (lineage_x in c("B", "E")) {
    celltype_order = cell_type_order_l[[lineage_x]]
    for (marker_x in c("me1", "me2", "me3")) {

	d_tran_sub = d_transition[
	    trans %in% c("Active->Un->Repressive", "Active->Bivalent->Repressive", "Active->Repressive") & 
		grepl(marker_x, marker) & lineage %in% c(lineage_x), .(lineage, trans, gene)]

	    d_long = melt(d)

	    d_plot = merge(d_long, d_tran_sub, by = "gene")
	    d_plot$variable = factor(d_plot$variable, levels = rev(celltype_order))
	    d_plot = na.omit(d_plot)
	    # d_plot = d_plot[, .(median = median(value)), by = .(lineage, trans, variable)]

	    p = ggplot(d_plot, aes(x = variable, y = log1p(value), color = trans, group = variable)) +
		geom_boxplot(outliers = F) +
		# geom_line() +
		theme_classic() +
		scale_color_manual(values = 
		    c("Active->Bivalent->Repressive" = "#ecb682", "Active->Un->Repressive" = "#a4cce5", "Active->Repressive" = "#d88c80")) +
		facet_wrap(~trans, nrow = 1, scales = "free_x") +
		theme(
		    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		) + labs(x = "Cell type", y = "Median expression") + ggtitle(paste0(lineage_x, " lineage - ", marker_x)) +
		coord_flip()
	    
		
		print(p)
    }
}
dev.off()


## Gene ontology analysis {{{
# library(clusterProfiler)
# library(org.Hs.eg.db)
#
# genes = d_tran_sub[lineage == "B" & trans == "Active->Un->Repressive", gene]
# genes_entrez = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# genes_entrez[genes_entrez$SYMBOL == "MPO", ]
# ego <- enrichGO(
#     gene 	  = genes_entrez$ENTREZID,
#     OrgDb         = org.Hs.eg.db,
#     keyType       = 'ENTREZID',
#     ont           = "BP",
#     pAdjustMethod = "BH",
#     pvalueCutoff  = 0.01,
#     qvalueCutoff  = 0.05
# )
#
#
# genes = d_tran_sub[lineage == "E" & trans == "Active->Bivalent->Repressive", gene]
# genes_entrez = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego <- enrichGO(
#     gene 	  = genes_entrez$ENTREZID,
#     OrgDb         = org.Hs.eg.db,
#     keyType       = 'ENTREZID',
#     ont           = "BP",
#     pAdjustMethod = "BH",
#     pvalueCutoff  = 0.05,
#     qvalueCutoff  = 0.55
# )
#
# ego_result_abr = as.data.table(ego) }}}

