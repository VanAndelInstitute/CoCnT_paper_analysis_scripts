library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggpubr)
library(ggsci)

library(ArchR)

addArchRGenome("hg38")
addArchRThreads(threads = 8)

httpgd::hgd(port = 4323)

proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_HSPC/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_HSPC/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_HSPC/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_HSPC/")

d_k27 = as.data.table(proj_k27@cellColData)
d_h3k4me1 = as.data.table(proj_h3k4me1@cellColData)
d_h3k4me2 = as.data.table(proj_h3k4me2@cellColData)
d_h3k4me3 = as.data.table(proj_h3k4me3@cellColData)


d_k27$marker = "H3K27me3"
d_h3k4me1$marker = "H3K4me1"
d_h3k4me2$marker = "H3K4me2"
d_h3k4me3$marker = "H3K4me3"

d_all = rbind(
	d_k27,
	d_h3k4me1,
	d_h3k4me2,
	d_h3k4me3
)

## get UMAP
umap_k27 = as.data.table( getEmbedding( ArchRProj = proj_k27, embedding = "UMAP"))
names(umap_k27) = c("UMAP1", "UMAP2")
umap_k27$cell_id = d_k27$cell_id

d_plot = merge(d_all, umap_k27, by = "cell_id", all.x = TRUE)

d_plot$Clusters_HSPC = factor(d_plot$Clusters_HSPC, levels = c(
	"C1", "C2", "C6", "C7", "C8", "C3", "C5", "C4"
))

## Show nFrags on UMAP
p_list = lapply(
	unique(d_plot$marker),
	function(x) {
	d_sub = d_plot[marker == x]
	p = ggplot(d_sub, aes(x = UMAP1, y = UMAP2, color = log10(nFrags))) +
	    geom_point(size = 0.5, alpha = 0.5) +
	    scale_color_viridis_c() +
	    theme_classic() +
	    labs(
		title = paste0(x, ": nFrags on UMAP (H3K27me3 UMAP)"),
		x = "UMAP1",
		y = "UMAP2",
		color = "nFrags"
	    )
	return(p)
	}
)

ggarrange(plotlist = p_list, ncol = 2, nrow = 2)


## show nFrags with boxplot by cluster
p_list = lapply(
	unique(d_plot$marker),
	function(x) {
	d_sub = d_plot[marker == x]
	p = ggplot(d_sub, aes(x = Clusters_HSPC, y = nFrags)) +
	    geom_boxplot(outlier.size = 0.5, alpha = 0.7, outliers = F) +
	    theme_classic() +
	    labs(
		title = paste0(x, ": nFrags by cluster"),
		x = "Cluster",
		y = "nFrags"
	    ) +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	    guides(fill = guide_legend(override.aes = list(size=2, alpha = 1)))
	return(p)
	}
)
ggarrange(plotlist = p_list, ncol = 2, nrow = 2)

## Show donor assignment on UMAP
p_list = lapply(
	unique(d_plot$marker),
	function(x) {
	d_sub = d_plot[marker == x & assignment %in% c(0, 1, 2, 3)]
	p = ggplot(d_sub, aes(x = UMAP1, y = UMAP2, color = as.factor(assignment))) +
	    geom_point(size = 0.5, alpha = 0.5) +
	    scale_color_nejm() +
	    theme_classic() +
	    labs(
		title = paste0(x, ": donor assignment on UMAP (H3K27me3 UMAP)"),
		x = "UMAP1",
		y = "UMAP2",
		color = "donor assignment"
	    )
	return(p)
	}
)

ggarrange(plotlist = p_list, ncol = 2, nrow = 2)

## show donor assignment with boxplot by cluster
p_list = lapply(
	unique(d_plot$marker),
	function(x) {
	d_sub = d_plot[marker == x & assignment %in% c(0, 1, 2, 3)]
	p = ggplot(d_sub, aes(x = Clusters_HSPC, fill = as.factor(assignment))) +
	    geom_bar(position = "fill") +
	    scale_fill_nejm() +
	    theme_classic() +
	    labs(
		title = paste0(x, ": donor assignment by cluster"),
		x = "Cluster",
		y = "Proportion",
		fill = "donor assignment"
	    ) +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	    guides(fill = guide_legend(override.aes = list(size=2, alpha = 1)))
	return(p)
	}
)
ggarrange(plotlist = p_list, ncol = 2, nrow = 2)

d_plot[, .N, by = .(marker, Clusters_HSPC, assignment)][assignment %in% c(1, 2, 3, 4)][Clusters_HSPC %in% c("C1", "C6")]
d_plot[, .(singlet_rate = sum(status == "singlet") / .N), by = .(marker, Clusters_HSPC)] %>%
	dcast(marker ~ Clusters_HSPC, value.var = "singlet_rate")


## show antibody assignment on UMAP
p_list = lapply(
	unique(d_plot$marker),
	function(x) {
	d_sub = d_plot[marker == x]
	p = ggplot(d_sub, aes(x = UMAP1, y = UMAP2, color = as.factor(antibody_species))) +
	    geom_point(size = 0.5, alpha = 0.5) +
	    scale_color_nejm() +
	    theme_classic() +
	    labs(
		title = paste0(x, ": antibody species on UMAP (H3K27me3 UMAP)"),
		x = "UMAP1",
		y = "UMAP2",
		color = "antibody species"
	    )
	return(p)
	}
)
ggarrange(plotlist = p_list, ncol = 2, nrow = 2)

## show antibody species with boxplot by cluster
p_list = lapply(
	unique(d_plot$marker),
	function(x) {
	d_sub = d_plot[marker == x]
	p = ggplot(d_sub, aes(x = Clusters_HSPC, fill = as.factor(antibody_species))) +
	    geom_bar(position = "fill") +
	    scale_fill_nejm() +
	    theme_classic() +
	    labs(
		title = paste0(x, ": antibody species by cluster"),
		x = "Cluster",
		y = "Proportion",
		fill = "antibody species"
	    ) +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	    guides(fill = guide_legend(override.aes = list(size=2, alpha = 1)))
	return(p)
	}
)
ggarrange(plotlist = p_list, ncol = 2, nrow = 2)


