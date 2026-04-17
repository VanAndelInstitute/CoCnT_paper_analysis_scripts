library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(Cairo)

library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 11)
color = ArchRPalettes$stallion

httpgd::hgd(port = 4322)

set.seed(123)

proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_HSPC/")

d_meta_k27 = as.data.table(proj_k27@cellColData)
umap_k27 = as.data.table( getEmbedding( ArchRProj = proj_k27, embedding = "UMAP_Harmony"))
names(umap_k27) = c("UMAP1", "UMAP2")
d_meta = cbind(d_meta_k27, umap_k27)

## draw umap using ArchR color shceme
color_cluster = ArchRPalettes$stallion
names(color_cluster) = c(
	"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "Other_cells"
)
ggplot(d_meta, aes(x = UMAP1, y = UMAP2, color = Clusters)) +
    geom_point(alpha = 0.5, size = 0.3) +
    scale_color_manual(values = color_cluster) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.title = element_blank()
	) +
    ggtitle("UMAP of HSPC - colored by clusters")
ggsave("./figures/205_UMAP_HSPC_clusters.pdf", width = 5, height = 4)


## imputed gene score
m_h3k4me1 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me1_from_h3k27me3_HSPC.rds")
m_h3k4me2 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me2_from_h3k27me3_HSPC.rds")
m_h3k4me3 = read_rds("./tmp/matrix_imputed_gene_score_h3k4me3_from_h3k27me3_HSPC.rds")
m_h3k27 = read_rds("./tmp/matrix_imputed_gene_score_h3k27me3_HSPC.rds")

colnames(m_h3k4me1) %<>% sub(".*#", "", .)
colnames(m_h3k4me2) %<>% sub(".*#", "", .)
colnames(m_h3k4me3) %<>% sub(".*#", "", .)
colnames(m_h3k27) %<>% sub(".*#", "", .)

m_list = list(
    m_h3k4me1 = m_h3k4me1,
    m_h3k4me2 = m_h3k4me2,
    m_h3k4me3 = m_h3k4me3,
    m_h3k27 = m_h3k27
)

names(m_list) = c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27me3")

write_rds(m_list, "./tmp/list_matrix_imputed_gene_score_HSPC.rds")

m_list = read_rds("./tmp/list_matrix_imputed_gene_score_HSPC.rds")

color = c(
    H3K27me3 = '#3B8FC4', 
    H3K4me1 = "#8c8023",
    H3K4me1 = '#FFD000',
    H3K4me2 = '#F39C12', 
    H3K4me3 = '#E74C3C'
)

m_list_scale = lapply(m_list, function(m) {
	 m_scaled = t(scale(t(m)))
	 rownames(m_scaled) = rownames(m)
	 colnames(m_scaled) = colnames(m)
	 return(m_scaled)
 })

names(m_list_scale) = c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27me3") 

write_rds(m_list_scale, "./tmp/list_matrix_imputed_gene_score_scaled_HSPC.rds")
m_list_scale = read_rds("./tmp/list_matrix_imputed_gene_score_scaled_HSPC.rds")

d_meta$Clusters %<>% factor(levels = c("C1", "C2", "C6", "C7", "C8", "C3", "C5", "C4"))

plot_umap_gene_score = function(gene_name = "MECOM", m_list, color_v) {
    lapply(1:length(m_list), function(i) {
	color_i = color_v[names(m_list)[i]]
	m = m_list[[i]]
	gene_score_v = m[gene_name, ]

	## ceil the top by 10 percent by gene
	quants = quantile(gene_score_v, probs = c(0.95), na.rm = T)
	gene_score_v[gene_score_v > quants[1]] = quants[1]

	d_plot = d_meta[cell_id %in% names(gene_score_v)]
	d_plot$gene_score = gene_score_v[match(d_plot$cell_id, names(gene_score_v))]

	ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = gene_score)) +
	    geom_point(size = 0.2, shape = '.') +
	    scale_color_gradient(low = "white", high = color_i) +
	    theme_classic() +
	    theme(
		legend.position = "right",
		axis.title = element_blank()
		) +
	    ggtitle(str_glue("{names(m_list)[i]} - {gene_name}"))
    }) %>% ggarrange(plotlist = ., ncol = 2, nrow = 4)
}

Cairo::CairoPDF("./figures/205_imputed_umap_gene_score_example_HSPC.pdf", width = 9, height = 14)
gene_list = c("MECOM", "PAX5", "CEBPA", "GATA1", "CD34", "SPI1", "RUNX1", "KLF1", "EGR1", "IRF8")
for (gene_name in gene_list) {
    print(gene_name)
    print(plot_umap_gene_score(gene_name, m_list, color) )
}
dev.off()

plot_boxplot_gene_score = function(gene_name = "MECOM", m_list, color_v) {
    lapply(1:length(m_list), function(i) {
	color_i = color_v[names(m_list)[i]]
	m = m_list[[i]]
	gene_score_v = m[gene_name, ]

	d_plot = d_meta[cell_id %in% names(gene_score_v)]
	d_plot$gene_score = gene_score_v[match(d_plot$cell_id, names(gene_score_v))]

	ggplot(d_plot, aes(x = Clusters, y = gene_score, fill = Clusters)) +
	    geom_boxplot() +
	    scale_fill_manual(values = color_cluster) +
	    theme_classic() +
	    theme(
		legend.position = "none",
		axis.title.x = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1)
		) +
	    ggtitle(str_glue("{names(m_list)[i]} - {gene_name}"))
    }) %>% ggarrange(plotlist = ., ncol = 2, nrow = 4)
}

Cairo::CairoPDF("./figures/205_imputed_boxplot_gene_score_example_HSPC.pdf", width = 9, height = 14)
gene_list = c("MECOM", "PAX5", "CEBPA", "GATA1", "CD34", "SPI1", "RUNX1", "KLF1", "EGR1", "IRF8", "MPO", "ELANE")
for (gene_name in gene_list) {
	print(gene_name)
	print(plot_boxplot_gene_score(gene_name, m_list, color) )
}
dev.off()

