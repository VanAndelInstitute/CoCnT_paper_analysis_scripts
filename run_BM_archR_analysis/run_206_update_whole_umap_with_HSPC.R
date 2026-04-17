library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 11)

color = ArchRPalettes$stallion

proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")
proj_k27_HSPC = loadArchRProject("./tmp/H3K27me3_coCnT_HSPC")


d_meta_k27 = as.data.table(proj_k27@cellColData)
d_meta_k27_HSPC = as.data.table(proj_k27_HSPC@cellColData)
d_meta_k27$Clusters_HSPC = NULL

d_meta_k27 = merge(d_meta_k27, d_meta_k27_HSPC[, .(cell_id, Clusters_HSPC = Clusters)], by = "cell_id", all = T)


## add back to the project
proj_k27$Clusters_HSPC = d_meta_k27$Clusters_HSPC[match(proj_k27$cell_id, d_meta_k27$cell_id)]
proj_k27$Clusters_HSPC[is.na(proj_k27$Clusters_HSPC)] = "Other_cells"

# hspc_cluster2name = c(
#     "C1" = "HSC&MPP",
#     "C2" = "GMP",
#     "C3" = "MEP",
#     "C4" = "Erythroid progenitors1",
#     "C5" = "Erythroid progenitors1",
#     "C6" = "Pre-Pro-B cells",
#     "C7" = "Myelocyte/Classical Monocytes1",
#     "C8" = "Myelocyte/Classical Monocytes1"
#     )

# table(proj_k27$ct3)

ct3 = proj_k27$Clusters_HSPC #%>% revalue(replace = hspc_cluster2name)
ct2 = proj_k27$ct2

ct3[ct3 == "Other_cells"] = ct2[ct3 == "Other_cells"]

proj_k27$ct3 = ct3

proj_k27$ct3 %>% table

proj_k27$Clusters_HSPC %>% table
proj_k27_HSPC$Clusters_HSPC %>% table
# proj_k27_HSPC$Clusters_HSPC %<>% revalue(replace = hspc_cluster2name)

saveArchRProject(proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_final3/", load = F)
saveArchRProject(proj_k27_HSPC, outputDirectory = "./tmp/H3K27me3_coCnT_HSPC/", load = F)

httpgd::hgd(port = 4322)

plotEmbedding(
	ArchRProj = proj_k27_HSPC,
	colorBy = "cellColData",
	name = "Clusters_HSPC",
	embedding = "UMAP_Harmony"
	)


pdf("./figures/206_UMAP_HSPC_clusters.pdf", width = 5, height = 4)
plotEmbedding(
	ArchRProj = proj_k27,
	colorBy = "cellColData",
	name = "Clusters_HSPC",
	embedding = "UMAP_Harmony"
	)
dev.off()

plotTSSEnrichment(
	ArchRProj = proj_k27_HSPC,
	groupBy = "Clusters_HSPC"
) + theme(legend.position = "right")
