## Purpose:
## Re-cluster the HSPC-related subset from multi-mark ArchR projects and
## propagate the H3K27me3-derived HSPC cluster labels to the matched H3K4me1,
## H3K4me2, and H3K4me3 projects.
##
## Inputs:
## - ArchR projects under ./tmp/H3K27me3_coCnT_final3/, ./tmp/H3K4me1_coCnT_final3/,
##   ./tmp/H3K4me2_coCnT_final3/, and ./tmp/H3K4me3_coCnT_final3/
## - Helper functions from ../lib/R/lib_umap.R
##
## Main steps:
## 1. Select HSPC-related populations from each mark-specific ArchR project.
## 2. Re-run Harmony/UMAP/clustering on the H3K27me3 subset.
## 3. Transfer HSPC cluster labels to the matched cells in the other marks.
## 4. Generate UMAP marker-gene plots, export marker tables, and update metadata.
##
## Outputs:
## - Subset ArchR projects under ./tmp/*_coCnT_HSPC/
## - Figure: ./figure/201_figure_HSPC_marker_genes_umap_all_marks_HSPC_nFrag_bias_test.pdf
## - Marker table: ./tmp/table_HSPC_marker_genes_all_markers.tsv
## - Marker table: ./tmp/table_HSPC_marker_genes_K4me3_cluster1.tsv

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
# addArchRThreads(1)
color = ArchRPalettes$stallion

source("../lib/R/lib_umap.R")

## Load project
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final3/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final3/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final3/")

proj_k27$ct2  %>% table

ct_selected = c("HSPC", "Erythroid progenitors1", "Pre-Pro-B cells", "Myelocyte/Classical Monocytes1")

proj_k27_sub = subsetArchRProject(proj_k27,
    cells = proj_k27$cellNames[proj_k27$ct2 %in% ct_selected],
    outputDirectory = "./tmp/H3K27me3_coCnT_HSPC/",
    dropCells = T,
    force = T
)

proj_h3kme1_sub = subsetArchRProject(proj_h3k4me1,
    cells = proj_h3k4me1$cellNames[proj_h3k4me1$ct2 %in% ct_selected],
    outputDirectory = "./tmp/H3K4me1_coCnT_HSPC/",
    dropCells = T,
    force = T
)

proj_h3kme2_sub = subsetArchRProject(proj_h3k4me2,
    cells = proj_h3k4me2$cellNames[proj_h3k4me2$ct2 %in% ct_selected],
    outputDirectory = "./tmp/H3K4me2_coCnT_HSPC/",
    dropCells = T,
    force = T
)

proj_h3kme3_sub = subsetArchRProject(proj_h3k4me3,
    cells = proj_h3k4me3$cellNames[proj_h3k4me3$ct2 %in% ct_selected],
    outputDirectory = "./tmp/H3K4me3_coCnT_HSPC/",
    dropCells = T,
    force = T
)

proj_k27_sub = loadArchRProject("./tmp/H3K27me3_coCnT_HSPC/")
proj_h3kme1_sub = loadArchRProject("./tmp/H3K4me1_coCnT_HSPC/")
proj_h3kme2_sub = loadArchRProject("./tmp/H3K4me2_coCnT_HSPC/")
proj_h3kme3_sub = loadArchRProject("./tmp/H3K4me3_coCnT_HSPC/")


proj_k27_sub = run_cluster(proj_k27_sub, cluster_resolusion = 0.2, cluster_corcutoff = 0.5)
proj_h3kme1_sub = run_cluster(proj_h3kme1_sub)
proj_h3kme2_sub = run_cluster(proj_h3kme2_sub)
proj_h3kme3_sub = run_cluster(proj_h3kme3_sub)

proj_k27_sub = addHarmony(
	ArchRProj = proj_k27_sub,
	reducedDims = "IterativeLSI",
	name = "Harmony",
	groupBy = "library_id",
	force = T,
	seed = 123
)

proj_k27_sub = addUMAP(
	proj_k27_sub,
	reducedDims = "Harmony",
	name = "UMAP_Harmony",
	metric = "cosine",
	force = T,
	seed = 123
)

proj_k27_sub = addClusters(
    proj_k27_sub,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    filterBias = T,
    biasClusters = 0.20,
    biasCol = "nFrags",
    biasEnrich = 5,
    biasProportion = 0.10,
    biasPval = 0.05,
    force = T,
    seed = 123
)

proj_k27_sub = addImputeWeights(proj_k27_sub)
proj_h3kme3_sub = addImputeWeights(proj_h3kme3_sub)

## add Cluster infomation of K27 to other marks

proj_h3kme1_sub$Clusters_HSPC = proj_k27_sub$Clusters[match(proj_h3kme1_sub$cell_id, proj_k27_sub$cell_id)]
proj_h3kme2_sub$Clusters_HSPC = proj_k27_sub$Clusters[match(proj_h3kme2_sub$cell_id, proj_k27_sub$cell_id)]
proj_h3kme3_sub$Clusters_HSPC = proj_k27_sub$Clusters[match(proj_h3kme3_sub$cell_id, proj_k27_sub$cell_id)]
proj_k27_sub$Clusters_HSPC = proj_k27_sub$Clusters


pdf("./figures/201_figure_HSPC_marker_genes_umap_all_marks_HSPC_nFrag_bias_test.pdf", width = 7, height = 7)


plotEmbedding(proj_h3kme3_sub, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony")
plotEmbedding(proj_h3kme3_sub, colorBy = "cellColData", name = "Clusters_HSPC", embedding = "UMAP_Harmony")
plotEmbedding(proj_h3kme3_sub, colorBy = "cellColData", name = "ct2", embedding = "UMAP_Harmony")
plotEmbedding(proj_h3kme3_sub, colorBy = "GeneScoreMatrix", name = "MECOM", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_h3kme3_sub))
plotEmbedding(proj_h3kme3_sub, colorBy = "GeneScoreMatrix", name = "CD34", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_h3kme3_sub))
plotEmbedding(proj_h3kme3_sub, colorBy = "GeneScoreMatrix", name = "PAX5", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_h3kme3_sub))
plotEmbedding(proj_h3kme3_sub, colorBy = "GeneScoreMatrix", name = "CEBPA", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_h3kme3_sub))
plotEmbedding(proj_h3kme3_sub, colorBy = "GeneScoreMatrix", name = "GATA1", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_h3kme3_sub))
plotEmbedding(proj_h3kme3_sub, colorBy = "GeneScoreMatrix", name = "RAG1", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_h3kme3_sub))

plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "Clusters_HSPC", embedding = "UMAP_Harmony")

plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "antibody_species", embedding = "UMAP_Harmony")
plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "nFrags", embedding = "UMAP_Harmony")
plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony")


plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "ct2", embedding = "UMAP_Harmony")
plotEmbedding(proj_k27_sub, colorBy = "GeneScoreMatrix", name = "MECOM", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_k27_sub))
plotEmbedding(proj_k27_sub, colorBy = "GeneScoreMatrix", name = "CD34", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_k27_sub))
plotEmbedding(proj_k27_sub, colorBy = "GeneScoreMatrix", name = "PAX5", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_k27_sub))
plotEmbedding(proj_k27_sub, colorBy = "GeneScoreMatrix", name = "CEBPA", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_k27_sub))
plotEmbedding(proj_k27_sub, colorBy = "GeneScoreMatrix", name = "GATA1", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_k27_sub))
plotEmbedding(proj_k27_sub, colorBy = "GeneScoreMatrix", name = "RAG1", embedding = "UMAP_Harmony", imputeWeights = getImputeWeights(proj_k27_sub))

plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "nFrags", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_h3kme3_sub, colorBy = "cellColData", name = "nFrags", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "antibody_species", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_k27_sub, colorBy = "cellColData", name = "library_id", embedding = "UMAP_Harmony") + theme(legend.position = "right")

dev.off()

saveArchRProject(proj_k27_sub, outputDirectory = "./tmp/H3K27me3_coCnT_HSPC/", load = T)
saveArchRProject(proj_h3kme1_sub, outputDirectory = "./tmp/H3K4me1_coCnT_HSPC/", load = T)
saveArchRProject(proj_h3kme2_sub, outputDirectory = "./tmp/H3K4me2_coCnT_HSPC/", load = T)
saveArchRProject(proj_h3kme3_sub, outputDirectory = "./tmp/H3K4me3_coCnT_HSPC/", load = T)

## Using K27 based clusters
get_marker_genes <- function(proj) {
    proj = subsetCells(proj, cellNames = proj$cellNames[!is.na(proj$Clusters)])

    markersGS <- getMarkerFeatures(
	ArchRProj = proj, 
	useMatrix = "GeneScoreMatrix", 
	groupBy = "Clusters_HSPC",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
    )
    markersGS
}

# markersGS_k27_sub = get_marker_genes(proj_k27_sub)
# markersGS_k4me1_sub = get_marker_genes(proj_h3kme1_sub)
# markersGS_k4me2_sub = get_marker_genes(proj_h3kme2_sub)
markersGS_k4me3_sub = get_marker_genes(proj_h3kme3_sub)

markersGS_list_sub = list(
    K4me3 = markersGS_k4me3_sub
)

marker_v = names(markersGS_list_sub)

lapply(1:length(marker_v), function(j) {
    x = markersGS_list_sub[[marker_v[j]]]
    mk_up = getMarkers(x, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
    mk_down = getMarkers(x, cutOff = "FDR <= 0.01 & Log2FC <= -1.25")
    cluster_name_v = names(mk_up)
    lapply(1:length(cluster_name_v), function(i) {
	cluster_name = cluster_name_v[i]
	d_up = mk_up[[cluster_name]] %>% as.data.table
	d_up$Direction = "Up"
	d_up$marker = marker_v[j]
	d_down = mk_down[[cluster_name]] %>% as.data.table
	d_down$Direction = "Down"
	d_down$marker = marker_v[j]
	d = rbind(d_up, d_down)
	d$Cluster = cluster_name
	return(d)
	}) %>% rbindlist -> marker_cluster_dt
}) %>% rbindlist -> markerList_all_dt


write_tsv(markerList_all_dt, "./tmp/table_HSPC_marker_genes_all_markers.tsv")


getMarkers(markersGS_k4me3_sub, cutOff = "FDR < 1 & abs(Log2FC) > 0")[["C1"]] %>% as.data.frame %>% 
    write_tsv("./tmp/table_HSPC_marker_genes_K4me3_cluster1.tsv")

## Update the H3K27me3 proj cell metadata
ct2 = proj_k27$ct2
cell_id_cluster3 = proj_k27_sub$cell_id[proj_k27_sub$Clusters_HSPC == "C3"]
ct2[proj_k27$cell_id %in% cell_id_cluster3] = "HSPC_C3"
proj_k27$ct3 = ct2


saveArchRProject(ArchRProj = proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_final3/", load = F)


