##################################################
## Script: run_n_get_signature_genes_round2.R
## Purpose: Identify marker genes for clusters defined in ArchR projects for H3K27me3, H3K4me1, H3K4me2, and
## H3K4me3 datasets using the GeneScoreMatrix and Wilcoxon test
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## * ./tmp/H3K4me2_coCnT_final2/
## * ./tmp/H3K4me3_coCnT_final2/
##
## OUTPUT:
## - Table of marker genes for each cluster and direction (up/down) for each dataset:
## * ./tmp/table_marker_genes_all_clusters_both_directions_round2.tsv
## * ./tmp/table_marker_genes_all_clusters_both_directions_round2_match.tsv
## - RDS file containing marker gene results for each dataset:
## * ./tmp/markersGS_list_round2.rds
## * ./tmp/markersGS_list_round2_match.rds
##

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

## Load project
proj_h3k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final2/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final2/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final2/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final2/")

addGeneScoreMatrix(proj_h3k27, force = T)
addGeneScoreMatrix(proj_h3k4me1, force = T)
addGeneScoreMatrix(proj_h3k4me2, force = T)
addGeneScoreMatrix(proj_h3k4me3, force = T)

getAvailableMatrices(proj_h3k4me1)

get_marker_genes <- function(proj) {

    proj = subsetCells(proj, cellNames = proj$cellNames[!is.na(proj$Clusters_plus)])

    # if (!("GeneScoreMatrix" %in% getAvailableMatrices(proj))) {
	# addGeneScoreMatrix(proj)
    # }

    markersGS <- getMarkerFeatures(
	ArchRProj = proj, 
	useMatrix = "GeneScoreMatrix", 
	groupBy = "Clusters_plus",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
    )

    markersGS
}

## Gene signature based on H3K27me3 clusters
proj_h3k4me1@cellColData %>% names

markersGS_k27 <- get_marker_genes(proj_h3k27)
markersGS_k4me1 <- get_marker_genes(proj_h3k4me1)
markersGS_k4me2 <- get_marker_genes(proj_h3k4me2)
markersGS_k4me3 <- get_marker_genes(proj_h3k4me3)

markersGS_list = list(
    K27 = markersGS_k27,
    K4me1 = markersGS_k4me1,
    K4me2 = markersGS_k4me2,
    K4me3 = markersGS_k4me3
)


write_rds(markersGS_list, "./tmp/markersGS_list_round2.rds")

## Get marker genes
markersGS_list = read_rds("./tmp/markersGS_list_round2.rds")

names(markersGS_list)

marker_v = names(markersGS_list)

lapply(1:length(marker_v), function(j) {
    x = markersGS_list[[marker_v[j]]]
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

markerList_all_dt[Cluster == "C19"]

write_tsv(markerList_all_dt, "./tmp/table_marker_genes_all_clusters_both_directions_round2.tsv")


## Gene signature based on each maker's clusters
get_marker_genes <- function(proj) {

    proj = subsetCells(proj, cellNames = proj$cellNames[!is.na(proj$Clusters)])

    # if (!("GeneScoreMatrix" %in% getAvailableMatrices(proj))) {
	# addGeneScoreMatrix(proj)
    # }

    markersGS <- getMarkerFeatures(
	ArchRProj = proj, 
	useMatrix = "GeneScoreMatrix", 
	groupBy = "Clusters",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
    )
    
    markersGS
}

## Gene signature based on H3K27me3 clusters
proj_h3k4me1@cellColData %>% names

markersGS_k27 <- get_marker_genes(proj_h3k27)
markersGS_k4me1 <- get_marker_genes(proj_h3k4me1)
markersGS_k4me2 <- get_marker_genes(proj_h3k4me2)
markersGS_k4me3 <- get_marker_genes(proj_h3k4me3)

markersGS_list = list(
	K27 = markersGS_k27,
	K4me1 = markersGS_k4me1,
	K4me2 = markersGS_k4me2,
	K4me3 = markersGS_k4me3
)

write_rds(markersGS_list, "./tmp/markersGS_list_round2_match.rds")

## Get marker genes
markersGS_list = read_rds("./tmp/markersGS_list_round2_match.rds")

names(markersGS_list)

marker_v = names(markersGS_list)

lapply(1:length(marker_v), function(j) {
    x = markersGS_list[[marker_v[j]]]
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

write_tsv(markerList_all_dt, "./tmp/table_marker_genes_all_clusters_both_directions_round2_match.tsv")

d = fread("./tmp/table_marker_genes_all_clusters_both_directions_round2_match.tsv")



