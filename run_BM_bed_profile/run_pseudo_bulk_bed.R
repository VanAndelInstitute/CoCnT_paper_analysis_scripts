## Script: run_pseudo_bulk_bed.R
## Purpose: Generate pseudo-bulk BED files from filtered fragment data grouped by
## cell cluster and antibody target / antibody combination, then convert them to
## coverage tracks for downstream genome-browser-style visualization.
## Inputs:
## - ../data/sample_sheet_PR001798.tsv
## - ../data/sample_sheet_PR001799.tsv
## - ../data/sample_sheet_PR001855.tsv
## - ../data/sample_sheet_PR001856.tsv
## - ../run_BM_archR_analysis/tmp/table_cell_cluster_annotation_final_round2.tsv
## - ./tmp/all_fragments.rds 
## Outputs:
## - ./tmp/table_good_cells_fragments.tsv
## - ./tmp/pseudo_bulk_bed/*_pseudo_bulk.bed
## - derived sorted BED / bedGraph / BigWig files in the same output directory

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)
library(knitr)

sample_sheet = rbind(
    fread("../data/sample_sheet_PR001798.tsv"),
    fread("../data/sample_sheet_PR001799.tsv"),
    fread("../data/sample_sheet_PR001855.tsv"),
    fread("../data/sample_sheet_PR001856.tsv")
    )

sample_sheet[, .N, by = .(batch)]
sample_sheet[, .N, by = .(library_type, antibody_target, co_antibody_target, antibody_id, co_antibody_id)]

d_annot = fread("../run_BM_archR_analysis/tmp/table_cell_cluster_annotation_final_round2.tsv")

#+ eval=F
d_bed = read_rds("./tmp/all_fragments.rds")
d_bed_filtered = d_bed[cell_id %in% d_annot$cell_id]
write_tsv(d_bed_filtered, "./tmp/table_good_cells_fragments.tsv", col_names = T)
rm(d_bed)
gc()

#+ eval=T
d_bed_filtered = fread("./tmp/table_good_cells_fragments.tsv")

## for each cell_cluster, each omics, get the fragments, and write to bed file
dir_out = "./tmp/pseudo_bulk_bed/"
dir.create(dir_out, showWarnings = F, recursive = T)
cluster_v = unique(d_annot$cluster_k27) %>% na.omit()

sample_sheet_sub = sample_sheet[library_type == "CoCnT"]

process_antibody_target <- function(antibody_target_x) {
    sample_name_CT <- sample_sheet_sub[antibody_target == antibody_target_x, sample_name]
    d_bed_sub <- d_bed_filtered[sample_name %in% sample_name_CT]
    for (cluster_i in cluster_v) {
        bed_sub_sub <- d_bed_sub[cell_id %in% d_annot[cluster_k27== cluster_i, cell_id] & sample_name %in% sample_name_CT]
        output_dir <- file.path(dir_out, paste0(antibody_target_x, "_cluster_", cluster_i, "_pseudo_bulk.bed"))
        write_tsv(bed_sub_sub, output_dir, col_names = F)
    }
}

## Process all antibody targets
for (antibody_target_x in c("H3K27me3", "H3K4me3", "H3K4me2", "H3K4me1")) {
    process_antibody_target(antibody_target_x)
}


sample_sheet[, antibody_recombinaton := paste0(sort(unique(c(antibody_target, co_antibody_target))), collapse = "_"), by = .(library_id)] 

sample_sheet_sub = sample_sheet[library_type %in% c("CoCnt1", "CoCnt2")]

for (antibody_recombination_x in unique(sample_sheet_sub$antibody_recombinaton)) {
	sample_name_CT <- sample_sheet_sub[antibody_recombinaton == antibody_recombination_x, sample_name]
	d_bed_sub <- d_bed_filtered[sample_name %in% sample_name_CT]
	for (cluster_i in cluster_v) {
		bed_sub_sub <- d_bed_sub[cell_id %in% d_annot[cluster_k27 == cluster_i, cell_id] & sample_name %in% sample_name_CT]
		output_dir <- file.path(dir_out, paste0(antibody_recombination_x, "_cluster_", cluster_i, "_pseudo_bulk.bed"))
		write_tsv(bed_sub_sub, output_dir, col_names = F)
	}
}


## bed to bigwig
bed_files = dir(dir_out, pattern = "bulk.bed$", full = T)

for (i in seq_along(bed_files)) {
    print(i)
    bed_file = bed_files[i]
    sorted_bed_file = str_replace(bed_file, "bed$", "sorted.bed")
    bg_file = str_replace(bed_file, "bed$", "bedGraph")
    bw_file = str_replace(bed_file, "bed$", "bw")

    cmd_sort = str_glue("conda run -n ucsc bedtools sort -i {bed_file} > {sorted_bed_file}")
    system(cmd_sort)
    cmd_bg = str_glue("conda run -n ucsc bedtools genomecov -bg -i {sorted_bed_file} -g ./../../data/ref_genome/hg38.sizes > {bg_file}")
    system(cmd_bg)
    cmd_bw = str_glue("conda run -n ucsc bedGraphToBigWig {bg_file} ./../../data/ref_genome/hg38.sizes {bw_file}")
    system(cmd_bw)
}

## Bed to bigwig and normalize the cell number
# for (i in seq_along(bed_files)) {
#     print(i)
#     bed_file = bed_files[i]
#     bed_file_cell_num = fread(bed_file, header=F)[, uniqueN(V4)]
#     scale_factor = 2e4 / bed_file_cell_num
#     print(scale_factor)
#
#     sorted_bed_file = str_replace(bed_file, "bed$", "sorted.bed")
#     bg_file = str_replace(bed_file, "bed$", ".scale.bedGraph")
#     bw_file = str_replace(bed_file, "bed$", ".scale.bw")
#
#     cmd_sort = str_glue("conda run -n ucsc bedtools sort -i {bed_file} > {sorted_bed_file}")
#     system(cmd_sort)
#     cmd_bg = str_glue("conda run -n ucsc bedtools genomecov -bg -scale {scale_factor} -i {sorted_bed_file} -g ./../../data/ref_genome/hg38.sizes > {bg_file}")
#     system(cmd_bg)
#     cmd_bw = str_glue("conda run -n ucsc bedGraphToBigWig {bg_file} ./../../data/ref_genome/hg38.sizes {bw_file}")
#     system(cmd_bw)
# }

## Bed to bigwig and normalize the coveraged
scale_factor_v = c()
for (i in seq_along(bed_files)) {
    print(i)
    if (grepl("NA", bed_files[i])) {
	next
    }
    bed_file = bed_files[i]
    bed_file_cell_num = fread(bed_file, header=F)
    scale_factor = (2e10 / bed_file_cell_num[, sum(abs(V2 - V3))])
    # scale_factor = 2e4 / bed_file_cell_num
    print(scale_factor)
    scale_factor_v = c(scale_factor_v, scale_factor)

    sorted_bed_file = str_replace(bed_file, "bed$", "sorted.bed")
    bg_file = str_replace(bed_file, "bed$", ".scale_co.bedGraph")
    bw_file = str_replace(bed_file, "bed$", ".scale_co.bw")

    # cmd_sort = str_glue("conda run -n ucsc bedtools sort -i {bed_file} > {sorted_bed_file}")
    # system(cmd_sort)
    cmd_bg = str_glue("conda run -n ucsc bedtools genomecov -bg -scale {scale_factor} -i {sorted_bed_file} -g ./../../data/ref_genome/hg38.sizes > {bg_file}")
    system(cmd_bg)
    cmd_bw = str_glue("conda run -n ucsc bedGraphToBigWig {bg_file} ./../../data/ref_genome/hg38.sizes {bw_file}")
    system(cmd_bw)
}

scale_factor_v


## Scale by cell number
sanitize <- function(s) {
    s <- gsub(" ", "_", s)   # spaces -> underscores
    s <- gsub("/", "-", s)   # / -> -
    s <- gsub(":", "-", s)
    s <- gsub(",", "", s)
    s <- gsub(";", "", s)
    s
}


d_meta = fread("../run_BM_archR_analysis/tmp/table_cellColData_single_all.tsv")
d_meta$ct3 = d_meta$ct2 %>% sanitize

d_cell_count = d_meta[, .N, by = .(mark, ct3)] %>% na.omit

scale_factor_v = c()
for (i in seq_along(bed_files)) {
    print(i)
    if (grepl("NA", bed_files[i])) {
	next
    }
    bed_f = bed_files[i]

    cell_type = basename(bed_f) %>% str_match("cluster_(.*?)_pseudo_bulk") %>% .[,2]
    marker = basename(bed_f) %>% str_match("^(H3K27me3|H3K4me3|H3K4me2|H3K4me1|H3K27me3_H3K4me3|H3K27me3_H3K4me2|H3K27me3_H3K4me1)_") %>% .[,2]
    cell_num = d_cell_count[mark == marker & ct3 == cell_type, N]

    bed_file = bed_files[i]
    # bed_file_cell_num = fread(bed_file, header=F)
    # scale_factor = (2e10 / bed_file_cell_num[, sum(abs(V2 - V3))])
    scale_factor = 10000 / cell_num
    print(scale_factor)
    scale_factor_v = c(scale_factor_v, scale_factor)

    sorted_bed_file = str_replace(bed_file, "bed$", "sorted.bed")
    bg_file = str_replace(bed_file, "bed$", ".scale_cell_co.bedGraph")
    bw_file = str_replace(bed_file, "bed$", ".scale_cell_co.bw")

    # cmd_sort = str_glue("conda run -n ucsc bedtools sort -i {bed_file} > {sorted_bed_file}")
    # system(cmd_sort)
    cmd_bg = str_glue("conda run -n ucsc bedtools genomecov -bg -scale {scale_factor} -i {sorted_bed_file} -g ./../../data/ref_genome/hg38.sizes > {bg_file}")
    system(cmd_bg)
    cmd_bw = str_glue("conda run -n ucsc bedGraphToBigWig {bg_file} ./../../data/ref_genome/hg38.sizes {bw_file}")
    system(cmd_bw)
}

scale_factor_v

