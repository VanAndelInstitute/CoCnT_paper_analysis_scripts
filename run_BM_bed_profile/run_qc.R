## Script: run_qc.R
## Purpose: Perform fragment-level QC for BM CoCnT libraries by summarizing
## fragments per cell, merging with sample metadata, and generating diagnostic
## plots across libraries, antibody targets, donors, and library types.
## Inputs:
## - ../data/sample_sheet_PR001798.tsv
## - ../data/sample_sheet_PR001799.tsv
## - ../data/sample_sheet_PR001855.tsv
## - ../data/sample_sheet_PR001856.tsv
## - raw BED fragment files under ../data/Janssens_Lab_Data/BMMC_Beds/*/BED/
## Outputs:
## - ./tmp/all_fragments.rds
## - ./tmp/cell_fragments.tsv

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


#+ eval=F
bed_file1 = dir("../data/Janssens_Lab_Data/BMMC_Beds/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/", pattern = "bed.gz$", full = T)
bed_file2 = dir("../data/Janssens_Lab_Data/BMMC_Beds/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/", pattern = "bed.gz$", full = T)
bed_file3 = dir("../data/Janssens_Lab_Data/BMMC_Beds/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/", pattern = "bed.gz$", full = T)
bed_file4 = dir("../data/Janssens_Lab_Data/BMMC_Beds/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/", pattern = "bed.gz$", full = T)
bed_file = c(bed_file1, bed_file2, bed_file3, bed_file4)
sample_name =  str_replace(basename(bed_file), ".bed.gz", "")

d = lapply(seq_along(bed_file), function(i) {
    d = fread(bed_file[i])
    d$sample_name = sample_name[i]
    d
    }) %>% rbindlist()

names(d) = c("chr", "start", "end", "cell_id", "depth", "sample_name")

## check the size of d_l in memory
object.size(d) %>% print(units = "auto")

write_rds(d, "./tmp/all_fragments.rds")

d_plot = d[, .(nFrags = .N), .(cell_id, sample_name)]
write_tsv(d_plot, "./tmp/cell_fragments.tsv")

#+ eval=T
d_plot = fread("./tmp/cell_fragments.tsv")

d_plot = merge(d_plot, sample_sheet, by.x = "sample_name", by.y = "sample_name")

ggplot(d_plot[nFrags > 100, .N, by = .(library_id, antibody_target, library_type)]) + aes(x = library_id, y = N, fill = antibody_target) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_wrap(~ library_type) +
	theme_classic() +
	scale_fill_nejm() +
	scale_y_log10() +
	labs(title = "Number of Cells per Library", x = "Library ID", y = "Cell Count") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(d_plot[nFrags > 250, .N, by = .(library_id, antibody_target, library_type)]) + aes(x = library_id, y = N, fill = antibody_target) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_wrap(~ library_type) +
	theme_classic() +
	scale_fill_nejm() +
	labs(title = "Number of Cells per Library", x = "Library ID", y = "Cell Count") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))


## Plot Fragments Per Cell histogram
ggplot(d_plot[nFrags > 10, nFrags, .(cell_id, antibody_target, library_type, donor_id)]) + 
    geom_histogram(aes(x = nFrags, fill = antibody_target), position = "identity", alpha = 0.8, bins = 50) +
    facet_grid(donor_id ~library_type, scales = "free_y") +
    geom_vline(xintercept = 400, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 300, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 250, linetype = "dashed", color = "green") +
    geom_vline(xintercept = 100, linetype = "dashed", color = "black") +
    scale_x_log10() + scale_fill_nejm() +
    theme_classic() + 
    labs(title = "BM Fragments Per Cell histgram", x = "Fragments Per Cell (log10)", y = "Cell Count") 

## Plot Fragments Per Cell histogram
ggplot(d_plot[nFrags > 10 & library_type == "CoCnT", nFrags, .(cell_id, antibody_target, library_type, donor_id, antibody_species)]) + 
    geom_histogram(aes(x = nFrags, fill = antibody_target), position = "identity", alpha = 0.8, bins = 50) +
    facet_grid(donor_id ~antibody_species, scales = "free_y") +
    geom_vline(xintercept = 400, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 300, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 250, linetype = "dashed", color = "green") +
    geom_vline(xintercept = 100, linetype = "dashed", color = "black") +
    scale_x_log10() + scale_fill_nejm() +
    theme_classic() + 
    labs(title = "BM Fragments Per Cell histgram", x = "Fragments Per Cell (log10)", y = "Cell Count") 


d_plot

