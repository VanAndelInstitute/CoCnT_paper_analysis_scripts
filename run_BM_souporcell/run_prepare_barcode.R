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

sample_sheet = rbind(
    fread("../data/sample_sheet_PR001798.tsv"),
    fread("../data/sample_sheet_PR001799.tsv"),
    fread("../data/sample_sheet_PR001855.tsv"),
    fread("../data/sample_sheet_PR001856.tsv")
    )

arrow_file = dir("../run_BM_coCnT_create_arrow/", pattern = "arrow$", full = T)
sample_name =  str_replace(basename(arrow_file), ".arrow", "")

## create an ArchRProject
proj <- ArchRProject(
    ArrowFiles = arrow_file,
    outputDirectory = "./tmp/ArchRProject_qc",
    copyArrows = TRUE
    )

old_meta = proj@cellColData
rownames(old_meta)
new_meta = merge(old_meta, sample_sheet, by.x = "Sample", by.y = "sample_name")
colnames(new_meta)
rownames(new_meta) = rownames(old_meta)
rownames(old_meta)
proj@cellColData = new_meta

dim(new_meta)
new_meta$cell_id = rownames(new_meta) %>% sub(".*#", "", .)

write_tsv(data.frame(new_meta), "./tmp/cell_metadata_4.tsv")

d_plot = proj@cellColData %>% as.data.frame %>% data.table(keep.rownames = "cell_id")
d_plot$donor_id[is.na(d_plot$donor_id)] = "mix"

## Select cells for all samples
x = d_plot[nFrags > 250, sub(".*#", "", cell_id)]  %>% unique() #%>% writeLines("./tmp/cell_barcode_all_samples.tsv")
length(x)
writeLines(x, "./tmp/cell_barcode_PR001798_PR001799_PR001855_PR001856.tsv")



