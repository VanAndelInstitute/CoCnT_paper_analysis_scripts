## In this script the co-occupancy bed files from CoCnT1 and CoCnT2 were used to generate the Arrow files.
## Parameters:
## 	NO cell filtering was applied.

library(data.table)
library(magrittr)
library(readr)
library(stringr)
library(ArchR)

## get all the co-occupency bed files, and merge the coCnT1 and coCnT2 together
bed_dir1 = "../../data/Janssens_Lab_Data/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/"
bed_dir2 = "../../data/Janssens_Lab_Data/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/"
bed_dir3 = "../../data/Janssens_Lab_Data/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/"
bed_dir4 = "../../data/Janssens_Lab_Data/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/"

sample_sheet1 = fread("../../data/sample_sheet_PR001855.tsv")
sample_sheet2 = fread("../../data/sample_sheet_PR001856.tsv")
sample_sheet3 = fread("../../data/sample_sheet_PR001798.tsv")
sample_sheet4 = fread("../../data/sample_sheet_PR001799.tsv")

sample_sheet1$bed_path = paste0(bed_dir1, sample_sheet1$sample_name, ".bed.gz")
sample_sheet2$bed_path = paste0(bed_dir2, sample_sheet2$sample_name, ".bed.gz")
sample_sheet3$bed_path = paste0(bed_dir3, sample_sheet3$sample_name, ".bed.gz")
sample_sheet4$bed_path = paste0(bed_dir4, sample_sheet4$sample_name, ".bed.gz")

sample_sheet = rbind(sample_sheet1, sample_sheet2, sample_sheet3, sample_sheet4)

## for each library_id
library_id_v = unique(sample_sheet$library_id)

for (library_id_i in library_id_v) {
    bed_path_v = sample_sheet[library_id == library_id_i & library_type != "CoCnT", ]$bed_path
    output_file_name = sub("CoCnt\\d", "CoCnt1", bed_path_v[1]) %>% 
	sub(".bed.gz", ".bed", .) %>%
	basename
    d_bed_1 = fread(bed_path_v[1], header = F)
    d_bed_2 = fread(bed_path_v[2], header = F)
    write_tsv(rbind(d_bed_1, d_bed_2), str_glue("./tmp/{output_file_name}"), col_names = F)
}

## sort and compress the merged bed files
cmd = "for f in ./tmp/*CoCnt1_*.bed; do sort -k1,1 -k2,2n $f | bgzip > ${f}.gz; done"
system(cmd)

d = sample_sheet[library_type == "CoCnt1"]
d$bed_path = paste0("./tmp/", d$sample_name, ".bed.gz")

file.exists(d$bed_path) %>% table()
file.size(d$bed_path) / 1e6

addArchRGenome('hg38')
set.seed(1)
addArchRThreads(threads = 20)

## The arrow files will be created in the current working directory
arrow_file = createArrowFiles(
	inputFiles = d$bed_path,
	sampleNames = d$sample_name,
	minTSS = 0,
	nChunk = 1,
	minFrags = 0,
	maxFrags = Inf,
	addTileMat = F,
	addGeneScoreMat = F,
	excludeChr = c('chrM'),
	force = T
	)





