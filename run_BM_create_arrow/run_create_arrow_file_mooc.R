## Create Arrow files for BM uni-occupancy coCnT data, batches 1 and 2
## This script will create Arrow files for each sample, which will be used for downstream ArchR analysis.
## Parameters:
##   - Cells with < 100 fragments are filtered out
##   - TSS enrichment is not calculated at this step, as it is not used for filtering cells in this dataset. It can be calculated later if needed.
##   - The Arrow files will be created in the current working directory, with the same name as the sample name (e.g., "sample_name.arrow").
##   - maxFrags is set to Inf to include all fragments, as there is no upper limit on the number of fragments per cell in this dataset.
##   - excludeChr is set to "chrM" to exclude mitochondrial reads, as they are not informative for chromatin accessibility analysis and can introduce noise.
##   - does not add tile matrix and gene score matrix at this step to save time and disk space.
##
library(data.table)
library(magrittr)

library(ArchR)

set.seed(1)
addArchRGenome("hg38")
addArchRThreads(threads = 20)

batch_info = data.table(
  sample_sheet_path = c(
    "../data/sample_sheet_PR001798.tsv",
    "../data/sample_sheet_PR001799.tsv",
    "../data/sample_sheet_PR001855.tsv",
    "../data/sample_sheet_PR001856.tsv"
  ),
  bed_dir = c(
    "../data/Janssens_Lab_Data/BMMC_Beds/PR001798_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/",
    "../data/Janssens_Lab_Data/BMMC_Beds/PR001799_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/",
    "../data/Janssens_Lab_Data/BMMC_Beds/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/",
    "../data/Janssens_Lab_Data/BMMC_Beds/PR001856_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED/"
  )
)

sample_sheets = lapply(seq_len(nrow(batch_info)), function(i) {
  sample_sheet = fread(batch_info$sample_sheet_path[i])
  sample_sheet$bed_path = paste0(batch_info$bed_dir[i], sample_sheet$sample_name, ".bed.gz")
  sample_sheet
})

sample_sheet = rbindlist(sample_sheets, use.names = TRUE, fill = TRUE)
head(sample_sheet)

## Uncomment to check if missing sample
# sample_name_exist = dir("./") %>% grep("arrow$", ., value = TRUE) %>%
#   sub("\\.arrow$", "", .)
# sample_sheet[!(sample_name %in% sample_name_exist), ] -> sample_sheet

file.exists(sample_sheet$bed_path) %>% table()

file.size(sample_sheet$bed_path) / 1e6

arrow_file = createArrowFiles(
  inputFiles = sample_sheet$bed_path,
  sampleNames = sample_sheet$sample_name,
  minTSS = 0,
  nChunk = 1,
  minFrags = 100,
  maxFrags = Inf,
  addTileMat = FALSE,
  addGeneScoreMat = FALSE,
  excludeChr = c("chrM"),
  force = TRUE
)

dir("./", "arrow$")
