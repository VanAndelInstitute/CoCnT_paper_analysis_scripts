#####################################################
### Script: run_n_gene_cpg_island_promoter.R
### Purpose: Identify genes with promoters overlapping CpG islands in hg19
### OUTPUT:
###   - Text file with gene symbols of genes with promoters overlapping CpG islands in hg19
###   * ./tmp/table_gene_cpg_island_promoter_hg19.txt

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationHub)
library(IRanges)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# promoters: change upstream/downstream as you like
prom <- promoters(genes(txdb), upstream = 2000, downstream = 200)

# CpG islands from UCSC via AnnotationHub
ah <- AnnotationHub()

cpg_q <- query(ah, c("UCSC", "hg19", "cpgIslandExt"))
cpg <- cpg_q[[1]]   # GRanges

hits <- findOverlaps(prom, cpg, ignore.strand = TRUE)
gene_ids <- unique(names(prom)[queryHits(hits)])

# map Entrez -> Symbol
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = gene_ids,
                                      keytype = "ENTREZID", column = "SYMBOL",
                                      multiVals = "first")
gene_symbols <- unique(na.omit(gene_symbols))
gene_symbols

writeLines(gene_symbols, "./tmp/table_gene_cpg_island_promoter_hg19.txt")
