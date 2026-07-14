library(Signac)
library(Seurat)

library(GenomicRanges)
library(ggplot2)
library(patchwork)

setwd("/datg/xuxiaopeng/sc_eQTL/ATAC")

lung <- readRDS("lung_scATAC.rds")

# ZKSCAN1
# regions <- c(
# "chr7-100012504-100012504",
# "chr7-100016690-100016690",
# "chr7-100020229-100020229",
# "chr7-100035359-100035359",
# "chr7-100041777-100041777",
# "chr7-100089250-100089250"
# )
# highlight <- StringToGRanges(regions = regions)
# highlight$color <- "red"

# STN1
regions <- c(
  "chr10-103835672-103835672",
  "chr10-103859920-103859920",
  "chr10-103889814-103889814",
  "chr10-103890055-103890055",
  "chr10-103894406-103894406",
  "chr10-103901228-103901228",
  "chr10-103901557-103901557",
  "chr10-103923736-103923736",
  "chr10-103933246-103933246"
)

highlight <- StringToGRanges(regions = regions)
highlight$color <- "red"

CoveragePlot(
  object = lung,
  region = c("chr10-103833672-103942246"),
  region.highlight = highlight,
  extend.upstream = 30000,
  extend.downstream = 10000
)

ggsave("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/STN1_ATAC_seq.pdf", 
       width = 20, height = 20, unit = "cm")
