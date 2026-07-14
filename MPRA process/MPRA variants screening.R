library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(data.table)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")

# colocation results
target_dir <- "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/"
file_list <- list.files(path = target_dir, pattern = "\\.txt$", full.names = TRUE)
data_list <- lapply(file_list, function(x){
  read.table(x, header = TRUE, sep = "\t")
})
results <- do.call("rbind", data_list)

# 
cell_types <- c("AT1", "AT2a", "AT2b", "Culb_1", "Culb_2", 
                "Transitional_AT2", "Goblet", "Ciliated", "Basal",
                "Adventitial_fibroblast", "Alveolar_fibroblast", "Fibroblast")

gwas_types <- c("FEV1", "FVC", "FEV1-FVC", "PEF", "COPD", "COPD-CG")

# 
filter_results <- results %>% 
  filter(PP.H4.abf >= 0.75) %>%
  filter(cell_type %in% cell_types) %>%
  filter(gwas_id %in% gwas_types)

# 
cleaned_df <- filter_results %>%
  distinct(cell_type, phenotype_id, .keep_all = TRUE)


setwd("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/MPRA")

for(i in 1:nrow(cleaned_df)){
  
  cell_type <- cleaned_df[i, 1]
  gene <- cleaned_df[i, 4]
  
  nominal <- fread(paste("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result/",
                         cell_type,".nominal.all.txt",sep=""), sep = "\t")
  snp_1mb <- nominal %>% 
    filter(phenotype_id == gene) %>%
    select("phenotype_id","variant_id","pval_nominal","slope", 
           "chrom", "pos", "slope_se","maf","n")
  
  write.table(snp_1mb, paste(cell_type, gene, "eQTL.txt", sep="_"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)    
  
}


# colocation results
target_dir <- "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/MPRA"
file_list <- list.files(path = target_dir, pattern = "\\.txt$", full.names = TRUE)
data_list <- lapply(file_list, function(x){
  read.table(x, header = TRUE, sep = "\t")
})
eqtl_results <- do.call("rbind", data_list)


filtered_eqtl_results <- eqtl_results %>% 
  filter(pval_nominal < 0.0001) %>% 
  distinct(variant_id, .keep_all = TRUE)



eqtl_df <- filtered_eqtl_results
peak_gr <- readRDS("/datg/xuxiaopeng/sc_eQTL/ATAC/peaks.rds")

eqtl_gr <- GRanges(
  seqnames = paste0("chr", eqtl_df$chrom),
  ranges = IRanges(start = eqtl_df$pos, end = eqtl_df$pos),
  strand = "*",
  mcols = eqtl_df[, c("phenotype_id", "variant_id", "pval_nominal", "slope")]
)
names(mcols(eqtl_gr)) <- c("phenotype_id", "variant_id", "pval_nominal", "slope")

overlaps <- findOverlaps(eqtl_gr, peak_gr)

eqtl_in_peaks <- eqtl_gr[queryHits(overlaps)]


ensg_short <- substr(unique(eqtl_in_peaks$phenotype_id), 1, 15)  
gene_map <- bitr(ensg_short, 
                 fromType = "ENSEMBL", 
                 toType = c("SYMBOL", "ENTREZID"), 
                 OrgDb = org.Hs.eg.db)

eqtl_in_peaks <- eqtl_gr[queryHits(overlaps)]
mcols(eqtl_in_peaks)$peak_info <- mcols(peak_gr[subjectHits(overlaps)])$peak_called_in

results_df <- as.data.frame(eqtl_in_peaks) %>%
  dplyr::select(
    seqnames, start, 
    variant_id, phenotype_id, 
    pval_nominal, slope, peak_info
  ) %>%
  dplyr::rename(chrom = seqnames, pos = start)

all_variants_df <- as.data.frame(eqtl_gr) %>%
  mutate(in_peak = row.names(.) %in% queryHits(overlaps)) %>%
  dplyr::select(seqnames, start, variant_id, phenotype_id, pval_nominal, slope, in_peak) %>%
  dplyr::rename(chrom = seqnames, pos = start)

write.table(results_df, "/datg/xuxiaopeng/sc_eQTL/MPRA/validation_eQTLs/eQTLs_results_2.txt", 
            sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE
            )

