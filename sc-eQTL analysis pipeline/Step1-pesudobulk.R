# Author
# Xiaopeng Xu

# Load all packages used in this program
library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(muscat)
library(purrr)
library(scater)
library(scran)
library(scMerge)
library(stringr)
library(data.table)
library(tibble)

# Define 34 cell type for sc-eQTL mapping
 cell_subtypes <- c('gCap','NKT cell','Memory CD4 T cell','Naive CD4 T cell','Non-classical monocytes','AT2a','CD8T cell','NK cell','cDC2','AT1',
                    'Classical monocytes','Aerocyte','Adventitial fibroblast','Arterial','Alveolar macrophage','Culb 1','Venous','Interstitial macrophages',
                    'Lymphatic','AT2b','Culb 2','Mast cell','Goblet','Ciliated','Transitional AT2','XCL1+ T cell','B cell','Neutrophils','Treg T cell', 
                    'Plasma cell','Alveolar fibroblast','Proliferating T cells','Fibroblast','Basal')


gene_anno <- fread("/share/home/xuxiaopeng/databases/GRCh38/gencode.v44.basic.annotation.genes.bed", 
                   select = c(1,2,3,4,7), 
                   col.names = c("Chr", "start", "end", "pid", "feature_id"))[, Chr := sub("^chr", "", get("Chr"))]
colnames(gene_anno) <- c("#chr", "start", "end", "phenotype_id", "feature_id")

for (cell_subtype in cell_subtypes){
  
  cell_subtype_2 <- gsub(" ", "_", cell_subtype)
  
  # read cell_subtype file
  setwd("/datg/xuxiaopeng/sc_eQTL/COPD/cell_subtype")
  X <- readRDS(paste0(cell_subtype, ".rds"))
  
  #######################################################################################################
  ################################## Filter donors fewer than 10 cells ##################################
  #######################################################################################################
  # filter cell number less than 10 in one donor
  data <- as.data.frame(table(X@meta.data$sample))
  result <- data %>% 
    filter(Freq < 10)
  filter_donors <- as.character(result$Var1)
  
  `%notin%` <- Negate(`%in%`)
  X <- subset(X, sample %notin% filter_donors)
  
  #######################################################################################################
  ################################ Generate aggregrated mean of cells ###################################
  #######################################################################################################
  # transform to singlecell object
  sce <- as.SingleCellExperiment(X)
  
  # filter genes for sc-eQTL mapping
  sce <- sce[(rownames(sce) %in% genes), ]
  
  # normalization and log transformation
  sce <- computeSumFactors(sce)
  sce <- logNormCounts(sce)
  
  # aggregate cell count using mean method
  library(scuttle)
  agg <- aggregateAcrossCells(sce, 
                              sce$sample, 
                              statistics = "mean",
                              use.assay.type = "logcounts"
  )
  agg <- logcounts(agg)
  agg <- as.data.frame(agg)
  agg_2 <- cbind(feature_id = rownames(agg), agg)
  merged_df <- inner_join(gene_anno, agg_2, by = "feature_id")
  merged_df <- merged_df[!`#chr` %chin% c("X", "Y")]
  merged_df_final <- merged_df %>% select(-feature_id)
  
  # output gene expression matrix
  setwd("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/phenotype")
  write.table(merged_df_final, file= paste0(cell_subtype_2, ".bed"), quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  #######################################################################################################
  ################################ Generate covariate file  #############################################
  #######################################################################################################
  samples <- c("id", unique(X@meta.data$sample))
  covariate_df <- setNames(
    data.frame(matrix(ncol = length(samples), nrow = 0)), 
    samples
  )
  
  ###################################### 1/ncells ######################################################
  cell_num <- as.data.frame(table(X@meta.data$sample))$Freq
  cell_num <- c("nCells", cell_num)
  new_row_df <- data.frame(t(cell_num), stringsAsFactors = FALSE)
  colnames(new_row_df) <- samples
  covariate_df <- rbind(covariate_df, new_row_df)
  
  ###################################### sex,age,disease status ################################
  # read other covariates 
  sample_info <- read.table("/datg/xuxiaopeng/sc_eQTL/COPD/Sample_info.txt", sep="\t", header=TRUE)
  # disease status
  sample_info$Disease <- ifelse(startsWith(sample_info$sample_id, "C"), "Case", 
                                ifelse(startsWith(sample_info$sample_id, "H"), "Control", NA))    
  test <- as.data.frame(t(sample_info))
  colnames(test) <- test[1, ]
  test <- test[-1, ]
  test <- cbind(id=rownames(test), test)
  
  covariate_df <- bind_rows(
    covariate_df,
    test %>% 
      add_column(!!!set_names(
        lapply(setdiff(colnames(covariate_df), colnames(test)), function(x) NA), 
        setdiff(colnames(covariate_df), colnames(test))
      )) %>% 
      select(all_of(colnames(covariate_df)))
  )
  
  #########################################  expression PCs  ################################
  # PCA
  col_var <- apply(agg, 1, var)
  agg_filtered <- agg[col_var > 0, ]
  pca_result <- prcomp(t(agg_filtered), center = TRUE, scale = TRUE)
  PCs <- as.data.frame(t(pca_result$x))
  PCs <- cbind(id = rownames(PCs), PCs)
  covariate_df <- rbind(covariate_df, PCs)
  
  ########################################## genotype PCs ################################
  # read genotype PCs
  genotype_info <- read.table("/datg/xuxiaopeng/sc_eQTL/05_QTLtools/genotype/COPD_all.pca", sep=" ", header=TRUE, row.names=1)
  id <- paste("gPC", 1:123, sep="")
  genotype_info <- cbind(id=id, genotype_info)    
  
  covariate_df <- rbind(
    covariate_df,
    genotype_info %>% 
      add_column(!!!set_names(
        lapply(setdiff(colnames(covariate_df), colnames(genotype_info)), function(x) NA), 
        setdiff(colnames(covariate_df), colnames(genotype_info))
      )) %>% 
      select(all_of(colnames(covariate_df)))
  )
  
  
  # output
  setwd("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/covariates")
  write.table(covariate_df, file=paste0(cell_subtype_2, ".covariates.txt"), quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}


