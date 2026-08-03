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
suppressMessages(library(DoubletFinder))
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)

##############################################################################################
####################################### load data ############################################
##############################################################################################

setwd("/datg/xuxiaopeng/sc_eQTL/COPD")
dir <- list.dirs("/datg/mazhuo/data/COPD/Matrix_after_SoupX/ourdata", recursive = F)

object_list <- list()

# load all sample 10x matrix data
for (filedir in dir[1:2]){
  scRNA_data <- Read10X(filedir)
  seurat_object <- CreateSeuratObject(counts = scRNA_data, 
                                      min.cells = 3, 
                                      min.features = 200)
  sample <- c(basename(filedir))
  seurat_object[["sample"]] <- sample
  
  object_list[[sample]] <- seurat_object
}
# merge all seurat objects into a single seurat object
COPD_sc_merged <- merge(object_list[[1]], y = object_list[-1], project = "sc_eQTL_COPD")

#############################################################################################
########################### marking doublets using DoubleFinder #############################
#############################################################################################

# doublet checking using DoubletFinder
doublet_rate_table <- data.frame(
  min_cells = c(500, seq(1000, 30000, 1000)),
  max_cells = c(1000, seq(2000, 31000, 1000)),
  rate = c(0.004, seq(0.008, 0.24, 0.008))
)

calculate_doublet_rate <- function(cell_count) {
  if (cell_count < 500) {
    return(0)
  }
  
  index <- which(cell_count >= doublet_rate_table$min_cells & cell_count < doublet_rate_table$max_cells)
  if (length(index) == 0) {
    return(max(doublet_rate_table$rate))
  }
  
  return(doublet_rate_table$rate[index])
}

# sample level doublets checking
object_list <- lapply(X = seq_along(object_list), FUN = function(i) {
  
  x <- object_list[[i]]
  
  # calculate mito percentage 
  x$mito_Ratio <- PercentageFeatureSet(x, pattern = "^MT-")
  x$mito_Ratio <- x@meta.data$mito_Ratio / 100
  
  # calculate hemoglobin gene percentage
  x$hb_Ratio <- PercentageFeatureSet(x, pattern = "^HB[^(P)]")
  x$hb_Ratio <- x@meta.data$hb_Ratio / 100
  
  x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & 
                nCount_RNA > 1000 & nCount_RNA < 50000 & 
                mito_Ratio < 0.15 & 
                hb_Ratio < 0.05)
  
  x <- NormalizeData(x)
  
  x <- CellCycleScoring(
    object = x,
    g2m.features = cc.genes$g2m.genes, 
    s.features = cc.genes$s.genes)
  
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)
  
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:20)
  
  ## pK Identification (no ground-truth)
  sweep.res <- paramSweep(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  #sample_name <- names(object_list)[i]
  #ggsave(paste0(sample_name, ".png"))
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  
  doublet_rate <- calculate_doublet_rate(nrow(x@meta.data))
  
  nExp_poi <- round(doublet_rate*nrow(x@meta.data))
  x <- doubletFinder(x, PCs = 1:20, pN = 0.25, pK = pK_choose, 
                     nExp = nExp_poi, 
                     reuse.pANN = FALSE, sct = FALSE)
  
  x@meta.data <- x@meta.data %>% 
    rename_with(~ "DF.classifications", starts_with("DF.classifications"))
  
  x <- subset(x)
  
})

# merge all seurat objects into a single seurat object
COPD_sc_merged <- merge(object_list[[1]], y = object_list[-1], project = "sc_eQTL_COPD")

COPD_sc_merged@meta.data <- COPD_sc_merged@meta.data %>% select(-starts_with("pANN"))

# save merged object, and release the memory
saveRDS(COPD_sc_merged, file = "COPD_sc_merged.rds")
rm(object_list)

#############################################################################################
############################# Four populations processing ###################################
#############################################################################################
`%notin%` <- Negate(`%in%`)
COPD_sc_merged <- subset(COPD_sc_merged, subset = subcluster_idents %notin% c("doublets"))

COPD_sc_merged <- NormalizeData(COPD_sc_merged)
COPD_sc_merged <- FindVariableFeatures(COPD_sc_merged, nfeatures = 2000)
COPD_sc_merged <- ScaleData(COPD_sc_merged, vars.to.regress = c("mito_Ratio", "S.Score", "G2M.Score"))

COPD_sc_merged <- RunPCA(COPD_sc_merged, npcs = 50)


COPD_sc_merged <- IntegrateLayers(COPD_sc_merged, method = HarmonyIntegration, 
                                  orig = "pca", new.reduction = "harmony", 
                                  k.anchor = 20, 
                                  dims = 1:30
)
COPD_sc_merged <- FindNeighbors(COPD_sc_merged, 
                                reduction = "harmony", 
                                dims = 1:30
)
COPD_sc_merged <- FindClusters(COPD_sc_merged, 
                               resolution = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0), 
                               algorithm = 1,
                               cluster.name = paste("harmony_clusters", c(0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0), sep = "_")
)
COPD_sc_merged <- RunUMAP(COPD_sc_merged,
                          dims = 1:30,
                          reduction = "harmony",
                          reduction.name = "umap"
)

cellpops <- as.character(COPD_sc_merged@meta.data$harmony_clusters_0.05)

cellpops[cellpops %in% c(1,5,6,12)] <- "Epithelial"
cellpops[cellpops %in% c(0,4,7,8,10,11)] <- "Immune"
cellpops[cellpops %in% c(2,9)] <- "Endothelial"
cellpops[cellpops %in% c(3)] <- "Mesenchymal"

COPD_sc_merged@meta.data$population <- cellpops

saveRDS(COPD_sc_merged, file = "COPD_merged_four_populations.rds")



























