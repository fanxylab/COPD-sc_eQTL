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


setwd("/datg/xuxiaopeng/sc_eQTL/COPD")
COPD_sc_merged <- readRDS("COPD_merged_four_populations_5.rds")

X <- subset(COPD_sc_merged, subset = population == "Endothelial")

X <- FindVariableFeatures(X, nfeatures = 2000)
X <- ScaleData(X, vars.to.regress = c("mito_Ratio", "S.Score", "G2M.Score"))

X <- RunPCA(X, npcs = 50)

X <- IntegrateLayers(X, 
                     method = HarmonyIntegration, 
                     orig = "pca", new.reduction = "harmony", 
                     k.anchor = 20, 
                     dims = 1:21
)

X <- FindNeighbors(X,  
                   reduction = "harmony", 
                   dims = 1:21)

X <- FindClusters(X, 
                  resolution = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5), 
                  algorithm = 1,
                  cluster.name = paste("Endo_harmony_clusters", c(0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5), sep = "_")
)

X <- RunUMAP(X, 
             dims = 1:21, 
             reduction = "harmony", 
             reduction.name = "umap"
)

X <- JoinLayers(X)

library(Azimuth)
library(patchwork)
X <- RunAzimuth(X, reference= "/share/home/xuxiaopeng/R/x86_64-pc-linux-gnu-library/4.3/lungref.SeuratData/azimuth")
X@meta.data <- X@meta.data %>%
  mutate(cell_subtype = ifelse(predicted.ann_finest_level %in% c("Smooth muscle", "Alveolar fibroblasts", "Pericytes", "Adventitial fibroblasts", "Multiciliated (non-nasal)", 
                                                                 "AT1", "Transitional Club-AT2", "Club (non-nasal)", "AT2"
  ), "doublets", cell_subtype))

# add one column
outliers <- ifelse(X@meta.data$harmony_clusters_1.5 %in% c(10,12,13), "outliers", "normal")

X@meta.data$Status <- outliers

# Annotation cell subtype
cell_subtype <- as.character(X@meta.data$Endo_harmony_clusters_0.2)

cell_subtype[cell_subtype %in% c(3)] <- "Lymphatic"
cell_subtype[cell_subtype %in% c(1)] <- "Aerocyte"
cell_subtype[cell_subtype %in% c(0,5)] <- "gCap"
cell_subtype[cell_subtype %in% c(4)] <- "Arterial"
cell_subtype[cell_subtype %in% c(2)] <- "Venous"

X@meta.data$cell_subtype <- cell_subtype

saveRDS(X, file = "COPD_merged_Endo_annotation_last.rds")
