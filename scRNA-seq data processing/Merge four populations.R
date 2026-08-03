# author
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


setwd("/datg/xuxiaopeng/sc_eQTL/COPD")
COPD_sc_merged <- readRDS("COPD_merged_four_populations_5.rds")
# Epi 
Epi <- readRDS("COPD_merged_Epi_annotation_last.rds")
# Immune
Immune <- readRDS("COPD_merged_Immune_annotation_last.rds")
# Endo
Endo <- readRDS("COPD_merged_Endo_annotation_last.rds")
# Mes
Mes <- readRDS("COPD_merged_Mes_annotation_last.rds")

# remove last doublets detection
outlier_1 <- rownames(Epi@meta.data[Epi@meta.data$cell_subtype == "doublets", ])
outlier_2 <- rownames(Immune@meta.data[Immune@meta.data$cell_subtype == "doublets", ])
outlier_3 <- rownames(Endo@meta.data[Endo@meta.data$cell_subtype == "doublets", ])
outlier_4 <- rownames(Mes@meta.data[Mes@meta.data$cell_subtype == "doublets", ])

outlier_all <- c(outlier_1, outlier_2, outlier_3, outlier_4)
COPD_sc_merged@meta.data$Status <- rownames(COPD_sc_merged@meta.data) %in% outlier_all

# remove doublets and outlier clusters according to manual selection
COPD_sc_merged <- subset(COPD_sc_merged, subset = 
                           Status == FALSE
)

`%notin%` <- Negate(`%in%`)
Epi <- subset(Epi, subset = cell_subtype %notin% c("doublets"))
Immune <- subset(Immune, subset = cell_subtype %notin% c("doublets"))
Endo <- subset(Endo, subset = cell_subtype %notin% c("doublets"))
Mes <- subset(Mes, subset = cell_subtype %notin% c("doublets"))

# re-project
liarary(Ragas)
liarary(Neaulosa)

subclusters <- list("Epi" = "Epi", "Immune" = "Immune", "Endo" = "Endo", "Mes" = "Mes")
merged.pi <- CreatePostIntegrationObject(object = COPD_sc_merged,
                                         child.object.list = subclusters, 
                                         keep.child.object.name = FALSE,
                                         rp.main.cluster.anno = "population",
                                         rp.subcluster.colname = "cell_subtype")


COPD_sc_merged <- merged.pi$seurat.obj
COPD_sc_merged@meta.data <- 
  COPD_sc_merged@meta.data[, !names(COPD_sc_merged@meta.data) %in% 
                             c("DF.classifications", "harmony_clusters_0.01", "harmony_clusters_0.025", "harmony_clusters_0.05",
                               "harmony_clusters_0.1", "harmony_clusters_0.2", "harmony_clusters_1", "harmony_clusters_1.5", "harmony_clusters_2",
                               "Outlier", "Status", "harmony_clusters_0.5", "B"
)]

COPD_sc_merged@meta.data$Population <- COPD_sc_merged@meta.data$population
COPD_sc_merged@meta.data$Cell_subtype <- COPD_sc_merged@meta.data$subcluster_idents
COPD_sc_merged@meta.data <- 
  COPD_sc_merged@meta.data[, !names(COPD_sc_merged@meta.data) %in% 
                             c("population", "subcluster_idents"
)]

reorderCluster = c( 
  "AT1", "Transitional AT2", "AT2a", "AT2a", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated", "Differentiating ciliated", "PNEC",   # 上皮
  
  "Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", # 免疫 淋巴细胞
  
  "Classical monocytes", "Non-classical monocytes", "cDC1", "cDC2", "DC Mature", "Alveolar macrophage", "Interstitial macrophages", "Proliferating macrophages", 
  
  "B cell", "Plasma cell",  "ILC", "Mast cell", "Neutrophils",  # 免疫
  
  "Adventitial fibroalast", "Alveolar fibroalast", "Fibroalast", "Myofibroalast", "Activated myofibroalast", "SMC 1", "SMC 2", "Pericyte", "Mesothelial",  # 间质
  
  "Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic"  # 内皮
)

COPD_sc_merged@meta.data$Cell_subtype <- factor(COPD_sc_merged@meta.data$Cell_subtype, levels=reorderCluster)

saveRDS(COPD_sc_merged, file = "COPD_scRNA.rds")



