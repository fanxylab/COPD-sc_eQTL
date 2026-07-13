
library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scCustomize)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")

setwd("/datg/xuxiaopeng/sc_eQTL/COPD")
X = readRDS("COPD_scRNA_slim.rds")

options(repr.plot.width = 12, repr.plot.height = 10)
DimPlot_scCustom(X,
                 repel = TRUE,
                 reduction = "rp",
                 group.by = "Population",
                 pt.size = 0.5,
                 label = TRUE,
                 label.size = 2) +
  scale_colour_discrete(c("#8E24AA", "#FFA726", "#29B6F6", "#1B5E20")) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  mytheme + 
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title.x=element_text(family="sans",size=7,margin=margin(1.5,0,0,0)),
        axis.title.y=element_text(family="sans",size=7,margin=margin(0,1.5,0,0)),
        axis.text.x=element_text(family="sans",size=6,margin=margin(1,0,0,0)),
        axis.text.y=element_text(family="sans",size=6,margin=margin(0,1,0,0))
  )

reorderCluster = c(
  "AT1", "Transitional AT2", "AT2a", "AT2b", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated", "Differentiating ciliated", "PNEC", # Epithelial
  "Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", 
  "Classical monocytes", "Non-classical monocytes", "cDC1", "cDC2", "DC Mature", "Alveolar macrophage", "Interstitial macrophages", "Proliferating macrophages",
  "B cell", "Plasma cell", "ILC", "Mast cell", "Neutrophils", # Immune
  "Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast", "Myofibroblast", "Activated myofibroblast", "SMC 1", "SMC 2", "Pericyte", "Mesothelial", # MesenchyMAL
  "Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic" # Endothelial
)
X$Cell_subtype <- factor(X$Cell_subtype, levels=reorderCluster)

options(repr.plot.width = 12, repr.plot.height = 10)

DimPlot_scCustom(X,
                 reduction = "rp",
                 group.by = "Cell_subtype",
                 pt.size = 0.5,
                 #colors_use = DiscretePalette_scCustomize(num_colors = 50, palette = "varibow"),
                 label = FALSE,
                 label.size = 6,
                 raster = FALSE,
                 alpha = 0.6
) + scale_colour_discrete(name = "Cell label", 
                          label = sprintf("%d_%s", seq_along(reorderCluster), reorderCluster), 
                          type = c(
                            c(
                              "#FFB74D",  # AT1
                              "#FF7043",  # Transitional AT2
                              "#FFD54F",  # AT2a
                              "#FFC107",  # AT2b
                              "#FFA726",  # Culb 1
                              "#FF5722",  # Culb 2
                              "#FFF176",  # Goblet
                              "#FFCC80",  # Basal
                              "#FFE082",  # Ciliated
                              "#FFD740",  # Differentiating ciliated
                              "#FFEB3B"   # PNEC
                            ),
                            c(
                              "#64B5F6",  # Treg T cell
                              "#42A5F5",  # Memory CD4 T cell
                              "#2196F3",  # Naive CD4 T cell
                              "#1E88E5",  # CD8T cell
                              "#1976D2",  # XCL1+ T cell
                              "#90CAF9",  # NKT cell
                              "#81D4FA",  # NK cell
                              "#4FC3F7",  # Proliferating T cells
                              "#29B6F6",  # Classical monocytes
                              "#26C6DA",  # Non-classical monocytes
                              "#00ACC1",  # cDC1
                              "#00BCD4",  # cDC2
                              "#0097A7",  # DC Mature
                              "#80DEEA",  # Alveolar macrophage
                              "#4DD0E1",  # Interstitial macrophages
                              "#26A69A",  # Proliferating macrophages
                              "#00ACC1",  # B cell
                              "#00838F",  # Plasma cell
                              "#006064",  # ILC
                              "#84FFFF",  # Mast cell
                              "#18FFFF"   # Neutrophils
                            ),
                            c(
                              "#A5D6A7",  # Adventitial fibroblast
                              "#81C784",  # Alveolar fibroblast
                              "#66BB6A",  # Fibroblast
                              "#4CAF50",  # Myofibroblast
                              "#43A047",  # Activated myofibroblast
                              "#388E3C",  # SMC 1
                              "#2E7D32",  # SMC 2
                              "#1B5E20",  # Pericyte
                              "#76FF03"   # Mesothelial
                            ),
                            c(
                              "#CE93D8",  # Aerocyte
                              "#BA68C8",  # gCap
                              "#AB47BC",  # Venous
                              "#9C27B0",  # Arterial
                              "#8E24AA"   # Lymphatic
                            )
                          )
) +
  theme(legend.position = "none")


setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("UMAP_all.png", width = 12, height = 10, unit="in")
