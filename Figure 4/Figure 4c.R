library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
library(Seurat)
library(scCustomize)
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(viridis)

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
X <- readRDS("Dynamic.rds")
X <- JoinLayers(X)

Idents(X) <- "cell_subtype_2"
DefaultAssay(X) <- "SCT"
X <- PrepSCTFindMarkers(object = X)
X.markers <- FindAllMarkers(X, only.pos = TRUE)

X.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

top5$cluster <- factor(top5$cluster, levels = c("AT2", "AT0", "Immature AT1", "AT1"))
df_sorted <- top5[order(top5$cluster), ]
genes <- df_sorted$gene

levels(X) <- c("AT2", "AT0", "Immature AT1", "AT1")

DotPlot_scCustom(X, features =  unique(genes), dot.scale = 2.5) +
  mytheme + 
  guides(
    color = guide_colorbar(
      ticks.colour = "white",
      ticks.linewidth = 0.1,
      theme = theme(
        legend.key.width  = unit(0.25, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.ticks.length = unit(0.05, "cm")
      )
    ),
    size = guide_legend(
      theme = theme(
        legend.key.width  = unit(0.1, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.key.spacing.y = unit(1, "pt")
      )
    )
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.title.x = element_text(size = 6), 
    
    legend.title = element_text(margin = margin(b = 2), size = 5),
    legend.text = element_text(margin = margin(l = 2), size = 5),
    legend.spacing = unit(0.1, "cm"),     
    legend.margin = margin(0, 0, 0, 0),         
    legend.position = "right",           
    legend.box.spacing = margin(5),
    
    axis.text.y = element_text(size = 5),
    axis.title.y = element_blank(),
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )
