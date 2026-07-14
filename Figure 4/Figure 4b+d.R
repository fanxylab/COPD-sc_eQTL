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

X@meta.data <- X@meta.data %>%
  mutate(
    Quantile = cut(
      monocle3PT,
      breaks = seq(0, 62.7, by = 7.8375),  
      include.lowest = TRUE,              
      right = FALSE,                      
      labels = paste0("Q", 1:8)          
    )
  )

X@meta.data <- X@meta.data %>%
  mutate(New_Quantile = case_when(
    Quantile %in% c("Q1", "Q2", "Q3") ~ "Q1",
    Quantile == "Q4" ~ "Q2",
    Quantile == "Q5" ~ "Q3",
    Quantile == "Q6" ~ "Q4",
    Quantile %in% c("Q7", "Q8") ~ "Q5"
  ))
X@meta.data$New_Quantile <- factor(X@meta.data$New_Quantile, levels=c("Q1","Q2","Q3","Q4","Q5"))

options(repr.plot.width = 10, repr.plot.height = 8)
# Idents(X) <- "harmony_clusters_1"
Idents(X) <- "New_Quantile"
DimPlot_scCustom(X, label = TRUE, reduction = "umap", label.size = 6, pt.size = 0.01,
                 colors_use = DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")[10:17]
) +
  mytheme +
  # guides(
  #   color = guide_legend(
  #       title = "Quantile",
  #       theme = theme(
  #           legend.key.width  = unit(0.5, "cm"),
  #           legend.key.height = unit(0.5, "cm"),
  #           legend.key.spacing.y = unit(0.5, "pt")
  #       )
  #   )
  # ) +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.15, "inches"), ends = 'last', type = "closed")
  ) + 
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    
    legend.title = element_text(margin = margin(b = 2), size = 15),
    legend.text = element_text(margin = margin(l = 2), size = 15),
    # # legend.spacing = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    # # legend.position = "right",
    # # legend.box.spacing = margin(5),
    
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )

options(repr.plot.width = 10, repr.plot.height = 8)
# Idents(X) <- "harmony_clusters_1"
Idents(X) <- "cell_subtype_2"
DimPlot_scCustom(X, label = TRUE, reduction = "umap", label.size = 6, pt.size = 0.01,
                 colors_use = DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")[10:17]
) +
  mytheme +
  # guides(
  #   color = guide_legend(
  #       title = "Quantile",
  #       theme = theme(
  #           legend.key.width  = unit(0.5, "cm"),
  #           legend.key.height = unit(0.5, "cm"),
  #           legend.key.spacing.y = unit(0.5, "pt")
  #       )
  #   )
  # ) +
  theme_dr(xlength = 0.2, ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.15, "inches"), ends = 'last', type = "closed")
  ) + 
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    
    legend.title = element_text(margin = margin(b = 2), size = 15),
    legend.text = element_text(margin = margin(l = 2), size = 15),
    # # legend.spacing = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),     
    # # legend.position = "right",
    # # legend.box.spacing = margin(5),
    
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )


