library(Seurat)
library(slingshot)
library(scCustomize)
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(lme4)
library(phateR)
library(dplyr)
library(tidyr)
library(viridis)
library(data.table)
library(Matrix)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(pdist)
library(monocle3)
library(tradeSeq)
library(BiocParallel)
library(pheatmap)
library(SeuratWrappers)
library(tidydr)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
X <- readRDS("Dynamic.rds")

Idents(X) <- "harmony_clusters_0.5"
DimPlot(X, label = TRUE, reduction = "umap")

rootCelltype <- "5"
oupDR <- Embeddings(X)
oupDR <- data.table(celltype = X$harmony_clusters_0.5, oupDR)
tmp <- oupDR[, lapply(.SD, mean), by = "celltype"]         # celltype centroid
tmp <- tmp[celltype != rootCelltype]
tmp <- data.frame(tmp[, -1], row.names = tmp$celltype)
oupDR$sampleID <- colnames(X)
oupDR <- oupDR[celltype == rootCelltype]
oupDR <- data.frame(oupDR, row.names = oupDR$sampleID)
oupDR <- oupDR[, colnames(tmp)]
tmp <- as.matrix(pdist(oupDR, tmp))
rownames(tmp) <- rownames(oupDR)
iTip <- grep(names(which.max(rowSums(tmp))), colnames(X))     # tip cell

# Apply Monocle3 on UMAP
cds <- as.cell_data_set(X)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = colnames(cds)[iTip])
X$monocle3PT <- pseudotime(cds)


plot_cells(
  cds, color_cells_by = "pseudotime", cell_size = 0.2, 
  trajectory_graph_color = "grey28",
  trajectory_graph_segment_size = 0.2,
  cell_stroke = 0.05,
  label_branch_points = F, label_roots = F, label_leaves = F) +
  scale_color_viridis(option = "D") + 
  guides(
    color = guide_colorbar(
      ticks.colour = "white",
      ticks.linewidth = 0.1,
      direction = "horizontal",
      theme = theme(
        legend.key.width  = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.ticks.length = unit(0.05, "cm")
      )
    )
  ) +
  mytheme +
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 6), 
    
    legend.title = element_blank(),
    legend.text = element_text(margin = margin(t = 2), size = 5),
    legend.spacing = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.position = c(0.2, 0.90),
    legend.box.spacing = margin(5),
    
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )


setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
ggsave("pseudotime.png", width = 6.5, height = 6, unit = "cm")
