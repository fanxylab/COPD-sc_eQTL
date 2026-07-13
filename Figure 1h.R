message(getRversion())
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(argparse)
  library(collapse)
  library(QTLExperiment)
  library(multistateQTL)
  library(Cairo)
  library(mashr)
  library(ashr)
  library(ComplexHeatmap)
  library(circlize)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")
msqe <- callSignificance(msqe, assay="lfsrs",  thresh=0.05, secondThresh=0.05)

sim_sig <- getSignificant(msqe)
sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
sim_top <- runPairwiseSharing(sim_top, assay = "betas", factor = 0.5, FUN = identity)

share <- sim_top@metadata$pairwiseSharing

cell_type_order <- c("Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", 
                     "Classical monocytes", "Non-classical monocytes", "cDC2", "Alveolar macrophage", "Interstitial macrophages",
                     "B cell", "Plasma cell", "Mast cell", "Neutrophils", 
                     "AT1", "Transitional AT2", "AT2a", "AT2b", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated",
                     "Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic",
                     "Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast"
)

cell_type_colors <- c(
  c(
    "AT1" = "#FFB74D",  # AT1
    "Transitional AT2" = "#FF7043",  # Transitional AT2
    "AT2a" = "#FFD54F",  # AT2a
    "AT2b" = "#FFC107",  # AT2b
    "Culb 1" = "#FFA726",  # Culb 1
    "Culb 2" = "#FF5722",  # Culb 2
    "Goblet" = "#FFF176",  # Goblet
    "Basal" = "#FFCC80",  # Basal
    "Ciliated" = "#FFE082"  # Ciliated
  ),
  c(
    "Treg T cell" = "#64B5F6",  # Treg T cell
    "Memory CD4 T cell" = "#42A5F5",  # Memory CD4 T cell
    "Naive CD4 T cell" = "#2196F3",  # Naive CD4 T cell
    "CD8T cell" = "#1E88E5",  # CD8T cell
    "XCL1+ T cell" = "#1976D2",  # XCL1+ T cell
    "NKT cell" = "#90CAF9",  # NKT cell
    "NK cell" = "#81D4FA",  # NK cell
    "Proliferating T cells" = "#4FC3F7",  # Proliferating T cells
    "Classical monocytes" = "#29B6F6",  # Classical monocytes
    "Non-classical monocytes" = "#26C6DA",  # Non-classical monocytes
    # "cDC1" = "#00ACC1",  # cDC1
    "cDC2" = "#00BCD4",  # cDC2
    # "DC Mature" = "#0097A7",  # DC Mature
    "Alveolar macrophage" = "#80DEEA",  # Alveolar macrophage
    "Interstitial macrophages" = "#4DD0E1",  # Interstitial macrophages
    "B cell" = "#00ACC1",  # B cell
    "Plasma cell" = "#00838F",  # Plasma cell
    # "ILC" = "#006064",  # ILC
    "Mast cell" = "#84FFFF",  # Mast cell
    "Neutrophils" = "#18FFFF"   # Neutrophils
  ),
  c(
    "Adventitial fibroblast" = "#A5D6A7",  # Adventitial fibroblast
    "Alveolar fibroblast" = "#81C784",  # Alveolar fibroblast
    "Fibroblast" = "#66BB6A"  # Fibroblast
    # "Myofibroblast" = "#4CAF50",  # Myofibroblast
    # "Activated myofibroblast" = "#43A047",  # Activated myofibroblast
    # "SMC 1" = "#388E3C",  # SMC 1
    # "Pericyte" = "#1B5E20"  # Pericyte
  ),
  c(
    "Aerocyte" = "#CE93D8",  # Aerocyte
    "gCap" = "#BA68C8",  # gCap
    "Venous" = "#AB47BC",  # Venous
    "Arterial" = "#9C27B0",  # Arterial
    "Lymphatic" = "#8E24AA"   # Lymphatic
  )
)

rownames(share) <- gsub("_", " ", rownames(share))
colnames(share) <- gsub("_", " ", colnames(share))

lineage <- as.data.frame(colnames(share))
colnames(lineage) <- c("Cell_type")
lineage <- lineage %>%
  mutate(Lineage = case_when(
    Cell_type %in% c("AT1", "Transitional AT2", "AT2a", "AT2b", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated") ~ "Epithelial",
    Cell_type %in% c("Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", 
                     "Classical monocytes", "Non-classical monocytes", "cDC2", "Alveolar macrophage", "Interstitial macrophages",
                     "B cell", "Plasma cell", "Mast cell", "Neutrophils") ~ "Immune",
    Cell_type %in% c("Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast") ~ "Mesenchymal",
    Cell_type %in% c("Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic") ~ "Endothelial"
  ))


options(repr.plot.width = 13, repr.plot.height = 10)

colAnn <-
  HeatmapAnnotation(
    df = lineage %>% select(-Cell_type),
    col = list(Lineage = c("Mesenchymal" = "#388E3D", "Endothelial" = "#9D27B1", "Immune" = "#208EEA", "Epithelial" = "#FF753F")),
    which = "column",
    annotation_width = unit(c(1, 4), 'cm'),
    gap = unit(1, 'mm'),
    annotation_label = "Lineage"
  )

ca1 <- columnAnnotation(
  cell_type = anno_points(
    x = rep(0.5, length(colnames(share))),
    border = FALSE,
    ylim = c(0,1),
    gp = gpar(
      col = cell_type_colors[colnames(share)],
      fill = cell_type_colors[colnames(share)]
    ),
    pch = 21,
    size = unit(5, "mm"),
    axis = FALSE
  ),
  show_annotation_name = FALSE
)

ca2 <- rowAnnotation(
  cell_type = anno_points(
    x = rep(0.5, length(colnames(share))),
    border = FALSE,
    ylim = c(0,1),
    gp = gpar(
      col = cell_type_colors[colnames(share)],
      fill = cell_type_colors[colnames(share)]
    ),
    pch = 21,
    size = unit(5, "mm"),
    axis = FALSE
  ),
  show_annotation_name = FALSE
)

# colors
library(ComplexHeatmap)
library(circlize)
# heatmap_colors <- colorRamp2(
#   c(0.39, 1), 
#   c("#0FA5A5", "#A50F15")
# )
heatmap_colors <- colorRamp2(
  c(0, 0.81, 1), 
  c("#A50F15", "white", "#0C8282")
)

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
pdf("eQTL_sharing_heatmap.pdf", width = 12, height = 8)

ht <- Heatmap(
  share,
  name = "mat",
  col = heatmap_colors,
  height = unit(6.8, "inch"),
  width = unit(6.8, "inch"),
  heatmap_legend_param = list(
    title = "Sharing",
    at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    direction = "vertical"
  ),
  row_dend_gp = gpar(lwd = 0.8),
  column_dend_gp = gpar(lwd = 0.8),
  
  row_split = lineage$Lineage,
  column_split = lineage$Lineage,
  row_title = NULL,
  column_title = NULL,
  row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), 
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_parent_dend_line = FALSE,
  
  top_annotation = colAnn,
  bottom_annotation = ca1,
  right_annotation = ca2,
  
  show_column_names = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.35)
)

draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", 
     annotation_legend_side = "right")

groups <- c("Mesenchymal", "Endothelial", "Immune", "Epithelial")

for (i in seq_along(groups)) {
  decorate_heatmap_body(
    "mat",
    row_slice = i,
    column_slice = i,
    {
      ro <- row_order(ht)[[i]]
      co <- column_order(ht)[[i]]
      
      grid.rect(
        x = unit(0, "npc"), 
        y = unit(1, "npc"),
        width = sum(lengths(co))/length(co) * unit(1, "npc") + unit(0.5, "mm"),
        height = sum(lengths(ro))/length(ro) * unit(1, "npc") + unit(0.5, "mm"),
        just = c("left", "top"),
        gp = gpar(
          lwd = 2, 
          col = "black",
          lty = 1
        )
      )
    }
  )
}

dev.off()
