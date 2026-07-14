
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(collapse)
  library(Cairo)
  library(UpSetR)
  library(dplyr)
  library(scales)
  library(ggplot2)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})         


msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/interaction/mashr_applied_significant.rds")
lfsrs <- assay(msqe, "lfsrs")
df <- as.data.frame(lfsrs)
df$gene <- sub("\\|.*", "", rownames(df))
df$rsID <- sub(".*\\|", "", rownames(df))


cell_types <- names(df)[1:34]
eGene_matrix <- df %>%
  select(all_of(cell_types)) %>%
  mutate(across(everything(), ~ ifelse(. < 0.05, 1, 0))) %>%
  as.data.frame()

eGene_matrix$gene <- df$gene

eGene_sets <- eGene_matrix %>%
  group_by(gene) %>%
  summarise(across(all_of(cell_types), ~ as.numeric(any(. == 1)))) %>%
  as.data.frame()

rownames(eGene_sets) <- eGene_sets$gene
eGene_sets$gene <- NULL

cell_type_colors <- c(
  "AT1" = "#FFB74D",  # AT1
  "Transitional AT2" = "#FF7043",  # Transitional AT2
  "AT2a" = "#FFD54F",  # AT2a
  "AT2b" = "#FFC107",  # AT2b
  "Culb 1" = "#FFA726",  # Culb 1
  "Culb 2" = "#FF5722",  # Culb 2
  "Goblet" = "#FFF176",  # Goblet
  "Basal" = "#FFCC80",  # Basal
  "Ciliated" = "#FFE082",  # Ciliated
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
  "cDC1" = "#00ACC1",  # cDC1
  "cDC2" = "#00BCD4",  # cDC2
  "DC Mature" = "#0097A7",  # DC Mature
  "Alveolar macrophage" = "#80DEEA",  # Alveolar macrophage
  "Interstitial macrophages" = "#4DD0E1",  # Interstitial macrophages
  "B cell" = "#00ACC1",  # B cell
  "Plasma cell" = "#00838F",  # Plasma cell
  "ILC" = "#006064",  # ILC
  "Mast cell" = "#84FFFF",  # Mast cell
  "Neutrophils" = "#18FFFF",  # Neutrophils
  "Adventitial fibroblast" = "#A5D6A7",  # Adventitial fibroblast
  "Alveolar fibroblast" = "#81C784",  # Alveolar fibroblast
  "Fibroblast" = "#66BB6A",  # Fibroblast
  "Myofibroblast" = "#4CAF50",  # Myofibroblast
  "Activated myofibroblast" = "#43A047",  # Activated myofibroblast
  "SMC 1" = "#388E3C",  # SMC 1
  "Pericyte" = "#1B5E20",  # Pericyte
  "Aerocyte" = "#CE93D8",  # Aerocyte
  "gCap" = "#BA68C8",  # gCap
  "Venous" = "#AB47BC",  # Venous
  "Arterial" = "#9C27B0",  # Arterial
  "Lymphatic" = "#8E24AA"   # Lymphatic
)


## pie plot
shared_level <- rowSums(eGene_sets)
total_cell_types <- ncol(eGene_sets)

category_data <- data.frame(
  level = factor(c("specific", "partial", "all"),
                 levels = c("specific", "partial", "all")),
  definition = c("只在1种细胞类型中显著", 
                 paste("在2-", total_cell_types-1, "种细胞类型中显著", sep=""),
                 paste("在所有", total_cell_types, "种细胞类型中显著"))
)

category_counts <- c(
  sum(shared_level == 1),
  sum(shared_level > 1 & shared_level < total_cell_types),
  sum(shared_level == total_cell_types)
)

category_data$count <- category_counts
category_data$percentage <- category_data$count / sum(category_data$count) * 100

ggplot(category_data, aes(x = "", y = count, fill = level)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(level, "\n", count, " (", round(percentage, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            color = "white", size = 1.5, fontface = "bold") +
  scale_fill_manual(values = c("specific" = "#E41A1C", 
                               "partial" = "#377EB8", 
                               "all" = "#4DAF4A")) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("COPD_interaction_pie.pdf", width = 6.1, height = 3.8, unit="cm")


## upset plot
specific_genes <- eGene_sets[shared_level == 1, ]

cell_type_counts <- colSums(specific_genes)

cell_type_df <- data.frame(
  cell_type = names(cell_type_counts),
  count = cell_type_counts
) %>%
  arrange(desc(count)) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type))

ggplot(cell_type_df, aes(x = reorder(cell_type, -count), y = count)) +
  geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
  labs(
    y = "Number of specific COPD-interaction genes"
  ) +
  mytheme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 6), 
    
    axis.text.y = element_text(size = 5),
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_text(aes(label = count), vjust = -0.5, size = 1.5)

ggsave("COPD_interaction_one.pdf", width = 8.2, height = 5.5, unit="cm")



