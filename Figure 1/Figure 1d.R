message(getRversion())
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(collapse)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")
lfsrs <- assay(msqe, "lfsrs")

df <- as.data.frame(lfsrs)
df$gene <- sub("\\|.*", "", rownames(df))
df$rsID <- sub(".*\\|", "", rownames(df))

eQTL_counts <- data.frame(
  Cell_type = character(),
  eQTL_count_1_20 = numeric(),
  eQTL_count_21_40 = numeric(),
  eQTL_count_41_60 = numeric(),
  eQTL_count_61_80 = numeric(),
  eQTL_count_81_100 = numeric(),
  eQTL_count_100_ = numeric(),
  stringsAsFactors = FALSE
)

for (Cell_type in colnames(df)[1:34]) {
  
  filtered <- df[df[[Cell_type]] <= 0.05, ]
  eGene_eQTL_count <- filtered %>%
    group_by(gene) %>%
    summarise(eQTL_count = n_distinct(rsID)) 
  
  result <- eGene_eQTL_count %>%
    mutate(group = case_when(
      eQTL_count >= 1 & eQTL_count <= 20 ~ "1-20",
      eQTL_count >= 21 & eQTL_count <= 40 ~ "21-40",
      eQTL_count >= 41 & eQTL_count <= 60 ~ "41-60",
      eQTL_count >= 61 & eQTL_count <= 80 ~ "61-80",
      eQTL_count >= 81 & eQTL_count <= 100 ~ "81-100",
      eQTL_count > 100 ~ ">100"
    )) %>%
    group_by(group) %>%
    summarise(gene_count = n()) %>%
    ungroup()
  
  eQTL_counts <- rbind(eQTL_counts, data.frame(
    Cell_type = Cell_type,
    eQTL_count_1_20 = result$gene_count[1],
    eQTL_count_21_40 = result$gene_count[2],
    eQTL_count_41_60 = result$gene_count[3],
    eQTL_count_61_80 = result$gene_count[4],
    eQTL_count_81_100 = result$gene_count[5],
    eQTL_count_100_ = result$gene_count[6]
  ))
  
}

eQTL_counts <- eQTL_counts %>%
  mutate(Cell_type = gsub("_", " ", Cell_type)) %>%
  mutate(Cell_population = case_when(
    Cell_type %in% c("AT1", "Transitional AT2", "AT2a", "AT2b", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated") ~ "Epithelial",
    Cell_type %in% c("Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", 
                     "Classical monocytes", "Non-classical monocytes", "cDC1", "cDC2", "DC Mature", "Alveolar macrophage", "Interstitial macrophages",
                     "B cell", "Plasma cell", "ILC", "Mast cell", "Neutrophils") ~ "Immune",
    Cell_type %in% c("Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast", "Myofibroblast", "Activated myofibroblast", "SMC 1", "Pericyte") ~ "Mesenchymal",
    Cell_type %in% c("Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic") ~ "Endothelial"
  ))

setwd("/datg/xuxiaopeng/sc_eQTL/tables")
eGene_counts <- read.table("cell_type_statistics_full.txt", sep = "\t", header = TRUE)

result_final <- inner_join(eQTL_counts, eGene_counts[,c("Cell_type","eGene_count")], by = "Cell_type")

df_final <- result_final %>%
  mutate(total_eQTL = eQTL_count_1_20 + eQTL_count_21_40 + eQTL_count_41_60 +
           eQTL_count_61_80 + eQTL_count_81_100 + eQTL_count_100_) %>%
  mutate(across(starts_with("eQTL_count"), ~ . / total_eQTL, .names = "prop_{.col}"))

population_order <- df_final %>%
  group_by(Cell_population) %>%
  summarise(total_eGene_count = sum(eGene_count)) %>%
  arrange(desc(total_eGene_count)) %>%
  pull(Cell_population)

df_final <- df_final %>%
  mutate(Cell_population = factor(Cell_population, levels = population_order)) %>%
  arrange(Cell_population, desc(eGene_count)) %>%
  mutate(Cell_type = factor(Cell_type, levels = unique(Cell_type)))

df_long <- df_final %>%
  pivot_longer(cols = starts_with("prop_eQTL_count"),
               names_to = "Group",
               values_to = "Proportion") %>%
  mutate(Group = factor(Group, levels = c("prop_eQTL_count_1_20", "prop_eQTL_count_21_40",
                                          "prop_eQTL_count_41_60", "prop_eQTL_count_61_80",
                                          "prop_eQTL_count_81_100", "prop_eQTL_count_100_"),
                        labels = c("1-20", "21-40", "41-60", "61-80", "81-100", ">100")))

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
    "cDC1" = "#00ACC1",  # cDC1
    "cDC2" = "#00BCD4",  # cDC2
    "DC Mature" = "#0097A7",  # DC Mature
    "Alveolar macrophage" = "#80DEEA",  # Alveolar macrophage
    "Interstitial macrophages" = "#4DD0E1",  # Interstitial macrophages
    "B cell" = "#00ACC1",  # B cell
    "Plasma cell" = "#00838F",  # Plasma cell
    "ILC" = "#006064",  # ILC
    "Mast cell" = "#84FFFF",  # Mast cell
    "Neutrophils" = "#18FFFF"   # Neutrophils
  ),
  c(
    "Adventitial fibroblast" = "#A5D6A7",  # Adventitial fibroblast
    "Alveolar fibroblast" = "#81C784",  # Alveolar fibroblast
    "Fibroblast" = "#66BB6A",  # Fibroblast
    "Myofibroblast" = "#4CAF50",  # Myofibroblast
    "Activated myofibroblast" = "#43A047",  # Activated myofibroblast
    "SMC 1" = "#388E3C",  # SMC 1
    "Pericyte" = "#1B5E20"  # Pericyte
  ),
  c(
    "Aerocyte" = "#CE93D8",  # Aerocyte
    "gCap" = "#BA68C8",  # gCap
    "Venous" = "#AB47BC",  # Venous
    "Arterial" = "#9C27B0",  # Arterial
    "Lymphatic" = "#8E24AA"   # Lymphatic
  )
)

custom_eqtl_color <- c(
  "1-20" = "#d0f0f0",
  "21-40" = "#a0dcdc",
  "41-60" = "#70c0c0",
  "61-80" = "#40a0a0",
  "81-100" = "#208080",
  ">100" = "#005050"
)

max_proportion <- max(df_long$Proportion, na.rm = TRUE)
max_eGene_count <- max(df_final$eGene_count, na.rm = TRUE)
norm_const <- max_proportion / max_eGene_count


ggplot() +
  geom_bar(
    data = df_long,
    aes(x = Cell_type, y = Proportion, fill = Group),
    stat = "identity"
  ) +
  geom_point(
    data = df_final,
    aes(x = Cell_type, y = eGene_count * norm_const),
    color = "red2", size = 0.1
  ) +
  scale_y_continuous(
    name = "Proportion of eGenes",
    expand = c(0.005, 0.005),
    sec.axis = sec_axis(
      trans = ~ . / norm_const,
      name = "Number of eGenes",
      breaks = c(0, 2500, 5000, 7500)
    )
  ) +
  scale_fill_manual(
    values = custom_eqtl_color,
    name = "eQTL/eGene",
    labels = c("1-20", "21-40", "41-60", "61-80", "81-100", expression("" > "100"))
  ) +
  # scale_color_manual(
  #   values = cell_type_colors
  # ) +
  labs(y = "Proportion of eGenes") +
  guides(
    fill = guide_legend(
      position = "inside",
      keywidth = 0.3,
      keyheight = 0.2,
      theme = theme(
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 3)
      )
    )
  ) +
  mytheme +
  theme(
    axis.text.x = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),  # X轴标签旋转
    axis.title.x = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.title = element_text(margin = margin(b = 2)),
    legend.text = element_text(margin = margin(l = 2)),
    legend.key.spacing.y = unit(1, "pt"),
    legend.justification = c("right", "top"),
    legend.background = element_rect(
      fill = alpha("grey90", 0.8),
      color = NA
    ),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.line.y.right = element_line(color = "red2", linewidth = 0.1),
    axis.ticks.y.right = element_line(linewidth = 0.1),
    axis.text.y.right = element_text(color = "red2"),
    axis.title.y.right = element_text(color = "red2"),
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("eGene_independent.pdf", width = 11.8, height = 4.6, unit="cm")

ggplot() +
  geom_point(
    data = df_circle,
    aes(x = Cell_type, y = y_circle, color = Cell_type),
    size = 1,
    shape = 16 
  ) +
  scale_color_manual(
    values = cell_type_colors
  ) +
  theme_void() +
  theme(legend.position = "none")

ggsave("eGene_independent_dot.pdf", width = 11.8, height = 4.6, unit="cm")
