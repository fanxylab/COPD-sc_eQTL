library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")

# colocation results
target_dir <- "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/"
file_list <- list.files(path = target_dir, pattern = "\\.txt$", full.names = TRUE)
data_list <- lapply(file_list, function(x){
  read.table(x, header = TRUE, sep = "\t")
})
results <- do.call("rbind", data_list)

# GWAS trait independent loci info
target_dir_gwas <- "/datg/xuxiaopeng/sc_eQTL/07_GWAS/processed_sumstats2"
file_list <- list.files(path = target_dir_gwas, pattern = "\\.txt$", full.names = TRUE)

results_gwas <- do.call(rbind, lapply(file_list, function(f) {
  pure_name <- tools::file_path_sans_ext(basename(f))
  gwas_id <- sub("_loci", "", pure_name)
  line_count <- length(readLines(f)) 
  
  data.frame(
    gwas_id = gwas_id,
    independent_loci_count = line_count
  )
}))

'%notin%' <- Negate('%in%')
filter_results <- results %>% 
  filter(PP.H4.abf >= 0.75) %>%
  filter(cell_type %notin% c("Pneumonia", "Pneumothorax", "PT")) %>%
  # filter(cell_type %notin% c("GTEX_v10_Lung_converter")) %>%
  group_by(cell_type, gwas_id, rsID) %>% 
  
  arrange(eQTL_pval, .by_group = TRUE) %>% 
  slice_head(n = 1) %>%
  ungroup()

results_gwas <- results_gwas %>%
  filter(gwas_id %notin% c("Pneumonia", "Pneumothorax", "PT"))

result_count <- filter_results %>%
  mutate(cell_type = gsub("_", " ", cell_type)) %>%
  group_by(cell_type, gwas_id) %>%
  summarise(coloc_count = n(), .groups = "drop") %>%
  arrange(cell_type, gwas_id) 

gwas_ids =c("Asthma", "Bronchiectasis", "Bronchitis", "Chronic_bronchitis", 
            "COPD", "ILD", "Lung_cancer", "PF", "Sarcoidosis", 
            "FEV1", "FVC", "FEV1-FVC", "PEF", "Asthma-CG", "COPD-CG", "IPF-CG"
)

# cell type information
cell_types = c('gCap','NKT_cell','Memory_CD4_T_cell','Naive_CD4_T_cell','Non-classical_monocytes','AT2a','CD8T_cell','NK_cell','cDC2','AT1',
               'Classical_monocytes','Aerocyte','Adventitial_fibroblast','Arterial','Alveolar_macrophage','Culb_1','Venous','Interstitial_macrophages',
               'Lymphatic','AT2b','Culb_2','Mast_cell','Goblet','Ciliated','Transitional_AT2','XCL1+_T_cell','B_cell','Neutrophils','Treg_T_cell', 
               'Plasma_cell','Alveolar_fibroblast','Proliferating_T_cells','Fibroblast','Basal',"GTEX_v10_Lung_converter")

full_grid <- expand.grid(
  cell_type = unique(cell_types),
  gwas_id = unique(gwas_ids),
  stringsAsFactors = FALSE
)

final_df <- full_grid %>%
  mutate(cell_type = gsub("_", " ", cell_type)) %>%
  left_join(result_count, by = c("cell_type", "gwas_id")) %>% 
  mutate(coloc_count = coalesce(coloc_count, 0L)) %>%
  left_join(results_gwas, by = "gwas_id") %>%
  select(cell_type, gwas_id, coloc_count, independent_loci_count)

final_df$proportion <- final_df$coloc_count/final_df$independent_loci_count

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
  ),
  c("GTEX_v10_Lung_converter" = "#727171")
)
NG_2021 = c("Asthma", "Bronchiectasis", "Bronchitis", "Chronic_bronchitis", "COPD", "ILD", "Lung_cancer", "PF", "Sarcoidosis")
CG_2022 = c("Asthma-CG", "COPD-CG", "IPF-CG")
NG_2023 = c("FEV1", "FVC", "FEV1-FVC", "PEF")


total_independent <- results_gwas %>%
  mutate(
    gwas_type = case_when(
      gwas_id %in% NG_2021 ~ "2021_NG",
      gwas_id %in% CG_2022 ~ "2022_CG",
      gwas_id %in% NG_2023 ~ "2023_NG",
      TRUE ~ NA_character_
    )
  ) %>% arrange(desc(independent_loci_count))


head(total_independent)

gwas_type_colors <- c("2021_NG" = "#1f77b4", "2022_CG" = "#ff7f0e", "2023_NG" = "#2ca02c")

total_coloc <- final_df %>%
  group_by(cell_type) %>%
  summarise(total_coloc = sum(coloc_count)) %>% 
  arrange(desc(total_coloc))



options(repr.plot.width = 13, repr.plot.height = 12)


gwas_type_colors <- c("2021_NG" = "#1f77b4", "2022_CG" = "#ff7f0e", "2023_NG" = "#2ca02c")

total_coloc <- final_df %>%
  group_by(cell_type) %>%
  summarise(total_coloc = sum(coloc_count)) %>% 
  arrange(desc(total_coloc))

ordered_cell_types <- total_coloc$cell_type
ordered_gwas_types <- total_independent$gwas_id

total_independent$gwas_id <- factor(total_independent$gwas_id, levels = ordered_gwas_types)
total_coloc$cell_type <- factor(total_coloc$cell_type, levels = ordered_cell_types)
final_df$cell_type <- factor(final_df$cell_type, levels = ordered_cell_types)
final_df$gwas_id <- factor(final_df$gwas_id, levels = rev(ordered_gwas_types))


bar_top <- ggplot(total_coloc, aes(x=cell_type, y=total_coloc, fill=cell_type)) +
  geom_col(width=0.8) +
  scale_fill_manual(values=cell_type_colors, guide = 'none') +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(y="Number of\ncolocalizations") +
  mytheme +
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    axis.ticks.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.ticks.y.left = element_line(colour="black",linewidth=0.1),
    axis.line.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.line.y.left = element_line(colour="black",linewidth=0.1),
    plot.margin=unit(c(1.2,0.5,-0.5,0.5),"mm")
  )

bar_left <- ggplot(total_independent, aes(x = as.numeric(factor(gwas_id)), y=independent_loci_count, fill=gwas_type)) +
  geom_col(width=0.7) +
  geom_text(aes(label=independent_loci_count), 
            color="black", size=1.7,
            hjust = 1.1
  ) +
  scale_fill_manual(values=gwas_type_colors, guide="none") +
  scale_y_continuous(expand = expansion(mult = c(0.28, 0.02)), 
                     sec.axis = dup_axis(), trans = "reverse") +
  scale_x_continuous(
    name = "gwas_id",
    breaks = seq_along(unique(total_independent$gwas_id)),
    labels = unique(total_independent$gwas_id),
    sec.axis = dup_axis(), expand = c(0.01, 0),
    trans = "reverse"
  ) +
  coord_flip() +
  labs(y="Independent Loci") +
  mytheme +
  theme(
    axis.text.y=element_blank(), 
    axis.title.y=element_blank(), axis.title.x.top = element_blank(),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 5),
    axis.line = element_blank(),
    axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), 
    axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank(), 
    axis.line.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.line.y.right = element_line(colour="black",linewidth=0.1),
    axis.ticks.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.ticks.y.right = element_line(colour="black",linewidth=0.1),
    plot.margin=unit(c(0,0.5,0.3,0.5),"mm")
  )

bubble_main <- ggplot(final_df, aes(x=cell_type, y=gwas_id)) +
  geom_point(
    data = subset(final_df, coloc_count == 0), 
    aes(size = proportion), 
    color = "grey70"
  ) +
  geom_point(
    data = subset(final_df, coloc_count > 0),
    aes(color = coloc_count, size = proportion)
  ) +
  scale_color_gradientn(
    # colors = c("#E0F3F8", "#4575B4", "#313695"),
    colors = c("#6DAED4", "#E63946"),
    name = "Colocalization Count",
    guide="none"
  ) +
  scale_size_continuous(range = c(0.2,2), guide="none") +
  scale_y_discrete(position = "right") +
  labs(x = "Cell Type", y = "Lung-related GWAS traits") +
  mytheme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 5),
    axis.title.x = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.line.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.line.y.right = element_line(colour="black",linewidth=0.1),
    axis.ticks.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.ticks.y.right = element_line(colour="black",linewidth=0.1),
    plot.margin=unit(c(0.0,0.5,0.3,0.5),"mm")
  )

spacer_plot <- plot_spacer() + theme_void()
top_combo <- (spacer_plot | bar_top) + 
  plot_layout(widths = c(0.275, 2))

layout <- top_combo / 
  (bar_left + bubble_main + plot_layout(widths = c(0.5, 2))) + 
  plot_layout(
    heights = c(1.25, 4),
    guides = "collect"
  )

print(layout)

ggsave("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/all_coloc.pdf", 
       width = 10, height = 10, unit="cm")





####### Figure 5b heatmap ########
copd_genes <- read.table(
  "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result/COPD_risk_genes.txt", 
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )
target_gwas <- c("FEV1", "FVC", "FEV1-FVC", "PEF", "COPD", "COPD-CG")


cell_lineage_info <- data.frame(
  cell_type = c(
    "AT1", "AT2a", "AT2b", "Transitional AT2", "Basal", "Ciliated", "Culb 1", "Culb 2", "Goblet",
    "Alveolar macrophage", "Interstitial macrophages", "Classical monocytes", "Non-classical monocytes", 
    "cDC2", "Mast cell", "Neutrophils", "B cell", "Plasma cell", "NK cell", "NKT cell", "CD8T cell", 
    "Memory CD4 T cell", "Naive CD4 T cell", "Treg T cell", "Proliferating T cells", "XCL1+ T cell",
    "Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast",
    "Aerocyte", "Arterial", "gCap", "Lymphatic", "Venous",
    "GTEX v10 Lung converter"
  ),
  Lineage = factor(
    rep(c("Epithelial", "Immune", "Stromal", "Endothelial", "Bulk lung"), times = c(9, 17, 3, 5, 1)),
    levels = c("Epithelial", "Immune", "Stromal", "Endothelial", "Bulk lung")
  )
)

ordered_cell_types <- cell_lineage_info$cell_type

heatmap_data <- filter_results %>%
  filter(gwas_id %in% target_gwas) %>%
  mutate(cell_type = gsub("_", " ", cell_type)) %>% 
  inner_join(copd_genes, by = c("phenotype_id" = "Gene_id")) %>%
  group_by(cell_type, Gene_symbol) %>%
  summarise(
    max_H4 = max(PP.H4.abf, na.rm = TRUE),
    gwas_count = n_distinct(gwas_id),
    .groups = "drop"
  )

valid_genes <- unique(heatmap_data$Gene_symbol)
valid_cells <- intersect(ordered_cell_types, unique(heatmap_data$cell_type))

unknown_cells <- setdiff(unique(heatmap_data$cell_type), ordered_cell_types)
valid_cells <- c(valid_cells, unknown_cells)

heatmap_data_complete <- heatmap_data %>%
  mutate(
    cell_type = factor(cell_type, levels = valid_cells),
    Gene_symbol = factor(Gene_symbol, levels = valid_genes) 
  ) %>%
  complete(cell_type, Gene_symbol) %>%
  filter(!is.na(cell_type) & !is.na(Gene_symbol))

p_main <- ggplot(heatmap_data_complete, aes(x = cell_type, y = Gene_symbol)) +
  geom_tile(aes(fill = max_H4), color = "white", size = 0.5) +
  scale_fill_gradient(
    low = "white", high = "red", na.value = "grey85", 
    name = "Max H4", limits = c(0.75, 1) 
  ) +
  geom_text(
    aes(label = ifelse(is.na(max_H4), "", sprintf("%.2f", max_H4))),
    size = 3, color = "black"
  ) +
  geom_text(
    aes(label = ifelse(is.na(gwas_count), "", paste0("n=", gwas_count))),
    size = 2.2, color = "black", fontface = "bold",
    nudge_x = 0.0, nudge_y = 0.35
  ) +
  mytheme +
  theme(
    axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "italic"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5) 
  ) +
  labs(y = "COPD Risk Genes")

anno_data <- cell_lineage_info %>% 
  filter(cell_type %in% valid_cells) %>%
  mutate(cell_type = factor(cell_type, levels = valid_cells))

if(length(unknown_cells) > 0) {
  anno_data <- bind_rows(anno_data, 
                         data.frame(cell_type = factor(unknown_cells, 
                                                       levels = valid_cells), 
                                    Lineage = NA
                                    )
                         )
}

p_anno <- ggplot(anno_data, aes(x = cell_type, y = 1, fill = Lineage)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(
    values = c("Epithelial" = "#E64B35", "Immune" = "#4DBBD5", 
               "Stromal" = "#00A087", "Endothelial" = "#3C5488", "Bulk lung" = "black"),
    na.value = "grey50",
    name = "Cell Lineage"
  ) +
  mytheme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5) 
  ) +
  labs(x = "Cell Type")

p_final <- p_main / p_anno + 
  plot_layout(heights = c(10, 0.5), guides = "collect")

print(p_final)

ggsave(
  "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/COPD_Genes_Coloc_Heatmap.pdf", 
  plot = p_final, width = 6.5, height = 9
  )

