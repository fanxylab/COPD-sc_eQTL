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

gene_counts <- data.frame(
  Cell_type = character(),
  eGene_count = numeric(),
  stringsAsFactors = FALSE
)

for (Cell_type in colnames(df)[1:34]) {
  filtered <- df[df[[Cell_type]] <= 0.05, ]
  eGene_count <- length(unique(filtered$gene))
  gene_counts <- rbind(gene_counts, data.frame(
    Cell_type = Cell_type,
    eGene_count = eGene_count
  ))
}

gene_counts <- gene_counts %>%
  mutate(Cell_type = gsub("_", " ", Cell_type))
setwd("/datg/xuxiaopeng/sc_eQTL/tables")
cell_type_stat <- read.table("cell_type_statistics.txt", sep = "\t", header = TRUE)
cell_type <- inner_join(gene_counts, cell_type_stat, by = "Cell_type")

reorderCluster = c(
  "AT1", "Transitional AT2", "AT2a", "AT2b", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated",  # Epithelial
  "Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", 
  "Classical monocytes", "Non-classical monocytes", "cDC1", "cDC2", "DC Mature", "Alveolar macrophage", "Interstitial macrophages",
  "B cell", "Plasma cell", "ILC", "Mast cell", "Neutrophils", # Immune
  "Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast", "Myofibroblast", "Activated myofibroblast", "SMC 1", "Pericyte", # MesenchyMAL
  "Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic" # Endothelial
)
cell_type$Cell_type <- factor(cell_type$Cell_type, levels=reorderCluster)
cell_type_colors <- c(
  c(
    "#FFB74D",  # AT1
    "#FF7043",  # Transitional AT2
    "#FFD54F",  # AT2a
    "#FFC107",  # AT2b
    "#FFA726",  # Culb 1
    "#FF5722",  # Culb 2
    "#FFF176",  # Goblet
    "#FFCC80",  # Basal
    "#FFE082"  # Ciliated
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
    "#1B5E20"  # Pericyte
  ),
  c(
    "#CE93D8",  # Aerocyte
    "#BA68C8",  # gCap
    "#AB47BC",  # Venous
    "#9C27B0",  # Arterial
    "#8E24AA"   # Lymphatic
  )
)

options(repr.plot.width = 20, repr.plot.height = 10)

ggplot(cell_type, aes(x = Donor_number, y = eGene_count)) +
  geom_smooth(method = "lm", se = TRUE, color = "lightblue", fill = "gray50", alpha = 0.3, linewidth = 0.2) +
  geom_point(aes(size = Cell_type_proportion, colour = Cell_type)) +
  scale_size_continuous(  
    range = c(1, 3),
    breaks = c(0.02, 0.07, 0.12)
  ) + 
  scale_colour_manual(values = cell_type_colors) +
  scale_y_continuous(limits = c(4200, 5600)) +
  labs(
    x = "Number of donors",
    y = "Number of eGenes"
  ) +
  mytheme + 
  guides(colour = "none") + 
  guides(size = guide_legend(
    title = "Cell type proportion",
    position = "inside",
    keywidth = 0.1, 
    keyheight = 0.1, 
    theme = theme(
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 4.5)
    )
  )
  ) + 
  annotate(
    "text", x = 90, y = 5000, label = expression("Pearson's " * italic(r) * " = 0.600"), 
    size = 1.8, hjust = 0
  ) +  
  annotate(
    "text", x = 90, y = 4500, label = expression(italic(P) * " value = 1.77E-4"), 
    size = 1.8, hjust = 0
  ) +  
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    axis.title.x=element_text(size=7),
    axis.title.y=element_text(size=7),
    legend.title = element_text(margin = margin(b = 2)),
    legend.text = element_text(margin = margin(l = 0.5)),
    legend.position.inside = c(0.20, 0.85)
  )

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("Donor_number_eGene_count.pdf", width = 6, height = 4, unit="cm")

# Pearson and p value
cor_test <- cor.test(cell_type$Donor_number, cell_type$eGene_count, method = "pearson")
r_value <- cor_test$estimate
r_squared <- r_value^2
p_value <- cor_test$p.value

cat("Pearson Correlation Coefficient (R):", r_value, "\n")
cat("R-squared (R^2):", r_squared, "\n")
cat("P-value:", p_value, "\n")








