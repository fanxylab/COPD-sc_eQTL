library(tidyr)
library(ggplot2)
library(dplyr)
library(readr)
library(viridis)
library(stringr)

source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
setwd("/datg/xuxiaopeng/sc_eQTL/08_motif")

data <- read_delim("centrimo.txt", delim = "\t", comment = "#") %>%
  head(20)

data$neg_log10_qvalue <- -log10(as.numeric(data$E_value))

df <- data %>%
  separate(motif_id, into = c("gene", "rest"), 
           sep = "\\.", extra = "merge", remove = FALSE
           )

df <- df %>%
  mutate(TF_family = case_when(
    gene %in% c('HIF1A', 'EPAS1', 'HIF3A', 'ARNT', 'ARNT2', 'NPAS4', 'AHR') ~ "bHLH-PAS",
    gene %in% c('GMEB1', 'GMEB2') ~ "GMEB",
    gene %in% c('CR3L2', 'CR3L4') ~ "CREB/ATF",
    gene %in% c('SRBP1', 'SRBP2') ~ "SREBP",
    gene %in% c('ZBTB2') ~ "ZBTB",
    gene %in% c('HES7') ~ "bHLH",
    gene %in% c('MYBB') ~ "MYB",
    gene %in% c('SP1') ~ "SP/KLF",
    gene %in% c('XBP1') ~ "bZIP",
    gene %in% c('HOMEZ') ~ "Homeobox",
    TRUE ~ "Other"
  ))



df_plot <- df %>%
  arrange(TF_family, motif_id) %>%
  mutate(motif_id = factor(motif_id, levels = unique(motif_id)))

ggplot(df_plot, aes(x = motif_id, y = TF_family)) +
  geom_point(aes(size = total_sites, color = neg_log10_qvalue)) +
scale_x_discrete(
  labels = function(x) {
    stringr::str_extract(x, "^[^.]+")
  }
) +
scale_size_continuous(
  name = "Total sites",
  range = c(1, 3.5),
  breaks = function(x) {
    unique(floor(pretty(seq(0, max(x), length.out = 5))))
  }
) +
  scale_color_gradient(
    name = "-log10(q value)",
    low = "#2166AC", 
    high = "#B2182B",
    breaks = function(x) {
      unique(floor(pretty(seq(0, max(x), length.out = 5))))
    }
  ) +
  labs(
    x = "Motif name",
    y = "TF family"
  ) +
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
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )

ggsave("/datg/xuxiaopeng/sc_eQTL/08_motif/motif_bubble_plot_2.pdf", 
       width = 9, height = 6, unit = "cm")
