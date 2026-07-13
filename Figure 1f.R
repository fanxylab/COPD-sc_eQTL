message(getRversion())
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(collapse)
  library(Cairo)
  library(ComplexHeatmap)
  library(circlize)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

Immune <- c("Treg T cell", "Memory CD4 T cell", "Naive CD4 T cell", "CD8T cell", "XCL1+ T cell", "NKT cell", "NK cell", "Proliferating T cells", "Classical monocytes", 
            "Non-classical monocytes", "cDC2", "Alveolar macrophage", "Interstitial macrophages", "B cell", "Plasma cell", "Mast cell", "Neutrophils")
Epi <- c("AT1", "Transitional AT2", "AT2a", "AT2b", "Culb 1", "Culb 2", "Goblet", "Basal", "Ciliated")
Endo <- c("Aerocyte", "gCap", "Venous", "Arterial", "Lymphatic")
Mes <- c("Adventitial fibroblast", "Alveolar fibroblast", "Fibroblast")

Immune <- gsub(" ", "_", Immune)
Epi <- gsub(" ", "_", Epi)
Endo <- gsub(" ", "_", Endo)
Mes <- gsub(" ", "_", Mes)

msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")
Immune_msqe <- msqe[, Immune]
Epi_msqe <- msqe[, Epi]
Endo_msqe <- msqe[, Endo]
Mes_msqe <- msqe[, Mes]
msqe_list <- list(Immune_msqe, Epi_msqe, Endo_msqe, Mes_msqe)

# get significant eQTL-eGENE pair in each lineage
for (i in 1:4) {
  mat <- assay(msqe_list[[i]], "lfsrs")
  filtered_mat <- mat[apply(mat, 1, function(row) {any(row <= 0.05)}), ]
  tmp <- msqe_list[[i]][rownames(filtered_mat), ]
  msqe_list[[i]] <- tmp
}

read_rs_data <- function(file_path) {
  dt <- fread(file_path, select = c(2, 4), 
              col.names = c("rsid", "position"))
  env <- new.env(hash = TRUE, size = nrow(dt)*1.2)
  for (i in 1:nrow(dt)) {
    rsid <- as.character(dt$rsid[i])
    position <- dt$position[i]
    env[[rsid]] <- position
  }
  return(env)
}

get_positions <- function(rs_env, rs_ids) {
  unlist(lapply(rs_ids, function(x) get0(x, envir = rs_env)))
}

rs_env <- read_rs_data("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/genotype/COPD_SNP_MAF0.1_R0.8.bim")

feature_snp_pos_all <- NULL
lineages <- c("Immune", "Epithelial", "Endothelial" , "Mesenchymal")
for (i in 1:4) {
  feature_SNP_pairs <- rownames(msqe_list[[i]])
  feature_id <- sub("\\|.*", "", feature_SNP_pairs)
  snp_id <- sub(".*\\|", "", feature_SNP_pairs)
  snp_position <- get_positions(rs_env, snp_id)
  feature_snp_pos <- data.frame("feature_id" = feature_id, 
                                "snp_id" = snp_id,
                                "snp_position" = snp_position,
                                "lineage_type" = lineages[i],
                                stringsAsFactors = FALSE
  )
  feature_snp_pos_all <- rbind(feature_snp_pos_all, feature_snp_pos)
}

gene_annotation_path <- "/share/home/xuxiaopeng/databases/GRCh38/gencode.v44.basic.annotation.genes.bed"
gene_table <- fread(gene_annotation_path, 
                    select = c(1,2,3,4,6), 
                    col.names = c("chrom", "start", "end", "feature_id","direction"))

result <- left_join(feature_snp_pos_all, gene_table, by = "feature_id")

result <- result %>%
  mutate(tss = case_when(
    direction == "+" ~ start,
    direction == "-" ~ end
  ))

result$distance <- result$snp_position - result$tss

library(ggplot2)
cols <-  c("Mesenchymal" = "#388E3D", "Endothelial" = "#9D27B1", "Immune" = "#208EEA", "Epithelial" = "#FF753F")
ggplot(result, aes(x = distance, color = lineage_type)) +
  geom_density(key_glyph = draw_key_path, 
               adjust = 2, linewidth = 0.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
  scale_x_continuous(
    name = "Distance from TSS (Kbp)",
    limits = c(-1000000, 1000000),
    breaks = seq(-1000000, 1000000, by = 500000),
    labels = c("-1000", "-500", "0", "500", "1000")
  ) +
  scale_y_continuous(
    name = "Density (10e-4)",
    limits = c(0, 1.998463e-06),
    breaks = c(0.0, 5e-07, 1e-06, 1.5e-06),
    labels = c("0.000", "0.005", "0.010", "0.015"),
    expand = c(2e-8, 2e-8)
  ) +
  scale_color_manual(values = cols) + 
  theme(legend.key = element_rect(fill = NA, color = NA)) +
  guides(color = guide_legend(override.aes = list(linetype = 1))) +
  mytheme + 
  guides(
    color = guide_legend(
      position = "inside",
      keywidth = 0.3,
      keyheight = 0.2,
      theme = theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 5)
      )
    )
  ) + 
  theme(
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 6), 
    
    legend.title = element_text(margin = margin(b = 2)),
    legend.text = element_text(margin = margin(l = 2)),
    legend.key.spacing.y = unit(3, "pt"),
    legend.justification = c("right", "top"),
    legend.position.inside = c(1.0, 0.98),
    
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("eSNP_density.pdf", width = 5, height = 3.8, unit="cm")

