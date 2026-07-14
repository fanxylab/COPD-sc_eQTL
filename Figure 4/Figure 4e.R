library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(data.table)
library(ComplexHeatmap)
library(circlize)

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
min_anova_df <- read.table("dynamic_eQTL_results_2.txt", sep = "\t", header = TRUE)

min_anova_df <- min_anova_df %>% 
  filter(FDR <= 0.05)

expr_zscore <- read.table("dynamic_eQTL_sig_2.txt", sep = "\t", header = TRUE)

eQTL_pair <- read.table("eQTLs_for_dynamic_mapping.txt", sep="\t", header = TRUE)

eQTL_pair_filtered <- eQTL_pair %>%
  inner_join(
    min_anova_df, 
    by = c("gene_symbol" = "phenotye", "rsID" = "variants")
  ) %>%
  distinct(gene_symbol, rsID, .keep_all = TRUE)



setDT(eQTL_pair_filtered)
result_df <- eQTL_pair_filtered[, .(gene_id = gene, gene_symbol, rsID)]

for (q_num in 1:5) {
  file_path <- sprintf(
    "/datg/xuxiaopeng/sc_eQTL/06_dynamic/quantile/nominal_result/Q%d.nominal.all.txt", q_num
    )
  
  q_data <- fread(file_path, 
                  select = c("phenotype_id", "variant_id", "slope"))
  
  setnames(q_data, 
           old = c("phenotype_id", "variant_id", "slope"),
           new = c("gene_id", "rsID", paste0("Q", q_num)))
  
  result_df <- merge(result_df, q_data, 
                     by = c("gene_id", "rsID"), 
                     all.x = TRUE)
}

print(result_df)


effect_matrix <- as.matrix(result_df[, .(Q1, Q2, Q3, Q4, Q5)])
rownames(effect_matrix) <- result_df$gene_symbol
colnames(effect_matrix) <- paste0("Q", 1:5)

effect_zscore <- t(scale(t(effect_matrix)))

common_genes <- intersect(rownames(expr_zscore), rownames(effect_zscore))
expr_aligned <- expr_zscore[common_genes, ]
effect_aligned <- effect_zscore[common_genes, ]



up_gene   <- "RGCC"
down_gene <- "SPTLC3"
highlight_genes <- c(up_gene, down_gene)
gene_indices <- which(rownames(effect_aligned) %in% highlight_genes)


right_annotation <- rowAnnotation(
  link = anno_mark(
    at = gene_indices, 
    labels = rownames(effect_aligned)[gene_indices],
    labels_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
    link_gp = gpar(lwd = 1, col = "red")
  )
)

heatmap_colors <- colorRamp2(
  c(-2, 0, 2), 
  c("#0C8282", "white", "#A50F15")
)


ht_expr <- Heatmap(
  as.matrix(expr_aligned),  
  name = "expression_z",
  col = heatmap_colors,
  
  heatmap_legend_param = list(
    title = "Mean gene expression (Z score)",  
    title_position = "topcenter",
    direction = "horizontal"
  ),
  row_dend_gp = gpar(lwd = 0.5),
  column_dend_gp = gpar(lwd = 0.5),
  row_title = NULL,
  column_title = "Pseudotime window",  
  column_title_side = "bottom",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE
)

ht_effect <- Heatmap(
  as.matrix(effect_aligned),  
  name = "effect_z",
  col = heatmap_colors,  
  
  heatmap_legend_param = list(
    title = "Effect Size (Z score)",  
    title_position = "topcenter",
    direction = "horizontal"
  ),
  row_dend_gp = gpar(lwd = 0.5),
  column_dend_gp = gpar(lwd = 0.5),
  row_title = NULL,
  column_title = "Pseudotime window",  
  column_title_side = "bottom",
  cluster_rows = FALSE,
  show_row_dend = FALSE, 
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  right_annotation = right_annotation 
)

combined_ht <- ht_expr + ht_effect


pdf("Combined_Expression_Effect_Heatmap_FIXED.pdf", width = 5.5, height = 7.5) 
draw(
  combined_ht, 
  heatmap_legend_side = "bottom",
  merge_legends = TRUE
)
dev.off()


