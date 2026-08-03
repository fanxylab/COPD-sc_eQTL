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
  library(qvalue)
  library(RSQLite)
  library(parallel)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

setwd("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication")

read_GTEX_lung_slope <- function(file_path) {
  
  con <- dbConnect(SQLite(), dbname = "GTEX_lung_slope.db")
  
  dt <- fread(file_path, select = c(1, 2, 8), 
              col.names = c("gene_id", "variant_id", "slope"),
              nThread = 40)
  
  dt[, variant_id := sub("_b38$", "", variant_id)]
  dt[, feature_variant := paste(gene_id, variant_id, sep = "_")]
  
  dbWriteTable(con, "GTEX_data", dt, overwrite = TRUE)
  
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_feature_variant ON GTEX_data(feature_variant)")
  
  return(con)
}

database_lung <- read_GTEX_lung_slope("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Lung.allpairs.txt")

get_slope_batch <- function(con, feature_variants, batch_size = 10000) {
  batches <- split(feature_variants, ceiling(seq_along(feature_variants) / batch_size))
  
  results <- mclapply(batches, function(batch) {
    query <- paste0("SELECT feature_variant, slope FROM GTEX_data WHERE feature_variant IN ('", 
                    paste(batch, collapse = "','"), "')")
    dbGetQuery(con, query)
  }, mc.cores = 40)
  
  results <- do.call(rbind, results)
  return(results)
}



database_lung <- connect_existing_db("lung_slope")

get_slope_batch <- function(con, feature_variants, batch_size = 10000) {
  batches <- split(feature_variants, ceiling(seq_along(feature_variants) / batch_size))
  
  results <- mclapply(batches, function(batch) {
    query <- paste0("SELECT feature_variant, slope FROM GTEX_data WHERE feature_variant IN ('", 
                    paste(batch, collapse = "','"), "')")
    dbGetQuery(con, query)
  }, mc.cores = 40)
  
  results <- do.call(rbind, results)
  return(results)
}



msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")
setwd("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication")

all_results <- list()

cell_types <- colnames(msqe)[1:34]

for (Cell_type in cell_types) {
  tryCatch({
    tmp <- msqe[, Cell_type]
    mat <- assay(tmp, "lfsrs")
    mat <- as.data.frame(mat)
    
    split_names <- strsplit(rownames(mat), "\\|")
    genes <- sapply(split_names, `[`, 1)
    mat$gene <- genes
    mat$rowname <- rownames(mat)
    
    result <- mat %>%
      filter(mat[[Cell_type]] <= 0.05)
    
    feature_SNP_pairs <- result$rowname
    feature_id <- sub("\\|.*", "", feature_SNP_pairs)
    snp_id <- sub(".*\\|", "", feature_SNP_pairs)
    
    snp <- get_snp_info(rs_env, snp_id)
    feature_snp <- paste0(feature_id, "_", snp)
    
    mat_beta <- assay(tmp, "betas")
    mat_beta <- as.data.frame(mat_beta)
    mat_beta <- mat_beta[rownames(mat_beta) %in% result$rowname, , drop = FALSE]
    mat_beta$feature_variant <- feature_snp
    
    slope_result <- get_slope_batch(database_lung, feature_snp)
    
    merged_result <- inner_join(slope_result, mat_beta, by = "feature_variant")
    
    colnames(merged_result) <- c("feature_variant", "bulk_eQTL", "sc_eQTL")
    
    merged_result$cell_type <- Cell_type
    
    all_results[[Cell_type]] <- merged_result
    
  }, error = function(e) {
    message("Error in ", Cell_type, ": ", e$message)
  })
}

combined_data <- do.call(rbind, all_results)



# stacked plot
combined_data <- combined_data %>%
  mutate(
    direction_category = ifelse(
      sign(bulk_eQTL) == sign(sc_eQTL), 
      "Consistent", 
      "Inconsistent"
    )
  )

summary_data <- combined_data %>%
  group_by(cell_type, direction_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(cell_type) %>%
  mutate(
    total = sum(count),
    percentage = count / total * 100
  )

ggplot(summary_data, aes(x = cell_type, y = percentage, fill = direction_category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Consistent" = "#CD5C5C", "Inconsistent" = "#A9A9A9")) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    x = "Cell Type",
    y = "Percentage (%)"
  ) +
  mytheme +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 6), 
    
    legend.title = element_blank(),
    legend.text = element_text(margin = margin(l = 1, r = 3), size = 5),
    legend.spacing.x = unit(0.05, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0),    
    legend.position = "top",              
    legend.box.spacing = margin(2),
    legend.key.size = unit(0.15, "cm"),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.15, "cm"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6), 
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_blank(),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("eQTL_bulk_vs_sc_direction.pdf", width = 7, height = 4, unit="cm")
