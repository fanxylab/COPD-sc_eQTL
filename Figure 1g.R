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

read_rs_data <- function(file_path) {
  dt <- fread(file_path, select = c(1, 2, 4), 
              col.names = c("chrom", "rsid", "position"))
  dt[, chrom.position := paste(chrom, position, sep = ".")]
  env <- new.env(hash = TRUE, size = nrow(dt)*1.2)
  for (i in 1:nrow(dt)) {
    rsid <- as.character(dt$rsid[i])
    chrom.position <- dt$chrom.position[i]
    env[[rsid]] <- chrom.position
  }
  return(env)
}

get_chrom_positions <- function(rs_env, rs_ids) {
  unlist(lapply(rs_ids, function(x) get0(x, envir = rs_env)))
}

rs_env <- read_rs_data("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/genotype/COPD_SNP_MAF0.1_R0.8.bim")

msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")
setwd("/datg/xuxiaopeng/sc_eQTL/03_eQTL_enrichment")

for (Cell_type in colnames(msqe)[1:34]) {
  tmp <- msqe[, Cell_type]
  mat <- assay(tmp, "lfsrs")
  filtered_mat <- mat[mat <= 0.05, , drop = FALSE]
  tmp <- tmp[rownames(filtered_mat), ]
  
  feature_SNP_pairs <- rownames(tmp)
  feature_id <- sub("\\|.*", "", feature_SNP_pairs)
  snp_id <- sub(".*\\|", "", feature_SNP_pairs)
  chrom_position <- get_chrom_positions(rs_env, snp_id)
  
  feature_snp_pos <- data.frame("Locus" = feature_id, 
                                "Variant" = paste0("chr", chrom_position),
                                "EffectSize" = betas(tmp)[, Cell_type],
                                "StandardError" = errors(tmp)[, Cell_type],
                                stringsAsFactors = FALSE
  )
  
  write.table(feature_snp_pos, paste0(Cell_type, "_", "eQTL_statistics.txt"), 
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
  )
}

library(IRanges)
data <- read.table("/datg/xuxiaopeng/sc_eQTL/03_eQTL_enrichment/E096_15_coreMarks_hg38lift_mnemonics.bed", 
                   header = FALSE, stringsAsFactors = FALSE)
colnames(data) <- c("chr", "start", "end", "annotation")

annotation_list <- vector("list", 15)
names(annotation_list) <- paste0("class_", 1:15)

for (i in 1:nrow(data)) {
  chr <- data$chr[i]
  start <- data$start[i]
  end <- data$end[i]
  annotation <- data$annotation[i]
  
  class_id <- as.numeric(strsplit(annotation, "_")[[1]][1])
  
  if (is.null(annotation_list[[class_id]][[chr]])) {
    annotation_list[[class_id]][[chr]] <- IRanges()
  }
  
  annotation_list[[class_id]][[chr]] <- append(annotation_list[[class_id]][[chr]], IRanges(start, end))
}

query_snp <- function(chr, pos) {
  result <- logical(15)
  
  for (class_id in 1:15) {
    ranges <- annotation_list[[class_id]][[chr]]
    if (!is.null(ranges) && any(pos >= start(ranges) & pos <= end(ranges))) {
      result[class_id] <- TRUE
    }
  }
  
  return(result)
}

# snp_result <- query_snp("chr11", 2427484)
# print(snp_result)


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
eSNP_list <- list()
no_eSNP_list <- list()
# names(eSNP_list) <- paste0(c("Immune", "Epi", "Endo", "Mes"), "_", "eSNP")
data <- read.table("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/genotype/COPD_SNP_MAF0.1_R0.8.bim", 
                   header = FALSE, sep = "\t")
ALL_SNPs <- data[, 2]

for (i in 1:4) {
  mat <- assay(msqe_list[[i]], "lfsrs")
  filtered_mat <- mat[apply(mat, 1, function(row) {any(row <= 0.05)}), ]
  tmp <- msqe_list[[i]][rownames(filtered_mat), ]
  feature_SNP_pairs <- rownames(tmp)
  
  esnp_id <- unique(sub(".*\\|", "", feature_SNP_pairs))
  no_esnp_id <- setdiff(ALL_SNPs, esnp_id)
  
  eSNP_chrom_position <- get_chrom_positions(rs_env, esnp_id)
  no_eSNP_chrom_position <- get_chrom_positions(rs_env, no_esnp_id)
  eSNP_list[[i]] <- eSNP_chrom_position
  no_eSNP_list[[i]] <- no_eSNP_chrom_position
  
}

library(parallel)
process_snp <- function(snp) {
  chr <- paste0("chr", unlist(strsplit(snp, "[.]"))[1])
  pos <- as.numeric(unlist(strsplit(snp, "[.]"))[2])
  return(query_snp(chr, pos))
}

results <- list(data.frame(), data.frame(), data.frame(), data.frame())

for (i in 1:4) {
  
  esnp_results <- mclapply(eSNP_list[[i]], process_snp, mc.cores = 40)
  esnp_results <- as.data.frame(do.call(rbind, esnp_results))
  esnp_results$SNP_type <- "eSNP"
  
  no_esnp_results <- mclapply(no_eSNP_list[[i]], process_snp, mc.cores = 40)
  no_esnp_results <- as.data.frame(do.call(rbind, no_esnp_results))
  no_esnp_results$SNP_type <- "no_eSNP"
  
  merged_snp_results <- rbind(esnp_results, no_esnp_results)
  
  colnames(merged_snp_results) <- c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", 
                                    "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies", "SNP_type")
  
  results[[i]] <- merged_snp_results
}


lineage_types <- c("Immune", "Epithelial", "Endothelial", "Mesenchymal")
count_summary_list <- lapply(1:length(results), function(i) {
  count_summary <- results[[i]] %>%
    pivot_longer(cols = -SNP_type, names_to = "Feature", values_to = "Value") %>%
    group_by(Feature, SNP_type, Value) %>%
    summarise(Count = n(), .groups = "drop") %>%
    unite("Category", SNP_type, Value, sep = "_") %>%
    pivot_wider(names_from = Category, values_from = Count, values_fill = list(Count = 0)) %>%
    mutate(Lineage_type = lineage_types[i])
  return(count_summary)
})

final_count_summary <- bind_rows(count_summary_list)


library("epitools")
calculate_eSNP_TRUE_vs_no_eSNP_log_or_and_pvalue <- function(eSNP_FALSE, eSNP_TRUE, no_eSNP_FALSE, no_eSNP_TRUE) {

  contingency_table <- matrix(
    c(eSNP_TRUE, no_eSNP_TRUE, eSNP_FALSE, no_eSNP_FALSE),
    nrow = 2,
    byrow = TRUE
  )
  
  fisher_result <- epitools::oddsratio.fisher(contingency_table)
  
  log_or <- log(fisher_result$measure[2, "estimate"])
  p_value <- fisher_result$p.value[2, "fisher.exact"]
  
  return(c(log_or = log_or, p_value = p_value))
}

final_count_summary <- final_count_summary %>%
  rowwise() %>%
  mutate(
    log_or = calculate_eSNP_TRUE_vs_no_eSNP_log_or_and_pvalue(eSNP_FALSE, eSNP_TRUE, no_eSNP_FALSE, no_eSNP_TRUE)["log_or"],
    p_value = calculate_eSNP_TRUE_vs_no_eSNP_log_or_and_pvalue(eSNP_FALSE, eSNP_TRUE, no_eSNP_FALSE, no_eSNP_TRUE)["p_value"]
  )

final_count_summary$Feature <- factor(final_count_summary$Feature, 
                                      level = c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", 
                                                "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies", "SNP_type"))
final_count_summary <- final_count_summary %>%
  mutate(p_value_category = cut(
    p_value,
    breaks = c(-Inf, 1e-100, 0.05, Inf),
    labels = c("<1e-100", "<0.05", ">=0.05")
  ))

ggplot(final_count_summary, aes(x = Feature, y = Lineage_type)) +
  geom_point(aes(size = p_value_category, color = log_or)) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0
  ) +
  scale_size_manual(
    values = c(
      "<1e-100" = 2.2,
      "<0.05" = 1.5,
      ">=0.05" = 1
    ),
    name = "FDR"
  ) +
  labs(
    x = "Chromatin state",
    y = "eQTL",
    color = "Log OR",
    size = "FDR"
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
    axis.title.y = element_blank(),
    
    axis.line.y.left = element_line(linewidth = 0.1),
    axis.ticks.y.left = element_line(linewidth = 0.1),
    axis.ticks.x.bottom = element_line(linewidth = 0.1),
    axis.line.x.bottom = element_line(linewidth = 0.1)
  )

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("eSNP_enrichment.pdf", width = 6.1, height = 3.8, unit="cm")







