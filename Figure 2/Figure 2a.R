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

# STEP1 # bulk eQTL
read_rs_data <- function(file_path) {
  dt <- fread(file_path, select = c(1, 2, 3, 4, 5), 
              col.names = c("chrom", "pos", "ref", "alt", "rsid"))
  dt[, snp.info := paste(chrom, pos, ref, alt, sep = "_")]
  env <- new.env(hash = TRUE, size = nrow(dt)*1.2)
  for (i in 1:nrow(dt)) {
    rsid <- as.character(dt$rsid[i])
    snp.info <- dt$snp.info[i]
    env[[rsid]] <- snp.info
  }
  return(env)
}

get_snp_info <- function(rs_env, rs_ids) {
  unlist(lapply(rs_ids, function(x) get0(x, envir = rs_env)))
}

rs_env <- read_rs_data("/datg/xuxiaopeng/WGS/plink/merged_SNPs.txt")

# test
get_snp_info(rs_env, c("rs113141985", "rs2747966"))


# STEP2 # gene annotation
read_feature_data <- function(file_path) {
  dt <- fread(file_path, select = c(4, 7), 
              col.names = c("Ensembl", "Symbol"))
  env <- new.env(hash = TRUE, size = nrow(dt)*1.2)
  for (i in 1:nrow(dt)) {
    symbol <- as.character(dt$Symbol[i])
    ensembl <- dt$Ensembl[i]
    env[[symbol]] <- ensembl
  }
  return(env)
}

get_feature_info <- function(feature_env, feature_ids) {
  unlist(lapply(feature_ids, function(x) get0(x, envir = feature_env)))
}

feature_env <- read_feature_data("/share/home/xuxiaopeng/databases/GRCh38/gencode.v44.basic.annotation.genes.bed")

# test
get_feature_info(feature_env, c("WASH2P", "DDX11L1"))

# STEP3 # GTEX
setwd("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication")
read_GTEX_data <- function(file_path) {
  
  tissue <- gsub(".*_associations_([A-Za-z]+)\\.allpairs\\.txt", "\\1", file_path)
  clean_tissue <- gsub("[^A-Za-z0-9_]", "_", tissue)
  
  db_path <- file.path(getwd(), paste0("GETX_", clean_tissue, ".db"))
  
  con <- dbConnect(SQLite(), dbname = db_path)
  
  dt <- fread(file_path, select = c(1, 2, 7), 
              col.names = c("gene_id", "variant_id", "p_value"),
              nThread = 40)
  
  dt[, variant_id := sub("_b38$", "", variant_id)]
  dt[, feature_variant := paste(gene_id, variant_id, sep = "_")]
  
  dbWriteTable(con, "GTEX_data", dt, overwrite = TRUE)
  
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_feature_variant ON GTEX_data(feature_variant)")
  
  return(con)
}

con_lung <- read_GTEX_data("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Lung.allpairs.txt")
con_brain <- read_GTEX_data("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Brain_Cortex.allpairs.txt")
con_blood <- read_GTEX_data("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Whole_Blood.allpairs.txt")

get_p_value_batch <- function(con, feature_variants, batch_size = 10000) {
  
  batches <- split(feature_variants, ceiling(seq_along(feature_variants) / batch_size))
  
  results <- mclapply(batches, function(batch) {
    query <- paste0("SELECT p_value FROM GTEX_data WHERE feature_variant IN ('", 
                    paste(batch, collapse = "','"), "')")
    dbGetQuery(con, query)
  }, mc.cores = 40)
  
  results <- do.call(rbind, results)
  return(results)
}


setwd("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication")

connect_existing_db <- function(tissue_name) {
  db_path <- file.path(getwd(), paste0("GTEX_", tissue_name, ".db"))
  if (!file.exists(db_path)) {
    stop("数据库文件不存在: ", db_path)
  }
  con <- dbConnect(SQLite(), dbname = db_path)
  return(con)
}

con_lung <- connect_existing_db("Lung")
con_brain <- connect_existing_db("Brain_Cortex") 
con_blood <- connect_existing_db("Whole_Blood")


get_p_value_batch <- function(con, feature_variants, batch_size = 10000) {
  
  batches <- split(feature_variants, ceiling(seq_along(feature_variants) / batch_size))
  
  results <- mclapply(batches, function(batch) {
    query <- paste0("SELECT p_value FROM GTEX_data WHERE feature_variant IN ('", 
                    paste(batch, collapse = "','"), "')")
    dbGetQuery(con, query)
  }, mc.cores = 40)
  
  results <- do.call(rbind, results)
  return(results)
}

cat("肺组织数据库表:", dbListTables(con_lung), "\n")
cat("大脑皮层数据库表:", dbListTables(con_brain), "\n") 
cat("全血数据库表:", dbListTables(con_blood), "\n")

# lung GETX
msqe <- readRDS("/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")
setwd("/datg/xuxiaopeng/sc_eQTL/04_eQTL_replication")
for (Cell_type in colnames(msqe)[1:34]) {
  tryCatch({
    tmp <- msqe[, Cell_type]
    mat <- assay(tmp, "lfsrs")
    mat <- as.data.frame(mat)
    
    mat$rowname <- rownames(mat)
    
    result <- mat %>%
      filter(mat[[Cell_type]] <= 0.05)
    
    feature_SNP_pairs <- result$rowname
    
    split_pairs <- strsplit(feature_SNP_pairs, "\\|")
    gene_ids <- sapply(split_pairs, `[`, 1)
    rs_ids <- sapply(split_pairs, `[`, 2)
    
    snp_infos <- get_snp_info(rs_env, rs_ids)
    
    feature_snp <- paste(gene_ids, snp_infos, sep = "_")
    
    valid_indices <- !is.na(snp_infos)
    feature_snp <- feature_snp[valid_indices]
    
    if (length(feature_snp) == 0) {
      message("Cell_type ", Cell_type, ": 没有有效的SNP信息，跳过")
      next
    }
    
    p_value_result <- get_p_value_batch(con_lung, feature_snp)
    p_values <- p_value_result$p_value
    
    if (any(is.na(p_values)) || any(is.infinite(p_values))) {
      stop("p_values contains missing or infinite values")
    }
    
    # π1
    pi0_result <- pi0est(p_values)
    pi1 <- 1 - pi0_result$pi0
    print(paste("Cell_type:", Cell_type, "π1:", pi1))
    
  }, error = function(e) {
    message("Error in ", Cell_type, ": ", e$message)
  })
}

# brain GTEX
for (Cell_type in colnames(msqe)[1:34]) {
  tryCatch({
    tmp <- msqe[, Cell_type]
    mat <- assay(tmp, "lfsrs")
    mat <- as.data.frame(mat)
    
    mat$rowname <- rownames(mat)
    
    result <- mat %>%
      filter(mat[[Cell_type]] <= 0.05)
    
    feature_SNP_pairs <- result$rowname
    
    split_pairs <- strsplit(feature_SNP_pairs, "\\|")
    gene_ids <- sapply(split_pairs, `[`, 1)
    rs_ids <- sapply(split_pairs, `[`, 2)
    
    snp_infos <- get_snp_info(rs_env, rs_ids)
    
    feature_snp <- paste(gene_ids, snp_infos, sep = "_")
    
    valid_indices <- !is.na(snp_infos)
    feature_snp <- feature_snp[valid_indices]
    
    if (length(feature_snp) == 0) {
      message("Cell_type ", Cell_type, ": 没有有效的SNP信息，跳过")
      next
    }
    
    p_value_result <- get_p_value_batch(con_brain, feature_snp)
    p_values <- p_value_result$p_value
    
    if (any(is.na(p_values)) || any(is.infinite(p_values))) {
      stop("p_values contains missing or infinite values")
    }
    
    # π1
    pi0_result <- pi0est(p_values)
    pi1 <- 1 - pi0_result$pi0
    print(paste("Cell_type:", Cell_type, "π1:", pi1))
    
  }, error = function(e) {
    message("Error in ", Cell_type, ": ", e$message)
  })
}

# blood GTEX
for (Cell_type in colnames(msqe)[1:34]) {
  tryCatch({
    tmp <- msqe[, Cell_type]
    mat <- assay(tmp, "lfsrs")
    mat <- as.data.frame(mat)
    
    mat$rowname <- rownames(mat)
    
    result <- mat %>%
      filter(mat[[Cell_type]] <= 0.05)
    
    feature_SNP_pairs <- result$rowname
    
    split_pairs <- strsplit(feature_SNP_pairs, "\\|")
    gene_ids <- sapply(split_pairs, `[`, 1)
    rs_ids <- sapply(split_pairs, `[`, 2)
    
    snp_infos <- get_snp_info(rs_env, rs_ids)
    
    feature_snp <- paste(gene_ids, snp_infos, sep = "_")
    
    valid_indices <- !is.na(snp_infos)
    feature_snp <- feature_snp[valid_indices]
    
    if (length(feature_snp) == 0) {
      message("Cell_type ", Cell_type, ": 没有有效的SNP信息，跳过")
      next
    }
    
    p_value_result <- get_p_value_batch(con_blood, feature_snp)
    p_values <- p_value_result$p_value
    
    if (any(is.na(p_values)) || any(is.infinite(p_values))) {
      stop("p_values contains missing or infinite values")
    }
    
    # π1
    pi0_result <- pi0est(p_values)
    pi1 <- 1 - pi0_result$pi0
    print(paste("Cell_type:", Cell_type, "π1:", pi1))
    
  }, error = function(e) {
    message("Error in ", Cell_type, ": ", e$message)
  })
}




options(repr.plot.width = 15, repr.plot.height = 8)

# results
lung_data <- data.frame(
  Cell_type = c("Adventitial_fibroblast", "Aerocyte", "Alveolar_fibroblast", "Alveolar_macrophage", 
                "Arterial", "AT1", "AT2a", "AT2b", "B_cell", "Basal", "CD8T_cell", "cDC2", 
                "Ciliated", "Classical_monocytes", "Culb_1", "Culb_2", "Fibroblast", "gCap", 
                "Goblet", "Interstitial_macrophages", "Lymphatic", "Mast_cell", "Memory_CD4_T_cell", 
                "Naive_CD4_T_cell", "Neutrophils", "NK_cell", "NKT_cell", "Non-classical_monocytes", 
                "Plasma_cell", "Proliferating_T_cells", "Transitional_AT2", "Treg_T_cell", 
                "Venous", "XCL1+_T_cell"),
  π1 = c(0.415323171703791, 0.467993136907221, 0.465677175352036, 0.457979054018873, 
         0.508638551965176, 0.509391417119317, 0.481730689149544, 0.406450343830657, 
         0.469028967563823, 0.456298933274565, 0.444399262790264, 0.446268799907312, 
         0.400066559358764, 0.480053938717471, 0.468300713105705, 0.483400710917987, 
         0.484166700622349, 0.452707798460372, 0.474827335164875, 0.496151008215148, 
         0.555629886813878, 0.439600037860419, 0.431603196964044, 0.449194809403993, 
         0.430427177941855, 0.452408631866278, 0.421383109739577, 0.524583703787849, 
         0.462490722355722, 0.443423113020753, 0.442815955873262, 0.405024047981034, 
         0.454504984342943, 0.465156959110821),
  Organ = "Lung"
)

brain_data <- data.frame(
  Cell_type = c("Adventitial_fibroblast", "Aerocyte", "Alveolar_fibroblast", "Alveolar_macrophage", 
                "Arterial", "AT1", "AT2a", "AT2b", "B_cell", "Basal", "CD8T_cell", "cDC2", 
                "Ciliated", "Classical_monocytes", "Culb_1", "Culb_2", "Fibroblast", "gCap", 
                "Goblet", "Interstitial_macrophages", "Lymphatic", "Mast_cell", "Memory_CD4_T_cell", 
                "Naive_CD4_T_cell", "Neutrophils", "NK_cell", "NKT_cell", "Non-classical_monocytes", 
                "Plasma_cell", "Proliferating_T_cells", "Transitional_AT2", "Treg_T_cell", 
                "Venous", "XCL1+_T_cell"),
  π1 = c(0.0400256195633795, 0.130737341197471, 0.101944046242764, 0.0925973348749397, 
         0.108723038496567, 0.226276214357557, 0.103942263814084, 0.125007901138545, 
         0.0624051968488575, 0.0943213459104956, 0.112906768895547, 0.077124254770948, 
         0.173658245390001, 0.08812229720135, 0.0522129620536158, 0.0761469441269791, 
         0.285303562223223, 0.0409093524309757, 0.0516691244311922, 0.0736762750275849, 
         0.122108718022203, 0.015459469182655, 0.113867639183719, 0.0590449086286327, 
         0, 0.117183089708092, 0.0883011441562775, 0.0913855572630612, 
         0.00535034032214132, 0.110543875132059, 0.0210460044098599, 0.0319067561031985, 
         0.0809669356895413, 0.0735825004211649),
  Organ = "Brain"
)

blood_data <- data.frame(
  Cell_type = c("Adventitial_fibroblast", "Aerocyte", "Alveolar_fibroblast", "Alveolar_macrophage", 
                "Arterial", "AT1", "AT2a", "AT2b", "B_cell", "Basal", "CD8T_cell", "cDC2", 
                "Ciliated", "Classical_monocytes", "Culb_1", "Culb_2", "Fibroblast", "gCap", 
                "Goblet", "Interstitial_macrophages", "Lymphatic", "Mast_cell", "Memory_CD4_T_cell", 
                "Naive_CD4_T_cell", "Neutrophils", "NK_cell", "NKT_cell", "Non-classical_monocytes", 
                "Plasma_cell", "Proliferating_T_cells", "Transitional_AT2", "Treg_T_cell", 
                "Venous", "XCL1+_T_cell"),
  π1 = c(0.451790680537154, 0.342080107177766, 0.494378509152794, 0.504505982777369, 
         0.405976384693351, 0.459522092945285, 0.451460054041756, 0.425706872917369, 
         0.525576790346106, 0.442956211872178, 0.4318524042673, 0.481569019703448, 
         0.40408923370605, 0.472404374046266, 0.50139395492238, 0.496628305031047, 
         0.418935343953149, 0.395996225384294, 0.498658984416993, 0.509592219406703, 
         0.46880583148874, 0.484248398323095, 0.436761218240659, 0.442253468444704, 
         0.598902494097319, 0.434201074662358, 0.454290133791743, 0.493768614126113, 
         0.435637475611869, 0.436684888159399, 0.502117584673978, 0.505433139296956, 
         0.390013546258902, 0.501305962180901),
  Organ = "Blood"
)

combined_data <- bind_rows(lung_data, brain_data, blood_data)

cell_order <- lung_data %>% arrange(desc(π1)) %>% pull(Cell_type)
combined_data$Cell_type <- factor(combined_data$Cell_type, levels = cell_order)
combined_data$Organ <- factor(combined_data$Organ, levels = c("Lung","Blood","Brain"))

mean_values <- combined_data %>%
  group_by(Organ) %>%
  summarise(mean_pi1 = mean(π1, na.rm = TRUE))

organ_colors <- c("Lung" = "#228B22",
                  "Blood" = "#606060",
                  "Brain" = "#D3D3D3")

ggplot(combined_data, aes(x = Cell_type, y = π1, color = Organ)) +
  geom_line(aes(group = Cell_type), linetype = "solid", 
            linewidth = 0.4, alpha = 0.6, color = "gray") +
  scale_y_continuous(limits = c(0,0.8)) +
  geom_point(size = 1.4) +
  # geom_hline(data = mean_values, aes(yintercept = mean_pi1, color = Organ), 
  #            linetype = "dashed", linewidth = 0.5) + 
  scale_colour_manual(values = organ_colors) +
  labs(
    y = expression(pi[1]),
    color = "GTEX tissue"
  ) +
  mytheme + 
  guides(
    color = guide_legend(
      position = "inside",
      keywidth = 0.3,
      keyheight = 0.2,
      theme = theme(
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)
      )
    )
  ) + 
  theme(
    axis.text.x = element_text(size = 4.5, angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.title = element_text(margin = margin(b = 2)),
    legend.text = element_text(margin = margin(l = 2)),
    legend.key.spacing.y = unit(2, "pt"),
    legend.justification = c("right", "top"),
    legend.background = element_rect(
      fill = alpha("grey90", 0.8),
      color = NA
    ),
    axis.text.y = element_text(size = 5.5),
    axis.title.y = element_text(size = 7),
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
ggsave("eGene_replication.pdf", width = 11, height = 5.2, unit="cm")



