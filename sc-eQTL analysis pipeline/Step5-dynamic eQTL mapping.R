library(Seurat)
library(lme4)
library(data.table)
library(qvalue)
library(dplyr)
library(tidyr)
library(ggplot2)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
library(Seurat)
library(slingshot)
library(scCustomize)
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(lme4)
library(phateR)
library(dplyr)
library(tidyr)
library(viridis)
library(data.table)
library(Matrix)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(pdist)
library(monocle3)
library(tradeSeq)
library(BiocParallel)
library(pheatmap)

library(SeuratWrappers)
library(tidydr)

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
X <- readRDS("Dynamic.rds")
X <- JoinLayers(X)
expression_matrix <- GetAssayData(X, slot = "data")
genotype <- read.table("/datg/xuxiaopeng/WGS/plink/merged_ALL_SNPs_transform.txt", sep="\t", header=TRUE)

eQTL_pair <- read.table("eQTLs_for_dynamic_mapping.txt", sep="\t", header = TRUE)
eQTL_SNPs <- eQTL_pair$rsID
eQTL_GENEs <- eQTL_pair$gene_symbol

X@meta.data <- X@meta.data %>%
  mutate(
    Quantile = cut(
      monocle3PT,
      breaks = seq(0, 62.7, by = 7.8375),
      include.lowest = TRUE,
      right = FALSE,
      labels = paste0("Q", 1:8)         
    )
  )

X@meta.data <- X@meta.data %>%
  mutate(New_Quantile = case_when(
    Quantile %in% c("Q1", "Q2", "Q3") ~ "Q1",
    Quantile == "Q4" ~ "Q2",
    Quantile == "Q5" ~ "Q3",
    Quantile == "Q6" ~ "Q4",
    Quantile %in% c("Q7", "Q8") ~ "Q5"
  ))
X@meta.data$New_Quantile <- factor(X@meta.data$New_Quantile, levels=c("Q1","Q2","Q3","Q4","Q5"))

# transform to singlecell object
sce <- as.SingleCellExperiment(X, assay = "RNA")

# aggregate cell count using mean method
pb <- aggregateData(sce, 
                    assay = "logcounts", 
                    fun = "mean", 
                    by = c("New_Quantile", "sample")
                   )

saveRDS(pb, file = "Dynamic_sce_object.rds")

pb <- readRDS("Dynamic_sce_object.rds")

false_genes <- setdiff(eQTL_GENEs, rownames(expression_matrix))
eQTL_pair_filtered <- eQTL_pair %>% filter(gene_symbol %notin% false_genes)
eQTL_SNPs <- eQTL_pair_filtered$rsID
eQTL_GENEs <- eQTL_pair_filtered$gene_symbol

################################################################################
# 1. Import Necessary Parallel and Acceleration Packages
################################################################################
library(lme4)
library(dplyr)
library(foreach)
library(doMC)       # Shared-memory multi-processing for Linux
library(data.table) # For ultra-fast result merging

# ================= Optimize Thread Configuration =================
# Recommend using 8 or 10 threads first to ensure speed 
# while safely avoiding system Fork limits
threads <- 10 
registerDoMC(threads) 
# ============================================================

message(paste("🚀 Linux FORK shared-memory parallelization initiated, using", threads, "threads..."))

################################################################################
# 2. Extract Necessary Data
################################################################################
# Extract expression PCs (ePCs)
ePCs <- as.data.frame(Embeddings(X, reduction = "pca")[,1:8]) %>%
  select(ePC1=PC_1, ePC2=PC_2, ePC3=PC_3, ePC4=PC_4, 
         ePC5=PC_5, ePC6=PC_6, ePC7=PC_7, ePC8=PC_8)

# Extract cell-level metadata
cell_meta <- X@meta.data[, c("sample", "Diagnosis", "New_Quantile")]

# ==================== 🧹 Memory Cleanup (Critical Step) ====================
message("🧹 Destroying large Seurat object to release system memory...")
rm(X)         # Completely remove Seurat object
gc()          # Force R to release memory to the OS
gc()          # Second garbage collection to ensure clean state
# ==========================================================================

# Import external auxiliary data
sample_mapping <- read.table("/datg/xuxiaopeng/sc_eQTL/mapping/Others/sample_mapping_file.txt", sep="\t")

sample_info <- read.table("/datg/xuxiaopeng/sc_eQTL/COPD/Sample_info.txt", sep="\t", header=TRUE)
sample_info <- sample_info[,c("sample_id","Sex","Age")]

gPCs <- read.table("/datg/xuxiaopeng/sc_eQTL/05_QTLtools/genotype/COPD_all.pca", sep=" ", header=TRUE)
gPCs <- as.data.frame(t(gPCs))

gPCs <- gPCs[2:124,] %>%
  select(V1,V2,V3) %>%
  select(gPC1=V1, gPC2=V2, gPC3=V3) %>%
  mutate(sample_id = rownames(gPCs[2:124,]))

# Initialize progress log file
log_file <- "eqtl_progress.log"
writeLines(paste("Run started at:", Sys.time()), log_file)

################################################################################
# 3. Memory-Efficient Parallel Loop
################################################################################
total_pairs <- nrow(eQTL_pair_filtered)

system.time({
  anova_list <- foreach(i = 1:total_pairs, 
                        .errorhandling = "pass") %dopar% {
                          
                          gene_name <- eQTL_GENEs[i]
                          SNP_name <- eQTL_SNPs[i]
                          
                          # Write progress to log file every 100 pairs
                          if (i %% 100 == 0) {
                            log_msg <- sprintf("[%s] Completed %d / %d pairs (%.2f%%) | Currently processing: %s - %s", 
                                               Sys.time(), i, total_pairs, (i / total_pairs) * 100, gene_name, SNP_name)
                            write(log_msg, file = log_file, append = TRUE)
                          }
                          
                          # 1. Cell-level data collection
                          gene_expression <- data.frame(Expression = expression_matrix[gene_name, ])
                          cell_temp <- cbind(gene_expression, cell_meta, ePCs)
                          
                          # 2. Donor-level data collection
                          geno_temp <- genotype %>% filter(ID == SNP_name)
                          if (nrow(geno_temp) == 0) {
                            return(data.frame(phenotye = gene_name, variants = SNP_name, Pr = NA, status = "SNP_not_found"))
                          }
                          
                          geno_temp <- as.data.frame(t(geno_temp)[-1:-5,])
                          geno_temp$sample <- rownames(geno_temp)
                          colnames(geno_temp) <- c("SNP", "sample")
                          
                          geno <- geno_temp %>% 
                            left_join(sample_mapping, by = c("sample" = "V1")) %>% 
                            select(sample = V2, SNP)
                          
                          sample_temp <- geno %>% 
                            left_join(sample_info, by=c("sample" = "sample_id")) %>% 
                            left_join(gPCs, by=c("sample" = "sample_id"))
                          
                          # 3. Merge and type conversion
                          data_all <- cell_temp %>% left_join(sample_temp, by=c("sample" = "sample"))
                          
                          data_all$sample <- factor(data_all$sample)
                          data_all$Diagnosis <- factor(data_all$Diagnosis)
                          data_all$SNP <- factor(data_all$SNP)
                          data_all$Sex <- factor(data_all$Sex)
                          data_all$gPC1 <- as.numeric(data_all$gPC1)
                          data_all$gPC2 <- as.numeric(data_all$gPC2)
                          data_all$gPC3 <- as.numeric(data_all$gPC3)
                          
                          # 4. Model fitting (with optimization control parameters for convergence)
                          lmer_control <- lmerControl(
                            optimizer = "bobyqa", 
                            optCtrl = list(maxfun = 2e4),
                            calc.derivs = FALSE
                          )
                          
                          # Use tryCatch to prevent thread crash caused by data singularity
                          result <- tryCatch({
                            fit0 <- lmer(Expression ~ SNP + Age + Sex + Diagnosis + 
                                           ePC1 + ePC2 + ePC3 + ePC4 + ePC5 + ePC6 + ePC7 + ePC8 + gPC1 + gPC2 + gPC3 +
                                           (1 | sample) + New_Quantile + SNP*New_Quantile, 
                                         data = data_all, REML = FALSE, control = lmer_control)
                            
                            fit1 <- lmer(Expression ~ SNP + Age + Sex + Diagnosis + 
                                           ePC1 + ePC2 + ePC3 + ePC4 + ePC5 + ePC6 + ePC7 + ePC8 + gPC1 + gPC2 + gPC3 +
                                           (1 | sample) + New_Quantile, 
                                         data = data_all, REML = FALSE, control = lmer_control)
                            
                            lrt <- anova(fit1, fit0)
                            p_val <- lrt$`Pr(>Chisq)`[2]
                            
                            data.frame(phenotye = gene_name, variants = SNP_name, Pr = p_val, status = "OK")
                            
                          }, error = function(e) {
                            data.frame(phenotye = gene_name, variants = SNP_name, Pr = NA, status = paste("Error:", e$message))
                          })
                          
                          return(result)
                        }
})

################################################################################
# 4. Ultra-Fast Result Merging
################################################################################
message("📊 Merging calculation results...")
anova_list_cleaned <- anova_list[sapply(anova_list, is.data.frame)]
anova_df <- as.data.frame(data.table::rbindlist(anova_list_cleaned))

# Print execution status statistics
print(table(anova_df$status))


min_anova_df <- anova_df %>% 
  group_by(phenotye) %>% 
  slice_min(Pr, n = 1) %>% 
  ungroup()

min_anova_df$FDR <- p.adjust(min_anova_df$Pr, method = "bonferroni")

write.table(min_anova_df, "dynamic_eQTL_results_2.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names=FALSE)
write.table(anova_df, "dynamic_eQTL_results_ALL.txt", sep="\t", quote=FALSE, col.names = TRUE, row.names=FALSE)






















