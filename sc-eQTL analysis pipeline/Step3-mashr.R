
############################################################################################################
########################################## setp1: merge qtl ################################################
############################################################################################################
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(vroom)
  library(collapse)
  library(rhdf5)
})

dir <- "/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result"

files <- list.files(path = dir, pattern = "\\.txt$", 
                    full.names = TRUE, recursive = FALSE, 
                    ignore.case = TRUE)

nCelltypes <- length(files)

all_qtl <- vroom(files, show_col_types = FALSE, id="path",
                 col_select=list(celltype="path", feature_id = "phenotype_id", variant_id = "variant_id", 
                                 betas = "slope", error = "slope_se"))

all_qtl <- all_qtl %>%
  # drop duplicates that exist if non-biallelic variants were tested
  funique(cols=c("celltype", "feature_id", "variant_id")) %>% 
  fmutate(celltype = gsub("\\..*$", "", basename(celltype)),
          id = paste(feature_id, variant_id, sep="|"))
head(all_qtl)


betas <- all_qtl %>% pivot_wider(names_from = celltype, values_from = betas,
                                 id_cols = id) %>%
  tibble::column_to_rownames(var = "id") %>% qDF()

error <- all_qtl %>% pivot_wider(names_from = celltype, values_from = error,
                                 id_cols = id) %>%
  tibble::column_to_rownames(var = "id") %>% qDF()

message("Snapshot of merged data...")
n <- ifelse(ncol(betas) > 5, 5, ncol(betas))
message("betas:")
betas[1:5, 1:n]

message("beta standard error:")
error[1:5, 1:n]

### Save output
message("Saving merged results as an hdf5...")

save_file <- "/datg/xuxiaopeng/sc_eQTL/02_mashr/all_tests.h5"
h5createFile(save_file)

ncol <- ifelse(ncol(betas) >= 10, 10, ncol(betas))
nrow <- ifelse(nrow(betas) > 1e4, nrow(betas)/100, nrow(betas/10))

h5createDataset(file = save_file, dataset = "betas", dims = dim(betas), 
                chunk = c(1000, ncol))
h5createDataset(file = save_file, dataset = "error", dims = dim(error), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas), save_file, "betas")
h5write(as.matrix(error), save_file, "error")
h5write(rownames(betas), save_file, "rownames")
h5write(colnames(betas), save_file, "colnames")

message("Done!")

############################################################################################################
################################ setp2: get_strong_and_random_qtl ##########################################
############################################################################################################

message(getRversion())
suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(argparse)
  library(data.table)
  library(collapse)
  library(rhdf5)
  library(ggplot2)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

message("Loading merged h5 data...")

all_tests_h5 <- "/datg/xuxiaopeng/sc_eQTL/02_mashr/all_tests.h5"

# 01_merge_qtl
rownames <- h5read(all_tests_h5, "rownames")
colnames <- h5read(all_tests_h5, "colnames")

betas <- h5read(all_tests_h5, "betas")
row.names(betas) <- rownames
colnames(betas) <- colnames

error <- h5read(all_tests_h5, "error")
row.names(error) <- rownames
colnames(error) <- colnames

## strong tests, top cis-eQTL for ecah gene in each cell type and combined all the top cis-sQTLs as the strong tests
message("Get significant tests...")

dir <- "/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/cis_result"

files <- list.files(path = dir, pattern = "\\.txt$", 
                    full.names = TRUE, recursive = FALSE, 
                    ignore.case = TRUE)

sig_dirs <- list.dirs(args$dir, recursive = FALSE)
sig_dirs <- paste0(sig_dirs, "/top_qtl_results_all_FDR", args$fdr, ".txt")


sig_tests <- vroom(sig_dirs, show_col_types = FALSE, progress = FALSE,
                   col_select=list(feature_id = "feature_id",
                                   variant_id = "snp_id")) %>% 
  distinct() %>% fmutate(id = paste(feature_id, variant_id, sep="|"))


if(args$verbose){message("Saving ", nrow(sig_tests), " strong QTL.")}
betas_strong <- betas[sig_tests$id, ]
error_strong <- error[sig_tests$id, ]

save_strong <- gsub("all_tests", "strong", args$h5)
h5createFile(save_strong)
ncol <- ifelse(ncol(betas_strong) >= 10, 10, ncol(betas_strong))
nrow <- ifelse(nrow(betas_strong) > 1e4, nrow(betas_strong)/100, nrow(betas_strong/10))
h5createDataset(file = save_strong, dataset = "betas", dims = dim(betas_strong), 
                chunk = c(1000, ncol))
h5createDataset(file = save_strong, dataset = "error", dims = dim(error_strong), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas_strong), save_strong, "betas")
h5write(as.matrix(error_strong), save_strong, "error")
h5write(rownames(betas_strong), save_strong, "rownames")
h5write(colnames(betas_strong), save_strong, "colnames")


if(args$verbose){message("Saving ", args$nrandom, " random QTL.")}
random <- sample(1:nrow(betas), args$nrandom)
betas_random <- betas[random, ]
error_random <- error[random, ]

save_random <- gsub("all_tests", "random", args$h5)
h5createFile(save_random)
ncol <- ifelse(ncol(betas_random) >= 10, 10, ncol(betas_random))
nrow <- ifelse(nrow(betas_random) > 1e4, nrow(betas_random)/100, nrow(betas_random/10))
h5createDataset(file = save_random, dataset = "betas", dims = dim(betas_random), 
                chunk = c(1000, ncol))
h5createDataset(file = save_random, dataset = "error", dims = dim(error_random), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas_random), save_random, "betas")
h5write(as.matrix(error_random), save_random, "error")
h5write(rownames(betas_random), save_random, "rownames")
h5write(colnames(betas_random), save_random, "colnames")

message("Done!")

############################################################################################################
###################################### setp3: mashr_fit ####################################################
############################################################################################################
message(getRversion())
suppressPackageStartupMessages({
  library(rhdf5)
  library(mashr)
  library(ashr)
  library(data.table)
  library(tidyverse)
  library(argparse)
  library(matrixcalc)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

dir <- "/datg/xuxiaopeng/sc_eQTL/02_mashr"
strongFiles <- list.files(path = dir, pattern = ".*strong.h5", 
                          full.names = TRUE)
strongFiles

for(f in strongFiles){
  message("Merging: ", f)
  
  strongTmp <- fun_h5_2_mashr(f)
  
  frand <- gsub("strong", "random", f)
  randomTmp <- fun_h5_2_mashr(frand)
  
  strongTmp$Shat <- abs(strongTmp$Shat)
  randomTmp$Shat <- abs(randomTmp$Shat)
  
  if(!exists("data.strong")){
    data.strong <- strongTmp
    data.random <- randomTmp
  }else {
    data.strong <- mash_set_data(rbind(data.strong$Bhat, strongTmp$Bhat),
                                 rbind(data.strong$Shat, strongTmp$Shat)) 
    data.random <- mash_set_data(rbind(data.random$Bhat, randomTmp$Bhat),
                                 rbind(data.random$Shat, randomTmp$Shat))
  }
}

message("Total number of strong associations: ", nrow(data.strong$Bhat))
message("Total number of random associations: ", nrow(data.random$Bhat))

keep.random <- seq(1:nrow(data.random$Bhat))
keep.strong <- seq(1:nrow(data.strong$Bhat))

# Estimate the correlation structure in the null tests from the random data
message("Estimating covariates...")

# STEP1: Learn correlation structure among null tests using random test
Vhat <- estimate_null_correlation_simple(data.random)

random <- mash_set_data(data.random$Bhat[keep.random, ],
                        data.random$Shat[keep.random, ], 
                        V=Vhat) 

strong <- mash_set_data(data.strong$Bhat[keep.strong, ], 
                        data.strong$Shat[keep.strong, ], 
                        V=Vhat)


# STEP2: Learn data-driven covariance matrices using strong tests
if(ncol(data.strong$Bhat) < 5){
  n <- ncol(data.strong$Bhat)
} else{
  n <- 5
}
covar.pca <- cov_pca(strong, n)
covar.ed <- cov_ed(strong, covar.pca)
covar.c <- cov_canonical(random)


# STEP3: Fit the mashr model to the random tests, to learn the mixture weights on all the different covariance matrices and scaling coefficients
message("Fitting mashr...")

mashFit <- mash(random, Ulist = c(covar.c, covar.ed), outputlevel = 1) 

save <- paste0(dir, "/mashr_fit.rds")
saveRDS(mashFit, save)

message("Done!")

############################################################################################################
###################################### setp4: mashr_apply ##################################################
############################################################################################################

message(getRversion())
suppressPackageStartupMessages({
  library(ggplot2)
  library(rhdf5)
  library(mashr)
  library(ashr)
  library(data.table)
  #library(tidyverse)
  library(argparse)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

message("Loading arguments...")

model <- "/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_fit.rds"
data <- "/datg/xuxiaopeng/sc_eQTL/02_mashr/all_tests.h5"
max_missing <- 0
string <- ".*0.1.h5"

mashFit <- readRDS(model)

message("mixture proportions for different types of covariance matrix:")
print(get_estimated_pi(mashFit))

if(file.exists(data) && !dir.exists(data)){
  message("Applying fit mashr model to single input file...")
  
  if(grepl(".h5", data)){
    mashData <- fun_h5_2_mashr(data, max.missing=max_missing)
  }
  
  if(grepl(".rds", data)){
    mashData <- readRDS(data)
  }
  
} else if (dir.exists(data)) {
  
  listFiles <- list.files(path = data, pattern = string,
                          full.names = TRUE)
  message(length(listFiles), " files to merge.")
  
  for(f in listFiles){
    message("Merging: ", gsub(".*/", "", f))
    
    
    beta <- load_h5(f, "betas")
    betaSE <- load_h5(f, "betaSE")
    
    if(!exists("mashData")){
      mashData <- fun_h5_2_mashr(f)
    } else{
      tmp <- fun_h5_2_mashr(f)
      mashData <- mash_set_data(rbind(mashData$Bhat, tmp$Bhat),
                                rbind(mashData$Shat, tmp$Shat)) 
    }
  }
}

save_base <- "/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_chunk"


library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

features <- unique(gsub("\\|.*", "", row.names(mashData$Bhat)))
feature_chunks <- split(features, ceiling(seq_along(features)/50))

foreach(i = seq_along(feature_chunks), .packages = "mashr") %dopar% {
  chunk_features <- feature_chunks[[i]]
  
  pattern <- paste0("^(", paste(chunk_features, collapse = "|"), ")\\|")
  keep <- grepl(pattern, row.names(mashData$Bhat))
  
  mashData_chunk <- mash_set_data(
    mashData$Bhat[keep, ],
    mashData$Shat[keep, ],
    V = mashData$V
  )
  
  m2_chunk <- mash(mashData_chunk, 
                   g = get_fitted_g(mashFit), 
                   fixg = TRUE)
  
  save_path <- paste0(save_base, names(feature_chunks)[i], ".rds")
  saveRDS(m2_chunk, save_path)
}

stopCluster(cl)

############################################################################################################
################################ setp5: mashr_get_significant ##############################################
############################################################################################################
message(getRversion())
suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(tidyverse)
  library(argparse)
  library(collapse)
  library(QTLExperiment)
  library(multistateQTL)
  library(mashr)
  library(ashr)
  source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
})

res_list <- list.files("/datg/xuxiaopeng/sc_eQTL/02_mashr/", recursive = TRUE, full.names = TRUE)
res_list <- res_list[grepl("_applied_chunk", res_list)]

for(r in 1:length(res_list)){
  message("running chunk ", r)
  tmp <- readRDS(res_list[[r]])
  tmp <- QTLExperiment::mash2qtle(tmp, sep="\\|")
  # tmp <- multistateQTL::callSignificance(tmp, thresh = 0.05, assay = "lfsrs")
  # tmp <- multistateQTL::getSignificant(tmp)
  mat <- assay(tmp, "lfsrs")
  filtered_mat <- mat[apply(mat, 1, function(row) {any(row <= 0.05)}), ]
  tmp <- tmp[rownames(filtered_mat), ]
  message("# sig in ", basename(res_list[[r]]), ": ", nrow(tmp))
  
  if(r == 1){
    msqe <- tmp
  } else{
    duplicate_rows <- intersect(rownames(msqe), rownames(tmp))
    tmp_filtered <- tmp[!rownames(tmp) %in% duplicate_rows, ]
    msqe <- QTLExperiment::rbind(msqe, tmp_filtered)
  }
}

saveRDS(msqe, "/datg/xuxiaopeng/sc_eQTL/02_mashr/mashr_applied_significant.rds")

message("Done!")



































