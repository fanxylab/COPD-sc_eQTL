library(coloc)
library(snpStats)
library(locuscomparer)
library(data.table)
library(dplyr)


# initialize final results dataframe
all_results <- NULL

# gwas information
gwas_21_NG <- c("Asthma.tsv","Bronchiectasis.tsv","Bronchitis.tsv","Chronic_bronchitis.tsv","COPD.tsv","ILD.tsv",
                "Lung_cancer.tsv","PF.tsv","Pneumonia.tsv","Pneumothorax.tsv","PT.tsv","Sarcoidosis.tsv")
# N_cases <- c(51384,3129,5948,10430,17547,3313,8235,1692,24310,2439,8695,1938)
# N_controls <- c(574064,602066,644818,602066,617598,644534,663294,644534,634715,652285,647362,662622)

gwas_23_NG <- c("FEV1.tsv","FVC.tsv","FEV1-FVC.tsv","PEF.tsv")

gwas_22_CG <- c("Asthma-CG.tsv","COPD-CG.tsv","IPF-CG.tsv")
# N_cases <- c(153763,81568,8006)
# N_controls <- c(1647022,1310798,1246742)

gwas_list <- list(gwas_21_NG, gwas_23_NG, gwas_22_CG)

# cell type information
cell_types = c('gCap','NKT_cell','Memory_CD4_T_cell','Naive_CD4_T_cell','Non-classical_monocytes','AT2a','CD8T_cell','NK_cell','cDC2','AT1',
               'Classical_monocytes','Aerocyte','Adventitial_fibroblast','Arterial','Alveolar_macrophage','Culb_1','Venous','Interstitial_macrophages',
               'Lymphatic','AT2b','Culb_2','Mast_cell','Goblet','Ciliated','Transitional_AT2','XCL1+_T_cell','B_cell','Neutrophils','Treg_T_cell', 
               'Plasma_cell','Alveolar_fibroblast','Proliferating_T_cells','Fibroblast','Basal', "GTEX_v10_Lung_converter")

for(i in 2:2){
  
  gwas_files <- gwas_list[[i]]
  
  for(j in 3:3){
    gwas_file <- gwas_files[j]
    gwas_id <- strsplit(gwas_file, "\\.")[[1]][1]
    gwas <- fread(paste("/datg/xuxiaopeng/sc_eQTL/07_GWAS/processed_sumstats2/", gwas_id, ".tsv", sep=""), sep = "\t")
    
    if(i==1){
      gwas$N_CASE <- N_cases[j]
      gwas$N_CONTROL <- N_controls[j]           
    }
    
    #################################################
    ###################  step 1  ####################
    #################################################
    gwas_1e_7 <- gwas %>% 
      filter(P < 1e-7)
    pval <- gwas_1e_7[, c("CHR", "POS", "P")]
    colnames(pval) <- c("chr", "pos", "P")
    
    pval <- pval[!(chr == 6 & pos >= 25e6 & pos <= 35e6)] 
    setorder(pval, chr, P)
    
    gwas_loci <- pval[, {
      keep_pos <- numeric(0)
      candidate_pos <- .SD$pos
      candidate_p <- .SD$P
      
      while(length(candidate_pos) > 0) {
        selected_idx <- which.min(candidate_p)
        selected_pos <- candidate_pos[selected_idx]
        keep_pos <- c(keep_pos, selected_pos)
        exclude_range <- c(selected_pos - 5e5, selected_pos + 5e5)
        valid_idx <- which(candidate_pos < exclude_range[1] | candidate_pos > exclude_range[2])
        candidate_pos <- candidate_pos[valid_idx]
        candidate_p <- candidate_p[valid_idx]
      }
      .SD[which(pos %in% keep_pos), .(pos, P)]
    }, by = chr]
    
    setnames(gwas_loci, c("chr", "pos", "P"))
    setcolorder(gwas_loci, c("chr", "pos", "P"))
    setkey(gwas_loci, chr, pos)
    
    write.table(gwas_loci, paste("/datg/xuxiaopeng/sc_eQTL/07_GWAS/processed_sumstats2/", gwas_id,"_loci.txt",sep=""),
                sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    #################################################
    ###################  step 2  ####################
    #################################################
    gwas_1e_7 <- as.data.frame(gwas_1e_7)
    GWAS_pre <- gwas_1e_7[, c("rsID", "P")]
    colnames(GWAS_pre) <- c("variant_id", "p_value")
    
    for(k in 1:length(cell_types)){
      celltype <- cell_types[k]
      
      nominal_full <- fread(paste("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result/",celltype,".nominal.all.txt",sep=""), sep = "\t")
      nominal_full <- as.data.frame(nominal_full)
      
      nominal_sig <- nominal_full %>% filter(pval_nominal < 1e-3)
      
      common_snp <- merge(nominal_sig, GWAS_pre, by = "variant_id", all = FALSE)
      common_snp <- as.data.table(common_snp)
      
      common_snp <- common_snp[!(chrom == 6 & pos >= 25e6 & pos <= 35e6)]
      setorder(common_snp, chrom, p_value)
      
      common_loci <- common_snp[, {
        keep_pos <- numeric(0)
        candidate_pos <- .SD$pos
        candidate_p <- .SD$p_value
        while(length(candidate_pos) > 0) {
          selected_idx <- which.min(candidate_p)
          selected_pos <- candidate_pos[selected_idx]
          keep_pos <- c(keep_pos, selected_pos)
          exclude_range <- c(selected_pos - 5e5, selected_pos + 5e5)
          valid_idx <- which(candidate_pos < exclude_range[1] | candidate_pos > exclude_range[2])
          candidate_pos <- candidate_pos[valid_idx]
          candidate_p <- candidate_p[valid_idx]
        }
        .SD[which(pos %in% keep_pos), .(variant_id, phenotype_id, pos, p_value, pval_nominal)]
      }, by = chrom]
      
      #################################################
      ###################  step 3  ####################
      #################################################
      if(i==1 || i==3){
        GWAS_pre2 <- gwas[,c("rsID", "CHR", "POS", "MAF", "BETA", "SE", "P", "N_CASE", "N_CONTROL")]
        GWAS_pre2 <- as.data.frame(GWAS_pre2)                
      }
      if(i==2){
        GWAS_pre2 <- gwas[,c("rsID", "CHR", "POS", "MAF", "P", "N")]
        GWAS_pre2 <- as.data.frame(GWAS_pre2)        
      }
      
      if (nrow(common_loci)>0){
        
        temp_result <- as.data.frame(matrix(NA,nrow(common_loci),13))
        colnames(temp_result) <- c("cell_type","gwas_id","rsID","phenotype_id",
                                   "nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf",
                                   "eQTL_pval","gwas_pval", "best_H4_snp") 
        
        for(m in 1:nrow(common_loci)){
          common_loci <- as.data.frame(common_loci)
          
          temp_result[m,1] <- celltype
          temp_result[m,2] <- gwas_id
          temp_result[m,3] <- common_loci$variant_id[m]
          temp_result[m,4] <- common_loci$phenotype_id[m]
          temp_result[m,11] <- common_loci$pval_nominal[m]
          temp_result[m,12] <- common_loci$p_value[m]
          
          gene <- common_loci$phenotype_id[m]
          
          snp_1mb <- nominal_full %>%
            filter(phenotype_id == gene) %>%
            select("phenotype_id","variant_id","pval_nominal","slope","slope_se","maf","n")
          colnames(snp_1mb) <- c("phenotype_id","rsID","pval_nominal","slope","slope_se","maf","n")
          snp_1mb <- snp_1mb[order(snp_1mb$pval_nominal), ]
          
          input <- merge(snp_1mb, GWAS_pre2, by="rsID", all=FALSE, suffixes=c("_eqtl","_gwas"))
          input <- input[!duplicated(input$rsID), ]
          
          if(nrow(input)>10){
            
            if(i==2){
              dataset1 = list(snp=input$rsID, pvalues=input$pval_nominal, beta=input$slope, varbeta=(input$slope_se)^2, type="quant", N=input$n, MAF=input$maf)
              dataset2 = list(snp=input$rsID, pvalues=input$P, type="quant", N=input$N, MAF=input$MAF)                        
            }
            if(i==1 || i==3){
              dataset1 = list(snp=input$rsID, pvalues=input$pval_nominal, beta=input$slope, varbeta=(input$slope_se)^2, type="quant", N=input$n, MAF=input$maf)              
              dataset2 = list(snp=input$rsID, pvalues=input$P, beta=input$BETA, varbeta=(input$SE)^2, type="cc", 
                              s=N_cases[j]/(N_cases[j]+N_controls[j]), N=N_cases[j]+N_controls[j], MAF=input$MAF)                      
            }
            
            result <- coloc.abf(dataset1, dataset2)
            temp_result[m,5:10] <- t(as.data.frame(result$summary))[1,1:6]
            
            snp_results <- result$results
            best_idx <- which.max(snp_results$SNP.PP.H4)
            best_H4_snp <- snp_results$snp[best_idx]
            
            temp_result[m, 13] <- best_H4_snp 
          }
        }
      }else{
        temp_result <- NULL
      }
      all_results <- rbind(all_results, temp_result)
    }
  }
}

write.table(all_results, paste("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/", gwas_id, ".coloc.txt", sep=""),
            sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE
)
