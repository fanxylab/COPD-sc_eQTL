library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(ggplot2)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")

# read data
genotype <- read.table(
  "/datg/xuxiaopeng/WGS/plink/merged_ALL_SNPs_transform.txt", 
  sep="\t", header=TRUE)
pb <- readRDS("Dynamic_sce_object.rds")
sample_mapping <- read.table(
  "/datg/xuxiaopeng/sc_eQTL/mapping/Others/sample_mapping_file.txt", sep="\t")

# 
gene_name <- "RGCC"
SNP_name <- "rs12431090"

geno_temp <- genotype %>%
  filter(ID == SNP_name)
head(geno_temp)


#  Generate aggregrated mean of cells
Q1 <- as.data.frame(assay(pb, "Q1")[gene_name, ])
colnames(Q1) <- "Q1"
Q2 <- as.data.frame(assay(pb, "Q2")[gene_name, ])
colnames(Q2) <- "Q2"
Q3 <- as.data.frame(assay(pb, "Q3")[gene_name, ])
colnames(Q3) <- "Q3"
Q4 <- as.data.frame(assay(pb, "Q4")[gene_name, ])
colnames(Q4) <- "Q4"
Q5 <- as.data.frame(assay(pb, "Q5")[gene_name, ])
colnames(Q5) <- "Q5"
# Q6 <- as.data.frame(assay(pb, "Q6")[gene_name, ])
# colnames(Q6) <- "Q6"
# Q7 <- as.data.frame(assay(pb, "Q7")[gene_name, ])
# colnames(Q7) <- "Q7"
# Q8 <- as.data.frame(assay(pb, "Q8")[gene_name, ])
# colnames(Q8) <- "Q8"

expr <- cbind(Q1, Q2, Q3, Q4, Q5)
expr$sample_id <- rownames(expr)

geno_temp <- genotype %>%
  filter(ID == SNP_name)
geno_temp <- as.data.frame(t(geno_temp)[-1:-5,])
geno_temp$sample <- rownames(geno_temp)
colnames(geno_temp) <- c("SNP", "sample")
geno <- geno_temp %>% 
  left_join(sample_mapping, by = c("sample" = "V1")) %>% 
  select(sample = V2, SNP)
expr_geno <- geno %>% 
  inner_join(expr, by=c("sample" = "sample_id"))
expr_geno <- expr_geno %>%
  pivot_longer(
    cols = Q1:Q5,
    names_to = "Quantile",
    values_to = "Expression"
  )

expr_geno$SNP <- factor(expr_geno$SNP, levels = c("G/G", "G/A", "A/A"))
expr_geno$Quantile <- factor(expr_geno$Quantile)

head(expr_geno)



ggplot(expr_geno) +
  geom_boxplot(
    aes(x = Quantile, y = Expression, fill = SNP),
    linewidth = 0.05,
    fatten = 2,
    outliers = FALSE,
    position = position_dodge(width = 0.8)
  ) +
  geom_smooth(
    aes(x = X_dodge, y = Expression, group = Quantile),
    method = "lm",
    se = FALSE,
    color = "firebrick",  
    linewidth = 0.5,
    alpha = 0.9
  ) +
  geom_jitter(
    aes(x = Quantile, y = Expression, fill = SNP, group = SNP),
    position = position_jitterdodge(
      jitter.width = 0.2,   
      dodge.width = 0.8     
    ),
    shape = 21,  
    stroke = 0,
    size = 0.2
  ) +
  scale_fill_manual(values = c("#37b2cb", "#9baa66", "#ffa200")) +  
  scale_colour_manual(values = rep("gray10", 3)) +
  guides(
    fill = guide_legend(
      title = "Genotype",
      direction = "horizontal",
      theme = theme(
        legend.key.width  = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm"),
      )
    ), 
    colour = "none"  
  ) +  
  labs(
    x = SNP_name,
    y = bquote(paste(italic(.(gene_name)), " expression"))
  ) +
  mytheme +
  theme(
    axis.title.x = element_text(size = 7),
    strip.text = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 5.5, colour = "black"),
    
    legend.position = "top",  
    legend.justification = c(0, 1),  
    legend.margin = margin(0, 0, 0, 0),  
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, "pt"),
    legend.title = element_text(margin = margin(r = 2), size = 5.5),
    legend.text = element_text(margin = margin(l = 2), size = 5.5),
    legend.key.spacing.x = unit(0.2, "cm"),
    
    axis.text.y = element_text(size = 5.5, colour = "black"),
    axis.title.y = element_text(size = 7),
    
    axis.line.y.left = element_line(linewidth = 0.1, colour = "black"),
    axis.ticks.y.left = element_line(linewidth = 0.1, colour = "black"),
    axis.ticks.x.bottom = element_line(linewidth = 0.1, colour = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.1, colour = "black")
  )

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
ggsave("RGCC_pseudo_boxplot.pdf", width = 6.2, height = 2.6, unit = "cm")

