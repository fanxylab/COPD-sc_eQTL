library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
source("/datg/xuxiaopeng/sc_eQTL/mapping/code/project_functions.R")
# read genotype data
genotype <- fread("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/genotype/COPD_SNP.txt", sep="\t")

# define cell type, gene name and rsID
cell_type <- "Transitional AT2"
gene_id <- "ENSG00000198919.13"
gene_symbol <- "DZIP3"
rsID <- "rs6437823"

# select genotype for rsID
rsID_geno <- genotype %>% filter(ID == rsID)

# read phenotype data
phenotype <- fread(paste("/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/phenotype/", 
                         cell_type, ".bed", sep="" ), sep="\t")
gene_pheno <- phenotype %>% filter(phenotype_id == gene_id)

# 
common_samples <- intersect(
  names(rsID_geno)[6:ncol(rsID_geno)],
  names(gene_pheno)[5:ncol(gene_pheno)]
)

dt1_long <- melt(
  rsID_geno,
  id.vars = c("CHROM", "POS", "ID", "REF", "ALT"),
  measure.vars = common_samples,
  variable.name = "Sample",
  value.name = "Genotype"
)[, .(Sample, Genotype)]

dt2_long <- melt(
  gene_pheno,
  id.vars = c("#chr", "start", "end", "phenotype_id"),
  measure.vars = common_samples,
  variable.name = "Sample",
  value.name = "Expression"
)[, .(Sample, Expression)]

merged_dt <- merge(
  dt1_long,
  dt2_long,
  by = "Sample",
  all = FALSE
)

merged_dt[, Genotype := {
  alleles <- strsplit(Genotype, "\\|")[[1]]
  sorted_alleles <- sort(alleles)
  paste(sorted_alleles, collapse = "/")
}, by = Genotype]

ref_allele <- rsID_geno$REF[1]
alt_allele <- rsID_geno$ALT[1]

geno_ref_ref <- paste(sort(c(ref_allele, ref_allele)), collapse = "/")
geno_het     <- paste(sort(c(ref_allele, alt_allele)), collapse = "/")
geno_alt_alt <- paste(sort(c(alt_allele, alt_allele)), collapse = "/")

target_levels <- c(geno_ref_ref, geno_het, geno_alt_alt)
actual_levels <- target_levels[target_levels %in% unique(merged_dt$Genotype)]

merged_dt[, Genotype := factor(Genotype, levels = actual_levels)]

merged_dt[, Disease := ifelse(grepl("^C", Sample), "COPD", "Healthy")]


label_data <- merged_dt[, .(n = .N), by = .(Genotype, Disease)]
regression_data <- merged_dt[, {
  lm_fit <- lm(Expression ~ as.numeric(as.factor(Genotype)))
  data.table(
    Genotype = levels(as.factor(Genotype)),
    Predicted = predict(lm_fit, newdata = data.table(Genotype = levels(as.factor(Genotype))))
  )
}, by = Disease]

ggplot(merged_dt, aes(Genotype, Expression)) +
  geom_violin(
    aes(fill = Genotype), 
    alpha = 0.6, 
    trim = FALSE,
    scale = "width",
    width = 0.6,
    color = NA
  ) +
  geom_boxplot(
    aes(color = Genotype), 
    width = 0.20,
    alpha = 0.8,
    linewidth = 0.1, 
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = Genotype),
    position = position_jitterdodge(
      jitter.width = 0.15
    ),
    shape = 21,
    stroke = 0,
    fill = "gray10",
    size = 0.2
  ) +
  geom_text(
    data = label_data,
    aes(x = Genotype, y = -Inf, label = sprintf("(n=%d)", n)),
    vjust = -0.3,
    color = "black",
    size = 1.5,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = regression_data,
    aes(x = as.factor(Genotype), y = Predicted, group = Disease),
    color = "red",
    linewidth = 0.1,
    inherit.aes = FALSE
  ) +
  labs(
    x = rsID,
    y = bquote(paste(italic(.(gene_symbol)), " expression"))
  ) +
  scale_fill_manual(values = c("#37b2cb", "#9baa66", "#ffa200")) +
  scale_colour_manual(values = rep("gray50", 3)) +
  mytheme + theme(
    legend.position = "none",
    
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    
    axis.line.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.line.y.left = element_line(colour="black",linewidth=0.1),
    axis.ticks.x.bottom = element_line(colour="black",linewidth=0.1),
    axis.ticks.y.left = element_line(colour="black",linewidth=0.1),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 5, hjust = 0.5, 
                              face = "plain",margin = margin(t = 0, b = 0)),
    strip.placement = "outside"
  ) +
  facet_wrap(~ Disease, strip.position = "top")

ggsave(paste("/datg/xuxiaopeng/sc_eQTL/eQTL_boxplot/", gene_symbol, "_", 
             cell_type, "_disease.pdf", sep=""), 
       width = 5, height = 3.4, unit = "cm")

