
# load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")

min_anova_df <- read.table("dynamic_eQTL_results_2.txt", sep = "\t", header = TRUE)

min_anova_df <- min_anova_df %>% 
  filter(FDR <= 0.05)

gene_list <- min_anova_df$phenotye

# GO analysis
ego <- enrichGO(gene          = gene_list,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

ego_result <- as.data.frame(ego)
ego_result <- ego_result[order(ego_result$pvalue), ]

top10 <- head(ego_result, 10)


# plot
ggplot(top10, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "#69b3a2", width = 0.7) +
  geom_text(
    aes(label = Description, 
        y = 0.1),
    hjust = 0,
    color = "gray20",
    size = 1.7,
    fontface = "bold"
  ) +
  coord_flip() +
  labs(
    x = "GO term",
    y = "-log10(p-value)"
  ) +
  mytheme +
  theme(
    axis.title.x = element_text(size = 7),  
    axis.title.y = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.1, colour = "black"),
    axis.ticks.x.bottom = element_line(linewidth = 0.1, colour = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.1, colour = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)))


setwd("/datg/xuxiaopeng/sc_eQTL/06_dynamic")
ggsave("GO_enrichment_2.pdf", width = 5.5, height = 5.5, unit = "cm")
write.table(ego_result, "GO_enrichment.txt", 
            sep="\t", quote=FALSE, col.names = TRUE, row.names=TRUE)

