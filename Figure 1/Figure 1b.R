library(ggplot2)
library(tidyr)
library(dplyr)

sample_info <- read.table("/datg/xuxiaopeng/sc_eQTL/COPD/Sample_info.txt", 
                          sep="\t", header=TRUE)
sample_info$Disease <- ifelse(startsWith(sample_info$sample_id, "C"), "Case", 
                              ifelse(startsWith(sample_info$sample_id, "H"), "Control", NA))
sample_info$Ethnicity <- "East Asian"
sample_info$Smoke <- ifelse(startsWith(sample_info$Group, "HC"), "No", "Yes")

donors <- read.table("/datg/xuxiaopeng/sc_eQTL/COPD/donors_final.txt", 
                     sep="\t", header=FALSE)
sample_info <- sample_info %>% filter(sample_id %in% donors$V1)
sample <- sample_info[, c("sample_id", "Sex", "Age", 
                          "Disease", "Ethnicity", "Smoke")]
sample <- sample %>%
  mutate(Age = cut(Age, breaks = seq(20, 100, by = 20), right = FALSE, 
                   labels = c("20-39",  "40-59", "60-79",  "80-99")))

long_data <- sample %>%
  pivot_longer(cols = c(Sex, Age, Disease, Ethnicity, Smoke), 
               names_to = "Variable", 
               values_to = "Value")

count_data <- long_data %>%
  group_by(Variable, Value) %>%
  summarise(Count = n(), .groups = "drop")

count_data$Value <- factor(count_data$Value, 
                           levels=c("80-99", "60-79", "40-59", "20-39", 
                                    "Case", "Control", "East Asian", 
                                    "Male", "Female", "Yes", "No")
                           )

count_data$Variable <- factor(count_data$Variable, 
                              levels=c("Sex", "Age", "Disease", 
                                       "Smoke", "Ethnicity")
                              )

ggplot(count_data, aes(x = Variable, y = Count, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Value), 
            position = position_stack(vjust = 0.5),
            color = "black", size = 2) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(y = "Number of donors") +
  mytheme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
    # axis.text.x = element_text(angle = 0, hjust = 1)
  )

setwd("/datg/xuxiaopeng/sc_eQTL/Graph")
ggsave("Donor_info.pdf", width = 6.7, height = 5, unit="cm")
