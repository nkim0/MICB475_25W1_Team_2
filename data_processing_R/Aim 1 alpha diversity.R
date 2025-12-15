library(phyloseq)
library(tidyverse)
library(RColorBrewer)

load(file = "phylseq_files.RData")

### Plotting plco observed and shannon alpha diversity ###
rich_plot_plco <- plot_richness(PLCO_genus, x="disease_alc_status", measures=c("Observed", "Shannon")) +
  geom_boxplot(aes(fill = disease_alc_status)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Cancer and Alcohol Consumption",
    y = "Alpha Diversity Measure",
    fill = "Cancer and Alcohol Consumption") +
  scale_x_discrete(labels = label) +
  scale_fill_discrete(labels = label) +
  scale_fill_manual(values = palette, labels = label)
rich_plot_plco

### Kruskal wallis test on plco alpha diversity ###
alphadiv_plco <- estimate_richness(PLCO_genus, measures = c("Observed", "Shannon")) |>
  rownames_to_column(var = "sample.id") |>
  left_join(as.data.frame(sample_data(PLCO_genus)), by = "sample.id")
head(alphadiv_plco)
kruskal.test(Shannon ~ disease_alc_status, data = alphadiv_plco)


### Plotting  cpsII observed and shannon slpah diveresity ###
rich_plot_cpsII <- plot_richness(CPSII_genus, x="disease_alc_status", measures=c("Observed", "Shannon")) +
  geom_boxplot(aes(fill = disease_alc_status)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Cancer and Alcohol Consumption",
    y = "Alpha Diversity Measure",
    fill = "Cancer and Alcohol Consumption") +
  scale_x_discrete(labels = label) +
  scale_fill_discrete(labels = label) +
  scale_fill_manual(values = palette, labels = label)
rich_plot_cpsII


### Kruskal wallis test on cpsII alpha diversity ###
alphadiv_cpsII <- estimate_richness(CPSII_genus, measures = c("Observed", "Shannon")) |>
  rownames_to_column(var = "sample.id") |>
  left_join(as.data.frame(sample_data(CPSII_genus)), by = "sample.id")
head(alphadiv_cpsII)
kruskal.test(Shannon ~ disease_alc_status, data = alphadiv_cpsII)


ggsave("rp_plco.png"
       , rich_plot_plco
       , height = 5, width = 10)

ggsave("rp_cpsII.png"
       , rich_plot_cpsII
       , height = 5, width = 10)
