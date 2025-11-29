library(phyloseq)
library(tidyverse)
library(ggplot2)
library(pals)

# load objects
load(file = "phylseq_files.RData")

# get the name of each analysis group
print(levels(as.factor(sample_data(plco_rare)$disease_alc_status)))

# Agglomerate at genus level
plco_genus <- tax_glom(PLCO, taxrank = "Genus", NArm = FALSE)

# Convert to relative abundance, then melt
rel_abund <- plco_genus |>
  microbiome::transform('compositional') |>
  psmelt()

# get list of genus below 1%
sum_of_each_genus <- taxa_sums(plco_genus)
avg_of_each_genus <- sum_of_each_genus/sum(sum_of_each_genus)
below_1 <- avg_of_each_genus[avg_of_each_genus<0.01]
to_change <- names(below_1)


# Find the average abundance for each genus within disease_alc_status
average_abund <- rel_abund %>% 
  group_by(OTU,Genus,disease_alc_status) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  ungroup()

# rename Genus under 1% as "<1%
avg_modified <- average_abund
avg_modified$Genus[avg_modified$OTU %in% to_change] <- 'others (<1%)'

# Plot taxonomy graph at genus level with disease_alc_status on X-axis
tax_barplot <- avg_modified|> 
  mutate(Genus = str_remove(Genus, "^g__")) |>
  mutate(Genus = fct_relevel(Genus, "others (<1%)")) |>
  filter(!is.na(Genus)) |>
  ggplot(aes(disease_alc_status, Abundance, fill=Genus)) +
  geom_col(position='stack') +
  labs(
    x = "Cancer Status: Alcohol Consumption",
    y = "Relative Abundance",
    fill = "Genus") +
  theme_classic() +
  scale_fill_manual(values = as.vector(stepped()))
tax_barplot

# save plots
ggsave("plco_genus_barplot.png",
       tax_barplot,
       height = 5, width = 7)
