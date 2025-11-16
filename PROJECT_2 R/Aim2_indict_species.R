library(phyloseq)
library(indicspecies)
library(tidyverse)

# load objects
load(file = "phylseq_files.RData")

### AIM 2.2: Indicator Species Analysis ###
# Agglomerate at genus level and convert to relative abundance
PLCO_genus <- tax_glom(PLCO, taxrank = "Genus", NArm = FALSE)
PLCO_genus_RA <- transform_sample_counts(PLCO_genus, fun=function(x) x/sum(x))

CPSII_genus <- tax_glom(CPSII, taxrank = "Genus", NArm = FALSE)
CPSII_genus_RA <- transform_sample_counts(CPSII_genus, fun=function(x) x/sum(x))

# Run indicator species analysis using the donor_status variable and look at results
isa_plco <- multipatt(t(otu_table(PLCO_genus_RA)), cluster = sample_data(PLCO_genus_RA)$'disease_alc_status')
summary(isa_plco)

isa_cpsII <- multipatt(t(otu_table(CPSII_genus_RA)), cluster = sample_data(CPSII_genus_RA)$'disease_alc_status')
summary(isa_cpsII)

# Extract taxonomy table
plco_taxtable <- tax_table(plco_rare) |> as.data.frame() |> rownames_to_column(var="ASV")

cpsII_taxtable <- tax_table(cpsII_rare) |> as.data.frame() |> rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res <- isa_plco$sign |>
  rownames_to_column(var="ASV") |>
  left_join(plco_taxtable) |>
  filter(p.value < 0.05)
view(isa_res)

isa_cpsII_res <- isa_cpsII$sign |>
  rownames_to_column(var="ASV") |>
  left_join(cpsII_taxtable) |>
  filter(p.value < 0.05)
view(isa_cpsII_res)


# melt PLCO data and filter for significant genus
melted_plco <- psmelt(PLCO_genus_RA) |>
  filter(Genus %in% isa_res$Genus) |>
  filter(!str_detect(Genus, "^NA"))

# plot the graph facetted with Genus and relative abundance of each sample in each category
ISA_plco_plot <- ggplot(melted_plco, aes(x = disease_alc_status, y = Abundance)) +
  geom_boxplot(aes(fill = disease_alc_status), outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Genus, scales = "free_y") +
  labs(
    x = "Cancer and Alcohol Status",
    y = "Relative Abundance",
    fill = "Cancer and Alcohol Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) +
  scale_fill_manual(values = palette, labels = label)
ISA_plco_plot


# melt CPSII data and filter for significant genus
melted_cpsII <- psmelt(CPSII_genus_RA) |>
  filter(Genus %in% isa_cpsII_res$Genus) |>
  filter(!str_detect(Genus, "^NA"))

# plot the graph facetted with Genus and relative abundance of each sample in each category
ISA_cpsII_plot <- ggplot(melted_cpsII, aes(x = disease_alc_status, y = Abundance)) +
  geom_boxplot(aes(fill = disease_alc_status), outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Genus, scales = "free_y") +
  labs(
    x = "Cancer and Alcohol Status",
    y = "Relative Abundance",
    fill = "Cancer and Alcohol Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) +
  scale_fill_manual(values = palette, labels = label)
ISA_cpsII_plot

# save plots
ggsave("ISA_PLCO.png",
       ISA_plco_plot,
       height = 4, width = 6)
ggsave("ISA_CPSII.png",
       ISA_cpsII_plot,
       height = 4, width = 6)
