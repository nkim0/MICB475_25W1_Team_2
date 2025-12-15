library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# load objects
load(file = "phylseq_files.RData")

# get the name of each analysis group
print(levels(as.factor(sample_data(plco_rare)$disease_alc_status)))

### AIM 2.1: Core Microbiome ###
# set thersholds
prev_threshold <- 0.4
abund_threshold <- 0.0001

# 1. Agglomerate at genus level
plco_genus_unrare <- tax_glom(PLCO, taxrank = "Genus", NArm = FALSE)

# transform data to relative abundance
plco_rel <- transform_sample_counts(plco_genus_unrare, fun = function(x) x/sum(x))

# subset dataset into 4 groups from training data
cancer_yes_alc <- subset_samples(plco_rel, disease_alc_status == "Case:yes")
cancer_no_alc <- subset_samples(plco_rel, disease_alc_status == "Case:no")
control_yes_alc <- subset_samples(plco_rel, disease_alc_status == "Control:yes")
control_no_alc <- subset_samples(plco_rel, disease_alc_status == "Control:no")

#detect the core microbiome for each group
core_cancer_yes_alc <- core_members(cancer_yes_alc, detection = abund_threshold, prevalence = prev_threshold)
core_cancer_no_alc <- core_members(cancer_no_alc, detection = abund_threshold, prevalence = prev_threshold)
core_control_yes_alc <- core_members(control_yes_alc, detection = abund_threshold, prevalence = prev_threshold)
core_control_no_alc <- core_members(cancer_no_alc, detection = abund_threshold, prevalence = prev_threshold)

# plot the core microbiome
core_microbiome_plot <- ggVennDiagram(x=list(
  "Cancer, Yes" = core_cancer_yes_alc,
  "Cancer, No" = core_cancer_no_alc,
  "No Cancer, Yes" = core_control_yes_alc,
  "No Cancer, No" = core_control_no_alc
),label_alpha = 0, set_size=4) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  labs(title = "prevalence > 40% & abundance > 0.01%") +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
core_microbiome_plot

# save plots
ggsave("microbiome_plco.png",
       core_microbiome_plot,
       height = 5, width = 7)

