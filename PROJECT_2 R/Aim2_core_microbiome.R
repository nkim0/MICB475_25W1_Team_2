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
  ggplot2::scale_fill_distiller(direction = 1) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  labs(title = "Core Microbiome Venn Diagram (prevalence > 40% & abundance > 0.01%)")
core_microbiome_plot

# 1. Agglomerate at genus level
cpsII_genus_unrare <- tax_glom(CPSII, taxrank = "Genus", NArm = FALSE)

# transform data to relative abundance
cpsII_rel <- transform_sample_counts(cpsII_genus_unrare, fun = function(x) x/sum(x))

# subset dataset into 4 groups from testing data
test_cancer_yes <- subset_samples(cpsII_rel, disease_alc_status == "Case:yes")
test_cancer_no <- subset_samples(cpsII_rel, disease_alc_status == "Case:no")
test_control_yes <- subset_samples(cpsII_rel, disease_alc_status == "Control:yes")
test_control_no <- subset_samples(cpsII_rel, disease_alc_status == "Control:no")

#detect the core microbiome for each group
test_core_cancer_yes <- core_members(test_cancer_yes, detection = abund_threshold, prevalence = prev_threshold)
test_core_cancer_no <- core_members(test_cancer_no, detection = abund_threshold, prevalence = prev_threshold)
test_core_control_yes <- core_members(test_control_yes, detection = abund_threshold, prevalence = prev_threshold)
test_core_control_no <- core_members(test_control_no, detection = abund_threshold, prevalence = prev_threshold)

# plot the core microbiome
test_microbiome_plot <- ggVennDiagram(x=list(
  "Cancer, Yes" = test_core_cancer_yes,
  "Cancer, No" = test_core_cancer_no,
  "No Cancer, Yes" = test_core_control_yes,
  "No Cancer, No" = test_core_control_no
),label_alpha = 0, set_size=4) +
  ggplot2::scale_fill_distiller(direction = 1) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  labs(title = "Core Microbiome Venn Diagram (prevalence > 40% & abundance > 0.01%)")
test_microbiome_plot

# save plots
ggsave("microbiome_plco.png",
       core_microbiome_plot,
       height = 4, width = 6)
ggsave("microbiome_cpsII.png",
       test_microbiome_plot,
       height = 4, width = 6)
