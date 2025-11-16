library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# load objects
load(file = "phylseq_files.RData")

# get the name of each analysis group
print(levels(as.factor(sample_data(plco_rare)$disease_alc_status)))

### AIM 2.1: Core Microbiome ###
# set thersholds
prev_threshold <- 0.3
abund_threshold <- 0.001

# transform data to relative abundance
plco_rel <- transform_sample_counts(PLCO, fun = function(x) x/sum(x))

# subset dataset into 4 groups from training data
cancer_low_risk <- subset_samples(plco_rel, disease_alc_status == "Case:low_risk")
cancer_no_risk <- subset_samples(plco_rel, disease_alc_status == "Case:no_risk")
control_low_risk <- subset_samples(plco_rel, disease_alc_status == "Control:low_risk")
control_no_risk <- subset_samples(plco_rel, disease_alc_status == "Control:no_risk")

#detect the core microbiome for each group
core_cancer_low_risk <- core_members(cancer_low_risk, detection = abund_threshold, prevalence = prev_threshold)
core_cancer_no_risk <- core_members(cancer_no_risk, detection = abund_threshold, prevalence = prev_threshold)
core_control_low_risk <- core_members(control_low_risk, detection = abund_threshold, prevalence = prev_threshold)
core_control_no_risk <- core_members(control_no_risk, detection = abund_threshold, prevalence = prev_threshold)

# plot the core microbiome
core_microbiome_plot <- ggVennDiagram(x=list(
  "Cancer, Low Risk" = core_cancer_low_risk,
  "Cancer, No Risk" = core_cancer_no_risk,
  "No Cancer, Low Risk" = core_control_low_risk,
  "No Cancer, No Risk" = core_control_no_risk
),label_alpha = 0, set_size=4) +
  ggplot2::scale_fill_distiller(direction = 1) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0.2))
core_microbiome_plot

# Check number of samples for all 4 subsetted groups
print(paste("NoRisk_NoCancer n=", nsamples(cancer_low_risk)))
print(paste("NoRisk_Cancer n=", nsamples(cancer_no_risk)))
print(paste("HighRisk_NoCancer n=", nsamples(control_low_risk)))
print(paste("HighRisk_Cancer n=", nsamples(control_no_risk)))


# transform data to relative abundance for test data
cpsII_rel <- transform_sample_counts(CPSII, fun = function(x) x/sum(x))

# subset dataset into 4 groups from testing data
test_cancer_low_risk <- subset_samples(cpsII_rel, disease_alc_status == "Case:low_risk")
test_cancer_no_risk <- subset_samples(cpsII_rel, disease_alc_status == "Case:no_risk")
test_control_low_risk <- subset_samples(cpsII_rel, disease_alc_status == "Control:low_risk")
test_control_no_risk <- subset_samples(cpsII_rel, disease_alc_status == "Control:no_risk")

#detect the core microbiome for each group
test_core_cancer_low_risk <- core_members(test_cancer_low_risk, detection = abund_threshold, prevalence = prev_threshold)
test_core_cancer_no_risk <- core_members(test_cancer_no_risk, detection = abund_threshold, prevalence = prev_threshold)
test_core_control_low_risk <- core_members(test_control_low_risk, detection = abund_threshold, prevalence = prev_threshold)
test_core_control_no_risk <- core_members(test_control_no_risk, detection = abund_threshold, prevalence = prev_threshold)

# plot the core microbiome
test_microbiome_plot <- ggVennDiagram(x=list(
  "Cancer, Low Risk" = test_core_cancer_low_risk,
  "Cancer, No Risk" = test_core_cancer_no_risk,
  "No Cancer, Low Risk" = test_core_control_low_risk,
  "No Cancer, No Risk" = test_core_control_no_risk
),label_alpha = 0, set_size=4) +
  ggplot2::scale_fill_distiller(direction = 1) +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0.2))
test_microbiome_plot

# save plots
ggsave("microbiome_plco.png",
       core_microbiome_plot,
       height = 4, width = 6)
ggsave("microbiome_cpsII.png",
       test_microbiome_plot,
       height = 4, width = 6)
