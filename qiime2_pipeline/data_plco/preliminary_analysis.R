library(tidyverse)
library(dplyr)
library(ggplot2)

# import 2 metadata and merge them
metafp_1 <- "NCI_PLCO_cohort.csv"
meta_1 <- read.csv(metafp_1)
metafp_2 <- "ACS CPS-II cohort.csv"
meta_2 <- read.csv(metafp_2)

meta <- rbind(meta_1, meta_2)

# filter out NA and high alcohol consumption based on gender
# classify no_risk or low_risk based on alcohol consumption based on gender (53.8 or 26.9 for males and females)
meta_df <- meta |>
  filter(!is.na(host_ETHA_GRAMS_PER_DAY)) |>
  filter((hostgender == "female" & host_ETHA_GRAMS_PER_DAY <= 26.9) |
         (hostgender == "male" & host_ETHA_GRAMS_PER_DAY <= 53.8)) |>
  mutate(host_alc_category = case_when(
    host_ETHA_GRAMS_PER_DAY == 0 ~ "no_risk",
    hostgender == "male" & host_ETHA_GRAMS_PER_DAY > 0 & host_ETHA_GRAMS_PER_DAY <= 53.8 ~ "low_risk",
    hostgender == "female" & host_ETHA_GRAMS_PER_DAY > 0 & host_ETHA_GRAMS_PER_DAY <= 26.9 ~ "low_risk"))

# separate meta data for control and cancer patients
meta_ctrl <- meta_df |>
  filter(disease_status == "Control")

meta_cancer <- meta_df |>
  filter(disease_status == "Case")

# create matrix of alcohol and smoking status
ctrl_matrix <- table(meta_ctrl$host_smoke, meta_ctrl$host_alc_category)
cancer_matrix <- table(meta_cancer$host_smoke, meta_cancer$host_alc_category)

# create smoke_alc_status variable, and find number of samples and list corresponding
# to each category
ctrl_samples <- meta_ctrl |>
  mutate(smoke_alc_status = paste(host_smoke, host_alc_category, sep = ":")) |>
  group_by(smoke_alc_status) |>
  summarise(
    n = n(),
    samples = list(BioSample)
  )

cancer_samples <- meta_cancer |>
  mutate(smoke_alc_status = paste(host_smoke, host_alc_category, sep = ":")) |>
  group_by(smoke_alc_status) |>
  summarise(
    n = n(),
    samples = list(BioSample)
  )


  


