library(tidyverse)
library(phyloseq)
library(ggpicrust2)
devtools::install_github("cafferychen777/ggpicrust2")
source("ggpicrust2_errorbar_function_fixed.R")

# load objects
load(file = "phylseq_files.RData")

# Load data
meta_plco <- sample_data(PLCO_genus)

ko <- read.delim("PiCrust2_Data/pred_metagenome_unstrat.tsv", row.names = 1)

# Convert sample_data to data frame before filtering to avoid vec_slice() error
#3 combinations with Control:no, (1) Control:no and Case:yes, (2) Control:no and Case:no, (3) Control:no and Control:yes
meta_plco_filt_1 <- subset_samples(meta_plco, disease_alc_status %in% c("Control:no", "Case:yes"))
ko_filt_1 <- ko |> select(all_of(meta_plco_filt_1$sample.id))

meta_plco_filt_2 <- subset_samples(meta_plco, disease_alc_status %in% c("Control:no", "Case:no"))
ko_filt_2 <- ko |> select(all_of(meta_plco_filt_2$sample.id))

meta_plco_filt_3 <- subset_samples(meta_plco, disease_alc_status %in% c("Control:no", "Control:yes"))
ko_filt_3 <- ko |> select(all_of(meta_plco_filt_3$sample.id))

# perform pathway differential abundance analysis (DAA) using LinDA method
daa_results_df_1 <- pathway_daa(abundance = ko_filt_1,
                              metadata = meta_plco_filt_1,
                              group = "disease_alc_status",
                              daa_method = "LinDA",
                              select = NULL, reference = "Control:no")

daa_results_df_2 <- pathway_daa(abundance = ko_filt_2,
                                metadata = meta_plco_filt_2,
                                group = "disease_alc_status",
                                daa_method = "LinDA",
                                select = NULL, reference = "Control:no")

daa_results_df_3 <- pathway_daa(abundance = ko_filt_3,
                                metadata = meta_plco_filt_3,
                                group = "disease_alc_status",
                                daa_method = "LinDA",
                                select = NULL, reference = "Control:no")

# annotate pathway results using KO to KEGG conversion
data_annotated_results_df_1 <- pathway_annotation(pathway = "KO",
                                                daa_results_df = daa_results_df_1,
                                                ko_to_kegg = TRUE)

data_annotated_results_df_2 <- pathway_annotation(pathway = "KO",
                                                  daa_results_df = daa_results_df_2,
                                                  ko_to_kegg = TRUE)

data_annotated_results_df_3 <- pathway_annotation(pathway = "KO",
                                                  daa_results_df = daa_results_df_3,
                                                  ko_to_kegg = TRUE)

# Creating Graph
peb_1 <- pathway_errorbar(abundance = ko_filt_1,
                        daa_results_df = data_annotated_results_df_1,
                        Group = meta_plco_filt_1$disease_alc_status,
                        p_values_threshold = 0.05,
                        order = "pathway_class",
                        p_value_bar = TRUE,
                        x_lab = "pathway_name")

peb_2 <- pathway_errorbar(abundance = ko_filt_2,
                          daa_results_df = data_annotated_results_df_2,
                          Group = meta_plco_filt_2$disease_alc_status,
                          p_values_threshold = 0.05,
                          order = "pathway_class",
                          p_value_bar = TRUE,
                          x_lab = "pathway_name")

peb_3 <- pathway_errorbar(abundance = ko_filt_3,
                          daa_results_df = data_annotated_results_df_3,
                          Group = meta_plco_filt_3$disease_alc_status,
                          p_values_threshold = 0.05,
                          order = "pathway_class",
                          p_value_bar = TRUE,
                          x_lab = "pathway_name")
peb_1
peb_2
peb_3


combined_data_annotated_results_df <- bind_rows(data_annotated_results_df_1, data_annotated_results_df_3)


ggsave("func_analysis.png",
       plot = peb_1, width = 18, height = 10)
ggsave("func_analysis.png",
       plot = peb_2, width = 18, height = 10)
ggsave("func_analysis.png",
       plot = peb_3, width = 18, height = 10)

