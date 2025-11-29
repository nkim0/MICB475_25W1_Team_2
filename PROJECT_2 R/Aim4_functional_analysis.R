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
meta_plco_filt <- subset_samples(meta_plco, disease_alc_status %in% c("Case:yes", "Control:no"))
ko_filt <- ko |> select(all_of(meta_plco_filt$sample.id))

# perform pathway differential abundance analysis (DAA) using LinDA method
daa_results_df <- pathway_daa(abundance = ko_filt,
                              metadata = meta_plco_filt,
                              group = "disease_alc_status",
                              daa_method = "LinDA",
                              select = NULL, reference = NULL)

# annotate pathway results using KO to KEGG conversion
data_annotated_results_df <- pathway_annotation(pathway = "KO",
                                                daa_results_df = daa_results_df,
                                                ko_to_kegg = TRUE)
# saveRDS(daa_annotated_results_df, "annotated_pathways.rds")
# 
# data_annotated_results_df <- readRDS("annotated_pathways.rds")


# find right p-value
# data_annotated_results_df |> 
#   filter(p_adjust<0.05, abs(log2FoldChange) > 2) |> 
#   nrow()


peb <- pathway_errorbar_fixed(abundance = ko_filt,
                        daa_results_df = data_annotated_results_df,
                        Group = meta_plco_filt$disease_alc_status,
                        wrap_label = T, wraplength = 60,
                        fc_cutoff = 2, order_by_log = F,
                        p_values_threshold = 0.05,
                        order = "pathway_class",
                        ko_to_kegg = T,
                        p_value_bar = F,
                        x_lab = "pathway_name")
ggsave("func_analysis.png",
       plot = peb, width = 18, height = 10)
