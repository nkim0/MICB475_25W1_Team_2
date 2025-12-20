library(tidyverse)
library(ggplot2)
library(phyloseq)
library(DESeq2)

# load objects
load(file = "phylseq_files.RData")

# get the name of each analysis group
print(levels(as.factor(sample_data(plco_rare)$disease_alc_status)))

### AIM 3: Differential Abundance ###
# transform PLCO to have +1 count and change to DESeq object
PLCO_plus1 <- transform_sample_counts(PLCO, function(x) x+1)
PLCO_deseq <- phyloseq_to_deseq2(PLCO_plus1, ~disease_alc_status)
DESEQ_PLCO <- DESeq(PLCO_deseq)

# get the DESeq result, while Control:no_risk is the reference
#3 combinations with Control:no, (1) Control:no and Case:yes, (2) Control:no and Case:no, (3) Control:no and Control:yes
DESeq_PLCO_res_1 <- results(DESEQ_PLCO, tidy=TRUE,
                     contrast = c("disease_alc_status", "Case:yes", "Control:no"))

DESeq_PLCO_res_2 <- results(DESEQ_PLCO, tidy=TRUE,
                            contrast = c("disease_alc_status", "Case:no", "Control:no"))

DESeq_PLCO_res_3 <- results(DESEQ_PLCO, tidy=TRUE,
                            contrast = c("disease_alc_status", "Control:yes", "Control:no"))
view(DESeq_PLCO_res_1)
view(DESeq_PLCO_res_2)
view(DESeq_PLCO_res_3)


# get table of significant ASVs
PLCO_sigASVs_1 <- DESeq_PLCO_res_1 |>
  filter(padj<0.01 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

PLCO_sigASVs_2 <- DESeq_PLCO_res_2 |>
  filter(padj<0.01 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

PLCO_sigASVs_3 <- DESeq_PLCO_res_3 |>
  filter(padj<0.01 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

# Get a vector of significant ASV names
PLCO_sigASVs_vec_1 <- PLCO_sigASVs_1 |>
  pull(ASV)

PLCO_sigASVs_vec_2 <- PLCO_sigASVs_2 |>
  pull(ASV)

PLCO_sigASVs_vec_3 <- PLCO_sigASVs_3 |>
  pull(ASV)

# Prune phyloseq obj for significant ASVs
PLCO_filt_1 <- prune_taxa(PLCO_sigASVs_vec_1, PLCO_plus1)

PLCO_filt_2 <- prune_taxa(PLCO_sigASVs_vec_2, PLCO_plus1)

PLCO_filt_3 <- prune_taxa(PLCO_sigASVs_vec_3, PLCO_plus1)

# add taxonomy up to Genus level on to DESeq results tables
PLCO_merged_results_1 <- tax_table(PLCO_plus1) |> 
  as.data.frame()|>
  rownames_to_column(var="ASV") %>%
  right_join(PLCO_sigASVs_1) %>%
  arrange(padj) |>
  slice_head(n = 25) |>
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) |>
  filter(!str_detect(Genus, "^NA"))

PLCO_merged_results_2 <- tax_table(PLCO_plus1) |> 
  as.data.frame()|>
  rownames_to_column(var="ASV") %>%
  right_join(PLCO_sigASVs_2) %>%
  arrange(padj) |>
  slice_head(n = 25) |>
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) |>
  filter(!str_detect(Genus, "^NA"))

PLCO_merged_results_3 <- tax_table(PLCO_plus1) |> 
  as.data.frame()|>
  rownames_to_column(var="ASV") %>%
  right_join(PLCO_sigASVs_3) %>%
  arrange(padj) |>
  slice_head(n = 25) |>
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) |>
  filter(!str_detect(Genus, "^NA"))

# Make DESeq volcano plots
vol_plot_1 <- DESeq_PLCO_res_1 |>
  mutate(Significance = case_when(padj < 0.01 & abs(log2FoldChange) > 2 ~ "Significant", TRUE ~ "Not Significant")) |>
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = Significance)) +
  theme_classic()

vol_plot_2 <- DESeq_PLCO_res_2 |>
  mutate(Significance = case_when(padj < 0.01 & abs(log2FoldChange) > 2 ~ "Significant", TRUE ~ "Not Significant")) |>
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = Significance)) +
  theme_classic()

vol_plot_3 <- DESeq_PLCO_res_3 |>
  mutate(Significance = case_when(padj < 0.01 & abs(log2FoldChange) > 2 ~ "Significant", TRUE ~ "Not Significant")) |>
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = Significance)) +
  theme_classic()

# Make DESeq bar plot
gg_PLCO_barplot_1 <- PLCO_merged_results_1 |>
  ggplot() +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  theme_classic() +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE,)) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=10)) +
  labs(y="Log2 Fold Change")

gg_PLCO_barplot_2 <- ggplot(PLCO_merged_results_2) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  theme_classic() +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=10)) +
  labs(y="Log2 Fold Change")

gg_PLCO_barplot_3 <- ggplot(PLCO_merged_results_3) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  theme_classic() +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=10)) +
  labs(y="Log2 Fold Change")



vol_plot_1
vol_plot_2
vol_plot_3

gg_PLCO_barplot_1
gg_PLCO_barplot_2
gg_PLCO_barplot_3

# save plots
  ggsave("DESeq_PLCO_ctno,csyes.png",
         gg_PLCO_barplot_1,
         height = 5, width = 7.5)
  ggsave("DESeq_PLCO_ctno,csno.png",
         gg_PLCO_barplot_2,
         height = 5, width = 7.5)
  ggsave("DESeq_PLCO_ctno,ctyes.png",
         gg_PLCO_barplot_3,
         height = 5, width = 7.5)
  
  ggsave("DESeq_PLCO_vol_ctno,csyes.png",
         vol_plot_1,
         height = 4, width = 6)
  ggsave("DESeq_PLCO_vol_ctno,csno.png",
         vol_plot_2,
         height = 4, width = 6)
  ggsave("DESeq_PLCO_vol_ctno,ctyes.png",
         vol_plot_3,
         height = 4, width = 6)

