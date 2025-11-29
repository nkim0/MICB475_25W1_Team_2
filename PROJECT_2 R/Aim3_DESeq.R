library(tidyverse)
library(ggplot2)
library(phyloseq)
library(DESeq2)

# load objects
load(file = "phylseq_files.RData")

# get the name of each analysis group
print(levels(as.factor(sample_data(cpsII_rare)$disease_alc_status)))

### AIM 3: Differential Abundance ###
# transform PLCO to have +1 count and change to DESeq object
PLCO_plus1 <- transform_sample_counts(PLCO, function(x) x+1)
PLCO_deseq <- phyloseq_to_deseq2(PLCO_plus1, ~disease_alc_status)
DESEQ_PLCO <- DESeq(PLCO_deseq)

# get the DESeq result, while Control:no_risk is the reference and Case:low_risk is treatment group
DESeq_PLCO_res <- results(DESEQ_PLCO, tidy=TRUE,
                     contrast = c("disease_alc_status", "Case:yes", "Control:no"))
view(DESeq_PLCO_res)


# get table of significant ASVs
PLCO_sigASVs <- DESeq_PLCO_res |>
  filter(padj<0.01 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

# Get a vector of significant ASV names
PLCO_sigASVs_vec <- PLCO_sigASVs |>
  pull(ASV)

# Prune phyloseq obj for significant ASVs
PLCO_filt <- prune_taxa(PLCO_sigASVs_vec, PLCO_plus1)

# add taxonomy up to Genus level on to DESeq results tables
PLCO_merged_results <- tax_table(PLCO_plus1) |> 
  as.data.frame()|>
  rownames_to_column(var="ASV") %>%
  right_join(PLCO_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) |>
  filter(!str_detect(Genus, "^NA"))

# Make DESeq bar plot
gg_PLCO_barplot <- ggplot(PLCO_merged_results) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  theme_classic() +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=7)) +
  labs(y="Log2 Fold Change")
gg_PLCO_barplot



# Perform Deseq on CPSII
# transform CPSII to have +1 count and change to DESeq object
CPSII_plus1 <- transform_sample_counts(CPSII, function(x) x+1)
CPSII_deseq <- phyloseq_to_deseq2(CPSII_plus1, ~disease_alc_status)
DESEQ_CPSII <- DESeq(CPSII_deseq)

# get the DESeq result, while Control:no_risk is the reference and Case:low_risk is treatment group
DESeq_CPSII_res <- results(DESEQ_CPSII, tidy=TRUE,
               contrast = c("disease_alc_status", "Case:yes", "Control:no"))
view(DESeq_CPSII_res)

# get table of significant ASVs
CPSII_sigASVs <- DESeq_CPSII_res |>
  filter(padj<0.01 & abs(log2FoldChange)>2) |>
  dplyr::rename(ASV=row)

# Get a vector of significant ASV names
CPSII_sigASVs_vec <- CPSII_sigASVs |>
  pull(ASV)

# Prune phyloseq obj for significant ASVs
CPSII_filt <- prune_taxa(CPSII_sigASVs_vec, CPSII_plus1)

# add taxonomy up to Genus level on to DESeq results tables
CPSII_merged_results <- tax_table(CPSII_plus1) |> 
  as.data.frame()|>
  rownames_to_column(var="ASV") %>%
  right_join(CPSII_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) |>
  filter(!str_detect(Genus, "^NA"))

# Make DESeq bar plot
gg_CPSII_barplot <- ggplot(CPSII_merged_results) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  theme_classic() +
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=7)) +
  labs(y="Log2 Fold Change")
gg_CPSII_barplot

# save plots
ggsave("DESeq_PLCO.png",
       gg_PLCO_barplot,
       height = 4, width = 6)
ggsave("DESeq_CPSII.png",
       gg_CPSII_barplot,
       height = 4, width = 6)
