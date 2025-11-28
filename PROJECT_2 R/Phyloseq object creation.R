library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(patchwork)
library(agricolae)
library(FSA)
library(rcompanion)


### Importing plco and cpsII metadata ###
meta_plco_dir <- "sd_plco/meta_export/metadata.txt"
meta_plco <- read_delim(meta_plco_dir, delim="\t")

meta_cpsII_dir <- "sd_cpsII/meta_export/metadata.txt"
meta_cpsII <- read_delim(meta_cpsII_dir, delim="\t")


### Importing plco and cpsII feature tables ###
otu_plco_dir <- "sd_plco/table_export/feature-table_plco.txt"
otu_plco <- read_delim(file = otu_plco_dir, delim="\t", skip=1)

otu_cpsII_dir <- "sd_cpsII/table_export/feature-table_cpsII.txt"
otu_cpsII <- read_delim(file = otu_cpsII_dir, delim="\t", skip=1)


### Importing plco and cpsII taxonomy tables ###
tax_plco_dir <- "sd_plco/taxa_export/taxonomy.tsv"
tax_plco <- read_delim(tax_plco_dir, delim="\t")

tax_cpsII_dir <- "sd_cpsII/taxa_export/taxonomy.tsv"
tax_cpsII <- read_delim(tax_cpsII_dir, delim="\t")


### Importing plco and cpsII phylogenetic trees ###
phylotree_plco_dir <- "sd_plco/tree_export/tree.nwk"
phylotree_plco <- read.tree(phylotree_plco_dir)

phylotree_cpsII_dir <- "sd_cpsII/tree_export/tree.nwk"
phylotree_cpsII <- read.tree(phylotree_cpsII_dir)


#### Cleaning Feature table and creation of phyloseq OTU table ####
otu_mat_plco <- as.matrix(otu_plco[,-1])
rownames(otu_mat_plco) <- otu_plco$`#OTU ID`
OTU_PLCO <- otu_table(otu_mat_plco, taxa_are_rows = TRUE) 

otu_mat_cpsII <- as.matrix(otu_cpsII[,-1])
rownames(otu_mat_cpsII) <- otu_cpsII$`#OTU ID`
OTU_CPSII <- otu_table(otu_mat_cpsII, taxa_are_rows = TRUE) 



#### Processing Metadata files ####
## filtering plco and cpsII metadata for alcohol and current smokers ##
# filtering for low risk alcohol consumption, females <= 26.9, males <= 53.8

meta_plco <- meta_plco |>
  filter(!is.na(host_ETHA_GRAMS_PER_DAY)) |>
  filter((hostgender == "female" & host_ETHA_GRAMS_PER_DAY <= 26.9) |
           (hostgender == "male" & host_ETHA_GRAMS_PER_DAY <= 53.8)) |>
  mutate(host_alc_category = case_when(
    host_ETHA_GRAMS_PER_DAY == 0 ~ "no",
    hostgender == "male" & host_ETHA_GRAMS_PER_DAY > 0 & host_ETHA_GRAMS_PER_DAY <= 53.8 ~ "yes",
    hostgender == "female" & host_ETHA_GRAMS_PER_DAY > 0 & host_ETHA_GRAMS_PER_DAY <= 26.9 ~ "yes")) |>
  filter(host_smoke == "Current") |>
  mutate(disease_alc_status = paste(disease_status, host_alc_category, sep = ":"))

meta_cpsII <- meta_cpsII |>
  filter(!is.na(host_ETHA_GRAMS_PER_DAY)) |>
  filter((hostgender == "female" & host_ETHA_GRAMS_PER_DAY <= 26.9) |
           (hostgender == "male" & host_ETHA_GRAMS_PER_DAY <= 53.8)) |>
  mutate(host_alc_category = case_when(
    host_ETHA_GRAMS_PER_DAY == 0 ~ "no",
  hostgender == "male" & host_ETHA_GRAMS_PER_DAY > 0 & host_ETHA_GRAMS_PER_DAY <= 53.8 ~ "yes",
    hostgender == "female" & host_ETHA_GRAMS_PER_DAY > 0 & host_ETHA_GRAMS_PER_DAY <= 26.9 ~ "yes")) |>
  filter(host_smoke == "Current") |>
  mutate(disease_alc_status = paste(disease_status, host_alc_category, sep = ":"))

## Formatting and creating metadata phyloseq file ##
samp_df_plco <- as.data.frame(meta_plco)
rownames(samp_df_plco)<- meta_plco$'sample-id'
SAMP_PLCO <- sample_data(samp_df_plco)


samp_df_cpsII <- as.data.frame(meta_cpsII)
rownames(samp_df_cpsII)<- meta_cpsII$'sample-id'
SAMP_CPSII <- sample_data(samp_df_cpsII)



### Formatting and creation of phyloseq taxonomy file ###
tax_mat_plco <- tax_plco %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat_plco <- tax_mat_plco[,-1]
rownames(tax_mat_plco) <- tax_plco$`Feature ID`
TAX_PLCO <- tax_table(tax_mat_plco)


tax_mat_cpsII <- tax_cpsII %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat_cpsII <- tax_mat_cpsII[,-1]
rownames(tax_mat_cpsII) <- tax_cpsII$`Feature ID`
TAX_CPSII <- tax_table(tax_mat_cpsII)


### Creation of Phylogenetic tree phyloseq object ###
phylotree_plco_dir <- "sd_plco/tree_export/tree.nwk"
phylotree_plco <- read.tree(phylotree_plco_dir)

phylotree_cpsII_dir <- "sd_cpsII/tree_export/tree.nwk"
phylotree_cpsII <- read.tree(phylotree_cpsII_dir)


### Creation of plco and cpsII phyloseq object ###
PLCO <- phyloseq(OTU_PLCO, SAMP_PLCO, TAX_PLCO, phylotree_plco)
CPSII <- phyloseq(OTU_CPSII, SAMP_CPSII, TAX_CPSII, phylotree_cpsII)


### Rarefy data ###
PLCO_filt <- subset_taxa(PLCO,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
PLCO_nolow <- filter_taxa(PLCO_filt, function(x) sum(x)>5, prune = TRUE)
rarecurve(t(as.data.frame(otu_table(PLCO_nolow))), cex=0.1)
plco_rare <- rarefy_even_depth(PLCO_nolow, rngseed = 1, sample.size = 4138)

CPSII_filt <- subset_taxa(CPSII,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
CPSII_nolow <- filter_taxa(CPSII_filt, function(x) sum(x)>5, prune = TRUE)
rarecurve(t(as.data.frame(otu_table(CPSII_nolow))), cex=0.1)
cpsII_rare <- rarefy_even_depth(CPSII_nolow, rngseed = 1, sample.size = 2664)

### Creating labels and palette list ###
palette <- c("#FC8D62FF", "#8DA0CBFF", "#E78AC3FF", "#A6D854FF")
label <- c("Case:yes" = "Cancer, Yes", 
           "Case:no" = "Cancer, No", 
           "Control:yes" = "No Cancer, Yes", 
           "Control:no" = "No Cancer, No")

save(plco_rare, cpsII_rare, CPSII, PLCO, palette, label, file = "phylseq_files.RData")
