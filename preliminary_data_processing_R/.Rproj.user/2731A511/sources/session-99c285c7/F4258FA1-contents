library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(patchwork)
library(agricolae)
library(FSA)
library(rcompanion)


### Importing plco metadata ###
meta_plco_dir <- "sd_plco/meta_export/metadata.txt"
meta_plco <- read_delim(meta_plco_dir, delim="\t")

### Importing plco feature tables ###
otu_plco_dir <- "sd_plco/table_export/feature-table_plco.txt"
otu_plco <- read_delim(file = otu_plco_dir, delim="\t", skip=1)

### Importing plco taxonomy tables ###
tax_plco_dir <- "sd_plco/taxa_export/taxonomy.tsv"
tax_plco <- read_delim(tax_plco_dir, delim="\t")

### Importing plco phylogenetic trees ###
phylotree_plco_dir <- "sd_plco/tree_export/tree.nwk"
phylotree_plco <- read.tree(phylotree_plco_dir)

#### Cleaning Feature table and creation of phyloseq OTU table ####
otu_mat_plco <- as.matrix(otu_plco[,-1])
rownames(otu_mat_plco) <- otu_plco$`#OTU ID`
OTU_PLCO <- otu_table(otu_mat_plco, taxa_are_rows = TRUE) 

#### Processing Metadata files ####
## filtering plco metadata for alcohol and current smokers ##
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

## Formatting and creating metadata phyloseq file ##
samp_df_plco <- as.data.frame(meta_plco)
rownames(samp_df_plco)<- meta_plco$'sample-id'
SAMP_PLCO <- sample_data(samp_df_plco)

### Formatting and creation of phyloseq taxonomy file ###
tax_mat_plco <- tax_plco %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat_plco <- tax_mat_plco[,-1]
rownames(tax_mat_plco) <- tax_plco$`Feature ID`
TAX_PLCO <- tax_table(tax_mat_plco)

### Creation of Phylogenetic tree phyloseq object ###
phylotree_plco_dir <- "sd_plco/tree_export/tree.nwk"
phylotree_plco <- read.tree(phylotree_plco_dir)

### Creation of plco phyloseq object ###
PLCO <- phyloseq(OTU_PLCO, SAMP_PLCO, TAX_PLCO, phylotree_plco)

### Rarefy data ###
PLCO_filt <- subset_taxa(PLCO,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
PLCO_nolow <- filter_taxa(PLCO_filt, function(x) sum(x)>5, prune = TRUE)
rarecurve(t(as.data.frame(otu_table(PLCO_nolow))), cex=0.1)
plco_rare <- rarefy_even_depth(PLCO_nolow, rngseed = 1, sample.size = 4138)

# Agglomerate at genus level and convert to relative abundance
PLCO_genus <- tax_glom(plco_rare, taxrank = "Genus", NArm = FALSE)
PLCO_genus_RA <- transform_sample_counts(PLCO_genus, fun=function(x) x/sum(x))

### Creating labels and palette list ###
palette <- c("#FC8D62FF", "#8DA0CBFF", "#E78AC3FF", "#A6D854FF")
label <- c("Case:yes" = "Cancer, Yes", 
           "Case:no" = "Cancer, No", 
           "Control:yes" = "No Cancer, Yes", 
           "Control:no" = "No Cancer, No")

save(plco_rare, PLCO, palette, label,PLCO_genus, file = "phylseq_files.RData")