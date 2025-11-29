library(phyloseq)
library(vegan)
library(microbiome)
library(RColorBrewer)

load(file = "phylseq_files.RData")


### Formating data and Conducting a weighted unifrac analysis on plco and cpsII data ###
plco_rel <- microbiome::transform(PLCO_genus, "compositional")
ord_plco <- ordinate(plco_rel, "PCoA", "unifrac", weighted = T)
weighted_unifraq_plco <- plot_ordination(plco_rel,
                                    ord_plco, color = "disease_alc_status") +
  theme_classic() +
  labs(color = "Cancer and Alcohol Consumption") +
  scale_color_manual(values = palette, labels = label)
weighted_unifraq_plco


cpsII_rel <- microbiome::transform(CPSII_genus, "compositional")
ord_cpsII <- ordinate(cpsII_rel, "PCoA", "unifrac", weighted = T)
weighted_unifraq_cpsII <- plot_ordination(cpsII_rel,
                                          ord_cpsII, color = "disease_alc_status") +
  theme_classic() +
  labs(color = "Cancer and Alcohol Consumption") +
  scale_color_manual(values = palette, labels = label)
weighted_unifraq_cpsII

### Conducting PERMANOVA test on weighted unifrac data ###
dist_plco <- phyloseq::distance(plco_rel, method = "unifrac", weighted = TRUE)
meta_plco <- data.frame(sample_data(plco_rel))
permanova_plco <- adonis2(dist_plco ~ disease_alc_status, 
                          data = meta_plco, 
                          permutations = 999)
print(permanova_plco)

dist_cpsII <- phyloseq::distance(cpsII_rel, method = "unifrac", weighted = TRUE)
meta_cpsII <- data.frame(sample_data(cpsII_rel))
permanova_cpsII <- adonis2(dist_cpsII ~ disease_alc_status, 
                          data = meta_cpsII, 
                          permutations = 999)
print(permanova_cpsII)

ggsave("wu_plco.png"
       , weighted_unifraq_plco
       , height = 4, width = 5)

ggsave("wu_cpsII.png"
       , weighted_unifraq_cpsII
       , height = 4, width = 5)
