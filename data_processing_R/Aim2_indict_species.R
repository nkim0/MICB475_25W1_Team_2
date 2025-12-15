library(phyloseq)
library(indicspecies)
library(tidyverse)
library(kableExtra)

# load objects
load(file = "phylseq_files.RData")

### AIM 2.2: Indicator Species Analysis ###
# Agglomerate at genus level and convert to relative abundance
PLCO_genus <- tax_glom(PLCO, taxrank = "Genus", NArm = FALSE)
PLCO_genus_RA <- transform_sample_counts(PLCO_genus, fun = function(x) x / sum(x))

# Run indicator species analysis using the donor_status variable and look at results
isa_plco <- multipatt(t(otu_table(PLCO_genus_RA)), cluster = sample_data(PLCO_genus_RA)$"disease_alc_status")
summary(isa_plco)

# Extract taxonomy table
plco_taxtable <- tax_table(PLCO) |>
  as.data.frame() |>
  rownames_to_column(var = "ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
isa_res <- isa_plco$sign |>
  rownames_to_column(var = "ASV") |>
  left_join(plco_taxtable) |>
  filter(p.value < 0.05)
view(isa_res)

# Format ISA results table
isa_plco_table <- isa_res |>
  filter(!is.na(Genus)) |>
  mutate(Genus = str_remove(Genus, "^g__")) |>
  mutate(Genus = str_replace_all(Genus, "_", " ")) |>
  select(Genus, 
    "Cancer, No"     = `s.Case:no`,
    "Cancer, Yes"    = `s.Case:yes`,
    "No Cancer, No"  = `s.Control:no`, 
    "No Cancer, Yes" = `s.Control:yes`, 
    p.value, stat)|>
  mutate(stat = round(stat, 3)) |>
  arrange(p.value) |>
  rename("Indicator Value" = stat,
         "P-value" = p.value)

# Create formatted table for PLCO
kable_plco <- isa_plco_table |>
  mutate(across(c("Cancer, No", "Cancer, Yes", "No Cancer, No", "No Cancer, Yes"), 
                ~ ifelse(. == 1, "✔", ""))) |>
  kbl(align = c("l", "c", "c", "c", "c", "c", "c", "c")) |>
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE,
                font_size = 12) |>
  add_header_above(c(" " = 1, "Cancer Status and Alcohol Consumption" = 4, " " = 2)) |>
  column_spec(1, bold = TRUE, italic = TRUE) |>
  row_spec(0, bold = TRUE, background = "#E8E8E8")

kable_plco

# Save as HTML
save_kable(kable_plco, file = "PLCO_ISA_table.html")
