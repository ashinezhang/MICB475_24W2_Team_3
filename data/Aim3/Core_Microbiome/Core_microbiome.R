# Load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

colombia_meta <- read_tsv("metadata.tsv")

# Convert to relative abundance
colombia_RA <- transform_sample_counts(colombia_final, fun=function(x) x/sum(x))

# Filter dataset by sex and diabetic status
d_f <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Diabetic_female")
pd_f <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Pre-diabetic_female")
nd_f <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Non-diabetic_female")

d_m <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Diabetic_male")
pd_m <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Pre-diabetic_male")
nd_m <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Non-diabetic_male")

# Set prevalence and abundance threshold
d_f_ASVs <- core_members(d_f, detection=0.001, prevalence = 0.8)
pd_f_ASVs <- core_members(pd_f, detection=0.001, prevalence = 0.8)
nd_f_ASVs <- core_members(nd_f, detection=0.001, prevalence = 0.8)

d_m_ASVs <- core_members(d_m, detection=0.001, prevalence = 0.8)
pd_m_ASVs <- core_members(pd_m, detection=0.001, prevalence = 0.8)
nd_m_ASVs <- core_members(nd_m, detection=0.001, prevalence = 0.8)


# Venn-diagram diabetic female vs. diabetic male
d_m_f_core_microbiome <- ggVennDiagram(x=list(
  "Diabetic Male" = d_m_ASVs, 
  "Diabetic Female" = d_f_ASVs)) 
d_m_f_core_microbiome

# Venn-diagram pre-diabetic female vs. pre-diabetic male
pd_m_f_core_microbiome <- ggVennDiagram(x=list(
  "Pre-diabetic Male" = pd_m_ASVs, 
  "Pre-diabetic Female" = pd_f_ASVs)) 
pd_m_f_core_microbiome

# Venn-diagram non-diabetic female vs. non-diabetic male
nd_m_f_core_microbiome <- ggVennDiagram(x=list(
  "Non-diabetic Male" = nd_m_ASVs, 
  "Non-diabetic Female" = nd_f_ASVs)) 
nd_m_f_core_microbiome

ggsave("d_m_f_core_microbiome.png", d_m_f_core_microbiome, height = 12, width = 15)
ggsave("pd_m_f_core_microbiome.png", pd_m_f_core_microbiome, height = 12, width = 15)
ggsave("nd_m_f_core_microbiome.png", nd_m_f_core_microbiome, height = 12, width = 15)


#### Determine identity of ASVs and generate data frame ####
# Female
d_f_ASVs_v <- as.vector(d_f_ASVs) %>%
  prune_taxa(d_f)
d_f_species <- as.data.frame(tax_table(d_f_ASVs_v)) %>%
  select("Genus")
d_f_species$ASV <- rownames(sample_data(d_f_species))
rownames(d_f_species) <- NULL
colnames(d_f_species)[colnames(d_f_species) == "Genus"] <- "Diabetic_Female"

pd_f_ASVs_v <- as.vector(pd_f_ASVs) %>%
  prune_taxa(pd_f)
pd_f_species <- as.data.frame(tax_table(pd_f_ASVs_v)) %>%
  select("Genus")
pd_f_species$ASV <- rownames(sample_data(pd_f_species))
rownames(pd_f_species) <- NULL
colnames(pd_f_species)[colnames(pd_f_species) == "Genus"] <- "Pre-diabetic_Female"

nd_f_ASVs_v <- as.vector(nd_f_ASVs) %>%
  prune_taxa(nd_f)
nd_f_species <- as.data.frame(tax_table(nd_f_ASVs_v)) %>%
  select("Genus")
nd_f_species$ASV <- rownames(sample_data(nd_f_species))
rownames(nd_f_species) <- NULL
colnames(nd_f_species)[colnames(nd_f_species) == "Genus"] <- "Non-diabetic_Female"

# Male
d_m_ASVs_v <- as.vector(d_m_ASVs) %>%
  prune_taxa(d_m)
d_m_species <- as.data.frame(tax_table(d_m_ASVs_v)) %>%
  select("Genus")
d_m_species$ASV <- rownames(sample_data(d_m_species))
rownames(d_m_species) <- NULL
colnames(d_m_species)[colnames(d_m_species) == "Genus"] <- "Diabetic_Male"

pd_m_ASVs_v <- as.vector(pd_m_ASVs) %>%
  prune_taxa(pd_m)
pd_m_species <- as.data.frame(tax_table(pd_m_ASVs_v)) %>%
  select("Genus")
pd_m_species$ASV <- rownames(sample_data(pd_m_species))
rownames(pd_m_species) <- NULL
colnames(pd_m_species)[colnames(pd_m_species) == "Genus"] <- "Pre-diabetic_Male"

nd_m_ASVs_v <- as.vector(nd_m_ASVs) %>%
  prune_taxa(nd_m)
nd_m_species <- as.data.frame(tax_table(nd_m_ASVs_v)) %>%
  select("Genus")
nd_m_species$ASV <- rownames(sample_data(nd_m_species))
rownames(nd_m_species) <- NULL
colnames(nd_m_species)[colnames(nd_m_species) == "Genus"] <- "Non-diabetic_Male"


# Generating a combined data frame
combined_species <- merge(d_m_species, d_f_species, by = "ASV", all = TRUE) %>%
  merge(pd_m_species, by = "ASV", all = TRUE) %>%
  merge(pd_f_species, by = "ASV", all = TRUE) %>%
  merge(nd_m_species, by = "ASV", all = TRUE) %>%
  merge(nd_f_species, by = "ASV", all = TRUE)

combined_species_sorted <- combined_species[order(combined_species$Diabetic_Male), ]

write.csv(combined_species_sorted, "core_species_sorted.csv", row.names = FALSE)
