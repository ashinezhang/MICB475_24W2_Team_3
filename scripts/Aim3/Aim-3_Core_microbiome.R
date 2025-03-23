# Load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

colombia_meta <- read_tsv("data/fixed_colombia_diabetes_metadata.tsv")

load("./data/colombia_final.RData")

# Convert to relative abundance
colombia_RA <- transform_sample_counts(colombia_final, fun=function(x) x/sum(x))

# Filter dataset by sex and diabetic status
d_f <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Diabetic_female")
pd_f <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Pre-diabetic_female")
nd_f <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Non-diabetic_female")

d_m <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Diabetic_male")
pd_m <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Pre-diabetic_male")
nd_m <- subset_samples(colombia_RA, `diabetic_status_and_sex`=="Non-diabetic_male")

# Set prevalence and abudance threshold
d_f_ASVs <- core_members(d_f, detection=0.001, prevalence = 0.8)
pd_f_ASVs <- core_members(pd_f, detection=0.001, prevalence = 0.8)
nd_f_ASVs <- core_members(nd_f, detection=0.001, prevalence = 0.8)

d_m_ASVs <- core_members(d_m, detection=0.001, prevalence = 0.8)
pd_m_ASVs <- core_members(pd_m, detection=0.001, prevalence = 0.8)
nd_m_ASVs <- core_members(nd_m, detection=0.001, prevalence = 0.8)

# Venn-diagram pre-diabetic female vs. diabetic female
d_pd_f_core_microbiome <- ggVennDiagram(x=list(
  "Pre-diabetic Female" = pd_f_ASVs, 
  "Diabetic Female" = d_f_ASVs)) + coord_flip()
d_pd_f_core_microbiome

# Venn-diagram pre-diabetic male vs. diabetic male
d_pd_m_core_microbiome <- ggVennDiagram(x=list(
  "Pre-diabetic Male" = pd_m_ASVs, 
  "Diabetic Male" = d_m_ASVs)) + coord_flip()
d_pd_m_core_microbiome

# Venn-diagram diabetic female vs. diabetic male
d_m_f_core_microbiome <- ggVennDiagram(x=list(
  "Diabetic Male" = d_m_ASVs, 
  "Diabetic Female" = d_f_ASVs)) + coord_flip()
d_m_f_core_microbiome

# Venn-diagram pre-diabetic female vs. pre-diabetic male
pd_m_f_core_microbiome <- ggVennDiagram(x=list(
  "Pre-diabetic Male" = pd_m_ASVs, 
  "Pre-diabetic Female" = pd_f_ASVs)) + coord_flip()
pd_m_f_core_microbiome

# Venn-diagram non-diabetic female vs. non-diabetic male
nd_m_f_core_microbiome <- ggVennDiagram(x=list(
  "Non-diabetic Male" = nd_m_ASVs, 
  "Non-diabetic Female" = nd_f_ASVs)) + coord_flip()
nd_m_f_core_microbiome

ggsave("./data/Aim3/d_pd_f_core_microbiome.png", d_pd_f_core_microbiome, height = 12, width = 15)
ggsave("./data/Aim3/d_pd_m_core_microbiome.png", d_pd_m_core_microbiome, height = 12, width = 15)
ggsave("./data/Aim3/d_m_f_core_microbiome.png", d_m_f_core_microbiome, height = 12, width = 15)
ggsave("./data/Aim3/pd_m_f_core_microbiome.png", pd_m_f_core_microbiome, height = 12, width = 15)
ggsave("./data/Aim3/nd_m_f_core_microbiome.png", nd_m_f_core_microbiome, height = 12, width = 15)

# Determine identity of ASVs and plot - diabetic female vs. diabetic male
prune_taxa(d_f_ASVs,d_f) %>%
  plot_bar(fill="Genus") + 
  facet_wrap(.~`diabetic_status_and_sex`, scales = "free") %>%
prune_taxa(d_m_ASVs,d_m) %>%
  plot_bar(fill="Genus") + 
  facet_wrap(.~`diabetic_status_and_sex`, scales = "free")
  
# Determine identity of ASVs and plot - pre-diabetic female vs. pre-diabetic male
prune_taxa(pd_f_ASVs,pd_f) %>%
  plot_bar(fill="Genus") + 
  facet_wrap(.~`diabetic_status_and_sex`, scales = "free")
prune_taxa(pd_m_ASVs,pd_m) %>%
  plot_bar(fill="Genus") + 
  facet_wrap(.~`diabetic_status_and_sex`, scales = "free")

# Determine identity of ASVs and plot - non-diabetic female vs. non-diabetic male
prune_taxa(nd_f_ASVs,nd_f) %>%
  plot_bar(fill="Genus") + 
  facet_wrap(.~`diabetic_status_and_sex`, scales = "free")
prune_taxa(nd_m_ASVs,nd_m) %>%
  plot_bar(fill="Genus") + 
  facet_wrap(.~`diabetic_status_and_sex`, scales = "free")

