library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Alpha diversity ######

shannon_chao1 <- plot_richness(colombia_rare, x = "diabetic_status_and_sex", measures = c("Shannon","Chao1"))

gg_richness <- shannon_chao1 +
  xlab("diabetic_status_and_sex") +
  geom_boxplot()
gg_richness

ggsave(filename = "Aim3_Shannon_plot_richness.png"
       , gg_richness
       , height=4, width=6)

richness_table <- as.data.frame(shannon_chao1)

mdl_richness <- lm("Shannon" ~ "diabetic_status_and_sex", data = shannon_chao1)
aov_mdl_richness <- aov(mdl_richness)

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(colombia_rare)), phy_tree(colombia_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(colombia_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(colombia_rare), aes(diabetic_status_and_sex, PD)) + 
  geom_boxplot() +
  xlab("Diabetic Status and Sex") +
  ylab("Phylogenetic Diversity")



ggsave(filename = "Aim3_Faiths_phylogenetic_diversity.png"
       , plot.pd
       , height=4, width=6)

# view plot
plot.pd


#### Beta diversity #####
wu_dm <- distance(colombia_rare, method="wunifrac")

pcoa_wu <- ordinate(colombia_rare, method="PCoA", distance=wu_dm)


sex_pcoa <- plot_ordination(colombia_rare, pcoa_wu, color = "sex", shape="diabetic_status") +
  labs(pch="Diabetic Status", col = "Sex")
sex_pcoa

status_pcoa <- plot_ordination(colombia_rare, pcoa_wu, color = "diabetic_status", shape="sex") +
  labs(pch="Sex", col = "Diabetic Status")
status_pcoa

ggsave("Aim3_sex_pcoa.png"
       , sex_pcoa
       , height=4, width=5)

ggsave("Aim3_status_pcoa.png"
       , status_pcoa
       , height=4, width=5)

#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(colombia_rare, fill="Phylum") 

# Convert to relative abundance
colombia_RA <- transform_sample_counts(colombia_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
colombia_phylum <- tax_glom(colombia_RA, taxrank = "Phylum", NArm=FALSE)

plot_bar(colombia_phylum, fill="Phylum") + 
  facet_wrap(.~subject, scales = "free_x")

gg_taxa <- plot_bar(colombia_phylum, fill="Phylum") + 
  facet_wrap(.~subject, scales = "free_x")
gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)
