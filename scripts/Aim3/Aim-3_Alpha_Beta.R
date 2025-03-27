getwd()
library(ggpubr)
library(ggplot2)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#Alpha and Beta Diversity
#### Load in RData ####
load("data/colombia_rare.RData")
load("data/colombia_final.RData")

#### Alpha diversity ######
=

gg_richness <- plot_richness(colombia_rare, x = "diabetic_status_and_sex", measures = c("Shannon", "Chao1")) +
  xlab("Diabetic Status and Sex") +
  geom_boxplot()
gg_richness

#stats
shannon_values <- estimate_richness(colombia_rare, measures = "Shannon")
metadata <- data.frame(sample_data(colombia_rare))
shannon_data <- cbind(metadata, shannon_values$Shannon)
wilcox_shannon<-pairwise.wilcox.test(shannon_data$`shannon_values$Shannon`, shannon_data$diabetic_status_and_sex,
                                     p.adjust.method = "BH")
print(wilcox_shannon)
#Shannon plot with stats
plot_shannon <- ggplot(sample_data(colombia_rare), aes(diabetic_status_and_sex,shannon_values$Shannon)) + xlab("") + ylab("Shannon Index") + geom_boxplot() + 
  stat_compare_means(comparisons = list(c("Diabetic_female","Diabetic_male"), c("Non-diabetic_female","Non-diabetic_male"), c("Pre-diabetic_female","Pre-diabetic_male")), method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_discrete(labels = function(x) gsub("_", " ", x))
plot_shannon

ggsave(filename = "data/Aim3/shannon_stats.png"
       , gg_richness
       , height=4, width=6)

#IDK how

#visualize
ggsave(filename = "data/aim2/plot_richness.png"
       , gg_richness
       , height=4, width=6)

estimate_richness(colombia_rare)

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(colombia_rare)), phy_tree(colombia_rare),
                 include.root=F) 

# add PD to metadata table and add stats
sample_data(colombia_rare)$PD <- phylo_dist$PD

metadata_PD <- data.frame(sample_data(colombia_rare))
faith_pd_data <- cbind(metadata_PD, phylo_dist$PD)
faith_pd_data
kruskal.test(phylo_dist$PD ~ diabetic_status_and_sex, data = faith_pd_data)

wilcox_results<-pairwise.wilcox.test(faith_pd_data$`phylo_dist$PD`, faith_pd_data$diabetic_status_and_sex,
                     p.adjust.method = "BH")
print(wilcox_results)
# plot any metadata category against the PD
plot_pd <- ggplot(sample_data(colombia_rare), aes(diabetic_status_and_sex, phylo_dist$PD)) + xlab("") + ylab("Phylogenetic Diversity") + geom_boxplot() + 
  stat_compare_means(comparisons = list(c("Diabetic_female","Diabetic_male"), c("Non-diabetic_female","Non-diabetic_male"), c("Pre-diabetic_female","Pre-diabetic_male")), method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_discrete(labels = function(x) gsub("_", " ", x))

plot_pd

ggsave(filename = "data/Aim3/plot_PD_stats.png", plot_pd)

### Beta diversity #####
bc_dm <- distance(colombia_rare, method="weighted_unifrac")

pcoa_bc <- ordinate(colombia_rare, method="PCoA", distance=bc_dm)
ordination_df <- as.data.frame(pcoa_bc$vectors)


plot_ordination(colombia_rare, pcoa_bc, color = "sex", shape="diabetic_Status")
