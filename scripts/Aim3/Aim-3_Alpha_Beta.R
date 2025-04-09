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
  stat_compare_means(comparisons = list(c("Diabetic_female","Diabetic_male"), c("Non-diabetic_female","Non-diabetic_male"), c("Pre-diabetic_female","Pre-diabetic_male")), method = "wilcox.test",label = "p.signif", p.adjust.method = "BH")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_discrete(labels = function(x) gsub("_", " ", x))
plot_shannon

ggsave(filename = "data/Aim3/shannon_stats.png"
       , plot_shannon
       , height=5, width=6)

#Chao1 with stats
#stats
chao1_values <- estimate_richness(colombia_rare, measures = "Chao1")
chao1_data <- cbind(metadata, chao1_values$Chao1)
wilcox_chao1<-pairwise.wilcox.test(chao1_data$`chao1_values$Chao1`, chao1_data$diabetic_status_and_sex,
                                     p.adjust.method = "BH")
print(wilcox_chao1)
#Shannon plot with stats
plot_chao1 <- ggplot(sample_data(colombia_rare), aes(diabetic_status_and_sex,chao1_values$Chao1)) + xlab("") + ylab("Chao1 Index") + geom_boxplot() + 
  stat_compare_means(comparisons = list(c("Diabetic_female","Diabetic_male"), c("Non-diabetic_female","Non-diabetic_male"), c("Pre-diabetic_female","Pre-diabetic_male")), method = "wilcox.test",label = "p.signif", p.adjust.method = "BH")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(NA, 475)
plot_chao1

ggsave(filename = "data/Aim3/chao1_stats.png"
       , plot_chao1
       , height=6, width=6)

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
  stat_compare_means(comparisons = list(c("Diabetic_female","Diabetic_male"), c("Non-diabetic_female","Non-diabetic_male"), c("Pre-diabetic_female","Pre-diabetic_male")), method = "wilcox.test",label = "p.signif", p.adjust.method = "BH")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(NA, 25)


plot_pd

ggsave(filename = "data/Aim3/plot_PD_stats.png", plot_pd, height=4, width=6)

##Note: Diabetic female and male are actually non-significant when BH is used to account for FDR, but stat compare means won't show it on the plot##

### Beta diversity #####
bc_dm <- distance(colombia_rare, method="weighted_unifrac")

pcoa_bc <- ordinate(colombia_rare, method="PCoA", distance=bc_dm)
ordination_df <- as.data.frame(pcoa_bc$vectors)


plot_ordination(colombia_rare, pcoa_bc, color = "sex", shape="diabetic_Status")

#Permanova diabetetic status and sex
dm_unifrac <- UniFrac(colombia_rare, weighted=TRUE) # Weighted UniFrac
data <- as(sample_data(colombia_rare), "data.frame")
adonis_results<-adonis2(dm_unifrac ~ diabetic_status*sex, data=data)
adonis2(dm_unifrac ~ diabetic_status, data=data)
R2_value <- round(adonis_results$R2[1], 3)
p_value <- adonis_results$`Pr(>F)`[1]


gg_pcoa <- plot_ordination(colombia_rare, pcoa_bc, color = "sex", shape="diabetic_status") +
  labs(pch="Diabetic Status #", col = "sex")+
  annotate(geom ="text", 
           x = min(ordination_df$Axis.1), 
           y = max(ordination_df$Axis.2), 
           label = paste("PERMANOVA: R² =", R2_value, "\n p =", p_value), 
           hjust = 0, size = 3)
