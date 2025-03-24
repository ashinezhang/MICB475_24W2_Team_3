library(ggpubr)
library(ggplot2)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(DESeq2)
library(indicspecies)
library(vegan)
library(reshape2)
library(pheatmap)
#Alpha and Beta Diversity
#### Load in RData ####
load("data/colombia_rare.RData")
load("data/colombia_final.RData")

#### Alpha diversity ######
plot_richness(colombia_rare) 

plot_richness(colombia_rare, measures = c("Shannon")) 

gg_richness <- plot_richness(colombia_rare, x = "diabetic_status", measures = c("Shannon")) +
  xlab("Diabetic Status") +
  geom_boxplot()
gg_richness

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
kruskal.test(phylo_dist$PD ~ diabetic_status, data = faith_pd_data)

pairwise.wilcox.test(faith_pd_data$`phylo_dist$PD`, faith_pd_data$diabetic_status,
                     p.adjust.method = "BH")

# plot any metadata category against the PD
plot_pd <- ggplot(sample_data(colombia_rare), aes(diabetic_status, phylo_dist$PD)) + xlab("Diabetic Status") + ylab("Phylogenetic Diversity") + geom_boxplot() + stat_compare_means(comparisons = list(c("Diabetic", "Pre-diabetic"), c("Pre-diabetic", "Non-diabetic"), c("Diabetic", "Non-diabetic")), method = "wilcox.test", label = "p.signif")

# view plot
plot_pd

kruskal.test( PD ~ diabetic_status, data=PD_meta)

ggsave(filename = "data/aim2/plot_PD_stats.png", plot_pd)

#### Beta diversity #####
bc_dm <- distance(colombia_rare, method="weighted_unifrac")

pcoa_bc <- ordinate(colombia_rare, method="PCoA", distance=bc_dm)
ordination_df <- as.data.frame(pcoa_bc$vectors)


plot_ordination(colombia_rare, pcoa_bc, color = "sex", shape="diabetic_Status")

#Permanova
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
gg_pcoa

ggsave("data/aim2/plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)
###Deseq###
colombia_deseq <- phyloseq_to_deseq2(colombia_final, ~`diabetic_status`)
DESEQ_mpt <- DESeq(colombia_deseq)

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
colombia_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
colombia_deseq <- phyloseq_to_deseq2(colombia_plus1, ~`diabetic_status`)
DESEQ_colombia <- DESeq(colombia_deseq)
res <- results(DESEQ_colombia, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("diabetic_status","Pre-diabetic","Diabetic"))
View(res)

##Volcano plot##
## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="data/aim2/vol_plot.png",vol_plot)

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
colombia_DESeq <- prune_taxa(sigASVs_vec,colombia_final)
sigASVs <- tax_table(colombia_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log_bar_plot<-ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

log_bar_plot

ggsave(filename="data/aim2/log_bar_plot.png",log_bar_plot)

#### Indicator Species/Taxa Analysis ####
# glom to Genus
colombia_genus <- tax_glom(colombia_final, "Genus", NArm = FALSE)
colombia_genus_RA <- transform_sample_counts(colombia_genus, fun=function(x) x/sum(x))

#ISA
isa_colombia <- multipatt(t(otu_table(colombia_genus_RA)), cluster = sample_data(colombia_genus_RA)$`diabetic_status`)
summary(isa_colombia)
taxtable <- tax_table(colombia_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_table<-isa_colombia$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 
isa_table
isa_table<-as.data.frame(isa_table)

write.csv(isa_table, file = "data/aim2/isa_table.csv", row.names = FALSE)


