library(ggplot2)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(DESeq2)

#Alpha and Beta Diversity
#### Load in RData ####
load("data/colombia_rare.RData")
load("data/colombia_final.RData")

#### Alpha diversity ######
plot_richness(colombia_rare) 

plot_richness(colombia_rare, measures = c("Shannon","Chao1")) 

gg_richness <- plot_richness(colombia_rare, x = "diabetic_status", measures = c("Shannon","Chao1")) +
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

# add PD to metadata table
sample_data(colombia_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(colombia_rare), aes(diabetic_status, PD)) + 
  geom_boxplot() +
  xlab("Diabetic Status") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd


#### Beta diversity #####
bc_dm <- distance(colombia_rare, method="bray")

pcoa_bc <- ordinate(colombia_rare, method="PCoA", distance=bc_dm)

plot_ordination(colombia_rare, pcoa_bc, color = "sex", shape="diabetic_Status")

gg_pcoa <- plot_ordination(colombia_rare, pcoa_bc, color = "sex", shape="diabetic_status") +
  labs(pch="Diabetic Status #", col = "sex")
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

ggsave(filename="data/aim2/log_bar_plot.png",log_bar_plot)

