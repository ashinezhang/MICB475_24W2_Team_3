# Load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)

###### DESeq ######
colombia_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
colombia_deseq <- phyloseq_to_deseq2(colombia_plus1, ~`diabetic_status_and_sex`)
DESEQ_colombia <- DESeq(colombia_deseq)



#### Diabetic male vs. Diabetic female ####
# Volcano plot
res_d_m_f <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Diabetic_male","Diabetic_female"))
vol_plot_d_m_f <- res_d_m_f %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_d_m_f.png",vol_plot_d_m_f)

# Table and bar plot of results
sigASVs_d_m_f <- res_d_m_f %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_d_m_f <- sigASVs_d_m_f %>%
  pull(ASV)

d_m_f_DESeq <- prune_taxa(sigASVs_vec_d_m_f,colombia_final)
sigASVs_d_m_f <- tax_table(d_m_f_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_d_m_f) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA" & !str_detect(Genus, "^g__uncultured"))

bar_plot_d_m_f <- ggplot(sigASVs_d_m_f) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_d_m_f.png", bar_plot_d_m_f, width = 15, height = 10)



#### Pre-diabetic male vs. Pre-diabetic female ####
# Volcano plot
res_pd_m_f <- results(DESEQ_colombia, tidy=TRUE, 
                     contrast = c("diabetic_status_and_sex","Pre-diabetic_male","Pre-diabetic_female"))
vol_plot_pd_m_f <- res_pd_m_f %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_pd_m_f.png",vol_plot_pd_m_f)

# Table and bar plot of results
sigASVs_pd_m_f <- res_pd_m_f %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_pd_m_f <- sigASVs_pd_m_f %>%
  pull(ASV)

pd_m_f_DESeq <- prune_taxa(sigASVs_vec_pd_m_f,colombia_final)
sigASVs_pd_m_f <- tax_table(pd_m_f_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_pd_m_f) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA.1" & !str_detect(Genus, "^g__uncultured"))

bar_plot_pd_m_f <- ggplot(sigASVs_pd_m_f) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_pd_m_f.png", bar_plot_pd_m_f, width = 15, height = 10)



#### Non-diabetic male vs. Non-diabetic female ####
# Volcano plot
res_nd_m_f <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Non-diabetic_male","Non-diabetic_female"))
vol_plot_nd_m_f <- res_nd_m_f %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_nd_m_f.png",vol_plot_nd_m_f)

# Table and bar plot of results
sigASVs_nd_m_f <- res_nd_m_f %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_nd_m_f <- sigASVs_nd_m_f %>%
  pull(ASV)

nd_m_f_DESeq <- prune_taxa(sigASVs_vec_nd_m_f,colombia_final)
sigASVs_nd_m_f <- tax_table(nd_m_f_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_nd_m_f) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  filter(Genus != "NA" & !str_detect(Genus, "^g__uncultured"))

bar_plot_nd_m_f <- ggplot(sigASVs_nd_m_f) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_nd_m_f.png", bar_plot_nd_m_f, width = 15, height = 10)
