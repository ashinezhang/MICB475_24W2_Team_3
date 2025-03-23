# Load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)

###### DESeq ######
colombia_plus1 <- transform_sample_counts(colombia_final, function(x) x+1)
colombia_deseq <- phyloseq_to_deseq2(colombia_plus1, ~`diabetic_status_and_sex`)
DESEQ_colombia <- DESeq(colombia_deseq)

#### Diabetic female vs. Pre-diabetic female ####
# Volcano plot
res_d_pd_f <- results(DESEQ_colombia, tidy=TRUE, 
               contrast = c("diabetic_status_and_sex","Diabetic_female","Pre-diabetic_female"))
vol_plot_d_pd_f <- res_d_pd_f %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_d_pd_f.png",vol_plot_d_pd_f)

# Table and bar plot of results
sigASVs_d_pd_f <- res_d_pd_f %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_d_pd_f <- sigASVs_d_pd_f %>%
  pull(ASV)

d_pd_f_DESeq <- prune_taxa(sigASVs_vec_d_pd_f,colombia_final)
sigASVs_d_pd_f <- tax_table(d_pd_f_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_d_pd_f) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_d_pd_f <- ggplot(sigASVs_d_pd_f) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_d_pd_f.png", bar_plot_d_pd_f, width = 15, height = 10)


#### Diabetic female vs. Non-diabetic female ####
# Volcano plot
res_d_nd_f <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Diabetic_female","Non-diabetic_female"))
vol_plot_d_nd_f <- res_d_nd_f %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_d_nd_f.png",vol_plot_d_nd_f)

# Table and bar plot of results
sigASVs_d_nd_f <- res_d_nd_f %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_d_nd_f <- sigASVs_d_nd_f %>%
  pull(ASV)

d_nd_f_DESeq <- prune_taxa(sigASVs_vec_d_nd_f,colombia_final)
sigASVs_d_nd_f <- tax_table(d_nd_f_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_d_nd_f) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_d_nd_f <- ggplot(sigASVs_d_nd_f) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_d_nd_f.png", bar_plot_d_nd_f, width = 15, height = 10)


#### Pre-diabetic female vs. Non-diabetic female ####
# Volcano plot
res_pd_nd_f <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Pre-diabetic_female","Non-diabetic_female"))
vol_plot_pd_nd_f <- res_pd_nd_f %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_pd_nd_f.png",vol_plot_pd_nd_f)

# Table and bar plot of results
sigASVs_pd_nd_f <- res_pd_nd_f %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_pd_nd_f <- sigASVs_pd_nd_f %>%
  pull(ASV)

pd_nd_f_DESeq <- prune_taxa(sigASVs_vec_pd_nd_f,colombia_final)
sigASVs_pd_nd_f <- tax_table(pd_nd_f_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_pd_nd_f) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_pd_nd_f <- ggplot(sigASVs_pd_nd_f) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_pd_nd_f.png", bar_plot_pd_nd_f, width = 15, height = 10)


#### Diabetic male vs. Pre-diabetic male ####
# Volcano plot
res_d_pd_m <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Diabetic_male","Pre-diabetic_male"))
vol_plot_d_pd_m <- res_d_pd_m %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_d_pd_m.png",vol_plot_d_pd_m)

# Table and bar plot of results
sigASVs_d_pd_m <- res_d_pd_m %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_d_pd_m <- sigASVs_d_pd_m %>%
  pull(ASV)

d_pd_m_DESeq <- prune_taxa(sigASVs_vec_d_pd_m,colombia_final)
sigASVs_d_pd_m <- tax_table(d_pd_m_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_d_pd_m) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_d_pd_m <- ggplot(sigASVs_d_pd_m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_d_pd_m.png", bar_plot_d_pd_m, width = 15, height = 10)


#### Diabetic male vs. Non-diabetic male ####
# Volcano plot
res_d_nd_m <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Diabetic_male","Non-diabetic_male"))
vol_plot_d_nd_m <- res_d_nd_m %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_d_nd_m.png",vol_plot_d_nd_m)

# Table and bar plot of results
sigASVs_d_nd_m <- res_d_nd_m %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_d_nd_m <- sigASVs_d_nd_m %>%
  pull(ASV)

d_nd_m_DESeq <- prune_taxa(sigASVs_vec_d_nd_m,colombia_final)
sigASVs_d_nd_m <- tax_table(d_nd_m_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_d_nd_m) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_d_nd_m <- ggplot(sigASVs_d_nd_m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_d_nd_m.png", bar_plot_d_nd_m, width = 15, height = 10)


#### Pre-diabetic male vs. Non-diabetic male ####
# Volcano plot
res_pd_nd_m <- results(DESEQ_colombia, tidy=TRUE, 
                       contrast = c("diabetic_status_and_sex","Pre-diabetic_male","Non-diabetic_male"))
vol_plot_pd_nd_m <- res_pd_nd_m %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_pd_nd_m.png",vol_plot_pd_nd_m)

# Table and bar plot of results
sigASVs_pd_nd_m <- res_pd_nd_m %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_pd_nd_m <- sigASVs_pd_nd_m %>%
  pull(ASV)

pd_nd_m_DESeq <- prune_taxa(sigASVs_vec_pd_nd_m,colombia_final)
sigASVs_pd_nd_m <- tax_table(pd_nd_m_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_pd_nd_m) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_pd_nd_m <- ggplot(sigASVs_pd_nd_m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_pd_nd_m.png", bar_plot_pd_nd_m, width = 15, height = 10)


#### Diabetic female vs. Diabetic male ####
# Volcano plot
res_d_f_m <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Diabetic_female","Diabetic_male"))
vol_plot_d_f_m <- res_d_f_m %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_d_f_m.png",vol_plot_d_f_m)

# Table and bar plot of results
sigASVs_d_f_m <- res_d_f_m %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_d_f_m <- sigASVs_d_f_m %>%
  pull(ASV)

d_f_m_DESeq <- prune_taxa(sigASVs_vec_d_f_m,colombia_final)
sigASVs_d_f_m <- tax_table(d_f_m_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_d_f_m) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_d_f_m <- ggplot(sigASVs_d_f_m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_d_f_m.png", bar_plot_d_f_m, width = 15, height = 10)


#### Pre-diabetic female vs. Pre-diabetic male ####
# Volcano plot
res_pd_f_m <- results(DESEQ_colombia, tidy=TRUE, 
                     contrast = c("diabetic_status_and_sex","Pre-diabetic_female","Pre-diabetic_male"))
vol_plot_pd_f_m <- res_pd_f_m %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_pd_f_m.png",vol_plot_pd_f_m)

# Table and bar plot of results
sigASVs_pd_f_m <- res_pd_f_m %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_pd_f_m <- sigASVs_pd_f_m %>%
  pull(ASV)

pd_f_m_DESeq <- prune_taxa(sigASVs_vec_pd_f_m,colombia_final)
sigASVs_pd_f_m <- tax_table(pd_f_m_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_pd_f_m) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_pd_f_m <- ggplot(sigASVs_pd_f_m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_pd_f_m.png", bar_plot_pd_f_m, width = 15, height = 10)


#### Non-diabetic female vs. Non-diabetic male ####
# Volcano plot
res_nd_f_m <- results(DESEQ_colombia, tidy=TRUE, 
                      contrast = c("diabetic_status_and_sex","Non-diabetic_female","Non-diabetic_male"))
vol_plot_nd_f_m <- res_nd_f_m %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="vol_plot_nd_f_m.png",vol_plot_nd_f_m)

# Table and bar plot of results
sigASVs_nd_f_m <- res_nd_f_m %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_vec_nd_f_m <- sigASVs_nd_f_m %>%
  pull(ASV)

nd_f_m_DESeq <- prune_taxa(sigASVs_vec_nd_f_m,colombia_final)
sigASVs_nd_f_m <- tax_table(nd_f_m_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_nd_f_m) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

bar_plot_nd_f_m <- ggplot(sigASVs_nd_f_m) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename="DESeq_nd_f_m.png", bar_plot_nd_f_m, width = 15, height = 10)
