library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)

#### Load data ####
metafp <- "colombia_export/metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "colombia_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "colombia_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "colombia_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata ####
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all into a phyloseq object
colombia<- phyloseq(OTU, SAMP, TAX, phylotree)

