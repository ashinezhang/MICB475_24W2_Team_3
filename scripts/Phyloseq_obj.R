library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)

#### Load data ####
metafp <- "data/colombia_export/metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "data/colombia_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "data/colombia_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "data/colombia_export/tree.nwk"
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

###Filtering Steps###
# Remove non-bacterial sequences, if any
colombia_filt <- subset_taxa(colombia,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
colombia_filt_nolow <- filter_taxa(colombia_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
colombia_final<- prune_samples(sample_sums(colombia_filt_nolow)>100, colombia_filt_nolow)

###Rarefaction###
rarecurve(t(as.data.frame(otu_table(colombia_final))), cex=0.1)
colombia_rare <- rarefy_even_depth(colombia_final, rngseed = 1, sample.size = 5900)

##### Saving #####
save(colombia_final, file="colombia_final.RData")
save(colombia_rare, file="colombia_rare.RData")


