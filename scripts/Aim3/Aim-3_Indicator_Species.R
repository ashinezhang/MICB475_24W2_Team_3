# Load packages
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(vegan)

#### Load in RData ####
load("data/colombia_final.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
colombia_genus <- tax_glom(colombia_final, "Genus", NArm = FALSE)
colombia_genus_RA <- transform_sample_counts(colombia_genus, fun=function(x) x/sum(x))

#ISA
isa_colombia <- multipatt(t(otu_table(colombia_genus_RA)), cluster = sample_data(colombia_genus_RA)$`diabetic_status_and_sex`)
summary(isa_colombia)
taxtable <- tax_table(colombia_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Make table
isa_table <- isa_colombia$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05, stat>0.7)

write.csv(isa_table, "sex_isa_table_sig.csv", row.names = FALSE)
