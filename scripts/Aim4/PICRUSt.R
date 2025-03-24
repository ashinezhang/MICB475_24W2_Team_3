##### Install packages #####
# Start by installing all necessary packages 
# When asked if you want to install from source, type Yes in the terminal below

# If you don't have BiocManager, here is the code to install it
# A lot of you probably already have this so you can skip
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# this can take very long to load, 
# when asked if you want to update, type "n" for none to make it faster

# After installing all of its above dependencies, install ggpicrust2
# install.packages("ggpicrust2")
### On Nov 25, 2024, ggpicrust2 was archived from cran, so it can't install
### https://github.com/cafferychen777/ggpicrust2?tab=readme-ov-file#installation
install.packages("devtools")
devtools::install_github("cafferychen777/ggpicrust2")

#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

#### Import files and preparing tables ####
abundance_file <- "data/Aim4/pathabun_exported/pathway_abundance.tsv"

### TROUBLESHOOTING - read_delim importing as a single column, not recognizing tabs, trying read.table instead
abundance_data <- read.table("data/Aim4/pathabun_exported/pathway_abundance.tsv", 
                             sep = "\t", header = FALSE, check.names = FALSE, quote = "", 
                             comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
# Remove the first row which just says # Constructed from biom file
abundance_data <- abundance_data[-1, ]
# Assign the second row as column names
colnames(abundance_data) <- as.character(abundance_data[1, ])
abundance_data <- abundance_data[-1, ]
# Change column 1 name to pathway
colnames(abundance_data)[colnames(abundance_data) == "#OTU ID"] <- "pathway"

# Convert to a data frame
abundance_data = as.data.frame(abundance_data, check.names = FALSE)

#Import your metadata file, no need to filter yet
metadata <- read_delim("data/colombia_export/metadata.tsv")

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$diabetic_status_and_sex),]

#can only compare 2 at a time, so total 9 analyses 
# 1. d_nd_f, 2. d_pd_f, 3. pd_nd_f, 4. d_nd_m, 5. d_pd_m, 6. pd_nd_m, 7. d_f_m, 8. pd_f_m, 9. nd_f_m

#### 1. diabetic female vs non-diabetic female
metadata_d_nd_f <- filter(metadata, diabetic_status_and_sex == "Diabetic_female" | diabetic_status_and_sex == "Non-diabetic_female")

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names_d_nd_f = metadata_d_nd_f$'sample-id'
sample_names_d_nd_f = append(sample_names_d_nd_f, "pathway")
abundance_data_filtered_d_nd_f = abundance_data[, colnames(abundance_data) %in% sample_names_d_nd_f] 

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered_d_nd_f = abundance_data_filtered_d_nd_f[, colSums(abundance_data_filtered_d_nd_f != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered_d_nd_f) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples_d_nd_f = rownames(t(abundance_data_filtered_d_nd_f[,-1])) #Getting a list of the sample names in the newly filtered abundance data
# metadata_d_nd_f = metadata[metadata$`sample-id` %in% abun_samples_d_nd_f,] #making sure the filtered metadata only includes these samples

### TROUBLESHOOTING
# For some reason there's one more in abundance than in metadata
# Get the sample IDs from both data sources
abundance_sample_ids <- colnames(abundance_data_filtered_d_nd_f)
metadata_sample_ids <- metadata_d_nd_f$`sample-id`
# Find the extra sample(s) in abundance data
extra_samples <- setdiff(abundance_sample_ids, metadata_sample_ids)
print(extra_samples) # should only have 'pathway' as an extra sample
### commented out a line of code because it was deleting sample ...383 from metadata_d_nd_f even though it was still in the abun_samples_d_nd_f

#### DESEq ####

### TROUBLESHOOTING - pathway DAA using DESEQ2 says that there are non-numeric values
abundance_data_filtered_d_nd_f <- abundance_data_filtered_d_nd_f %>%
  mutate(across(-pathway, as.numeric))
### Converted the abundance data to numeric values, but leave pathway column as non-numeric

#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df_d_nd_f <- pathway_daa(abundance = abundance_data_filtered_d_nd_f %>% column_to_rownames("pathway"), 
                                        metadata = metadata_d_nd_f, group = "diabetic_status_and_sex", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df_d_nd_f <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df_d_nd_f, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05_d_nd_f <- abundance_daa_results_df_d_nd_f %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc_d_nd_f = inner_join(feature_with_p_0.05_d_nd_f,metacyc_daa_annotated_results_df_d_nd_f, by = "feature")
feature_desc_d_nd_f$feature = feature_desc_d_nd_f$description
feature_desc_d_nd_f = feature_desc_d_nd_f[,c(1:7)]
colnames(feature_desc_d_nd_f) = colnames(feature_with_p_0.05_d_nd_f)

#Changing the pathway column to description for the abundance table
abundance_d_nd_f = abundance_data_filtered_d_nd_f %>% filter(pathway %in% feature_with_p_0.05_d_nd_f$feature)
colnames(abundance_d_nd_f)[1] = "feature"
abundance_desc_d_nd_f = inner_join(abundance_d_nd_f,metacyc_daa_annotated_results_df_d_nd_f, by = "feature")
abundance_desc_d_nd_f$feature = abundance_desc_d_nd_f$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc_d_nd_f = abundance_desc_d_nd_f[,-c(34:ncol(abundance_desc_d_nd_f))] 

# Generate a heatmap
pathway_heatmap_d_nd_f <- pathway_heatmap(abundance = abundance_desc_d_nd_f %>% column_to_rownames("feature"), metadata = metadata_d_nd_f, group = "diabetic_status_and_sex")



### TROUBLESHOOTING - metadata must contain a 'sample_name' column for pathway PCA plot
colnames(metadata_d_nd_f)[colnames(metadata_d_nd_f) == "sample-id"] <- "sample_name"
### new problem - there are constant/zero variance columns that PCA can't work with so need to remove?
# Generate pathway PCA plot
pathway_pca_d_nd_f <-pathway_pca(abundance = abundance_data_filtered_d_nd_f %>% column_to_rownames("pathway"), metadata = metadata_d_nd_f, group = "diabetic_status_and_sex")



# Generating a bar plot representing log2FC from the custom deseq2 function
# Go to the Deseq2 function script and update the metadata category of interest
# Lead the function in
source("scripts/Aim4/DESeq2_function.R")

# Run the function on your own data
res_d_nd_f =  DEseq2_function(abundance_data_filtered_d_nd_f, metadata_d_nd_f, "diabetic_status_and_sex")
res_d_nd_f$feature =rownames(res_d_nd_f)
res_desc_d_nd_f = inner_join(res_d_nd_f,metacyc_daa_annotated_results_df_d_nd_f, by = "feature")
res_desc_d_nd_f = res_desc_d_nd_f[, -c(8:13)]
View(res_desc_d_nd_f)

# Filter to only include significant pathways
sig_res_d_nd_f = res_desc_d_nd_f %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change

sig_res_d_nd_f <- sig_res_d_nd_f[order(sig_res_d_nd_f$log2FoldChange),]
Log2FC_d_nd_f <- ggplot(data = sig_res_d_nd_f, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# Save all plots to data folder 
ggsave(filename="data/Aim4/pathway_heatmap_d_nd_f.png", pathway_heatmap_d_nd_f, width = 15, height = 10)
ggsave(filename="data/Aim4/Log2FC_d_nd_f.png", Log2FC_d_nd_f, width = 15, height = 10)


#### 2. diabetic female vs pre-diabetic female
metadata_d_pd_f <- filter(metadata, diabetic_status_and_sex == "Diabetic_female" | diabetic_status_and_sex == "Pre-diabetic_female")

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names_d_pd_f = metadata_d_pd_f$'sample-id'
sample_names_d_pd_f = append(sample_names_d_pd_f, "pathway")
abundance_data_filtered_d_pd_f = abundance_data[, colnames(abundance_data) %in% sample_names_d_pd_f] 

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered_d_pd_f = abundance_data_filtered_d_pd_f[, colSums(abundance_data_filtered_d_pd_f != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered_d_pd_f) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples_d_pd_f = rownames(t(abundance_data_filtered_d_pd_f[,-1])) #Getting a list of the sample names in the newly filtered abundance data
# metadata_d_nd_f = metadata[metadata$`sample-id` %in% abun_samples_d_nd_f,] #making sure the filtered metadata only includes these samples

### TROUBLESHOOTING
# For some reason there's one more in abundance than in metadata
# Get the sample IDs from both data sources
abundance_sample_ids <- colnames(abundance_data_filtered_d_pd_f)
metadata_sample_ids <- metadata_d_pd_f$`sample-id`
# Find the extra sample(s) in abundance data
extra_samples <- setdiff(abundance_sample_ids, metadata_sample_ids)
print(extra_samples) # should only have 'pathway' as an extra sample
### commented out a line of code because it was deleting sample ...383 from metadata_d_nd_f even though it was still in the abun_samples_d_nd_f

#### DESEq ####

### TROUBLESHOOTING - pathway DAA using DESEQ2 says that there are non-numeric values
abundance_data_filtered_d_pd_f <- abundance_data_filtered_d_pd_f %>%
  mutate(across(-pathway, as.numeric))
### Converted the abundance data to numeric values, but leave pathway column as non-numeric

#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df_d_pd_f <- pathway_daa(abundance = abundance_data_filtered_d_pd_f %>% column_to_rownames("pathway"), 
                                               metadata = metadata_d_pd_f, group = "diabetic_status_and_sex", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df_d_pd_f <- pathway_annotation(pathway = "MetaCyc", 
                                                              daa_results_df = abundance_daa_results_df_d_pd_f, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05_d_pd_f <- abundance_daa_results_df_d_pd_f %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc_d_pd_f = inner_join(feature_with_p_0.05_d_pd_f,metacyc_daa_annotated_results_df_d_pd_f, by = "feature")
feature_desc_d_pd_f$feature = feature_desc_d_pd_f$description
feature_desc_d_pd_f = feature_desc_d_pd_f[,c(1:7)]
colnames(feature_desc_d_pd_f) = colnames(feature_with_p_0.05_d_pd_f)

#Changing the pathway column to description for the abundance table
abundance_d_pd_f = abundance_data_filtered_d_pd_f %>% filter(pathway %in% feature_with_p_0.05_d_pd_f$feature)
colnames(abundance_d_pd_f)[1] = "feature"
abundance_desc_d_pd_f = inner_join(abundance_d_pd_f,metacyc_daa_annotated_results_df_d_pd_f, by = "feature")
abundance_desc_d_pd_f$feature = abundance_desc_d_pd_f$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc_d_pd_f = abundance_desc_d_pd_f[,-c(34:ncol(abundance_desc_d_pd_f))] 

# Generate a heatmap
pathway_heatmap_d_pd_f <- pathway_heatmap(abundance = abundance_desc_d_pd_f %>% column_to_rownames("feature"), metadata = metadata_d_pd_f, group = "diabetic_status_and_sex")



### TROUBLESHOOTING - metadata must contain a 'sample_name' column for pathway PCA plot
colnames(metadata_d_pd_f)[colnames(metadata_d_pd_f) == "sample-id"] <- "sample_name"
### new problem - there are constant/zero variance columns that PCA can't work with so need to remove?
# Generate pathway PCA plot
pathway_pca_d_pd_f <-pathway_pca(abundance = abundance_data_filtered_d_pd_f %>% column_to_rownames("pathway"), metadata = metadata_d_pd_f, group = "diabetic_status_and_sex")



# Generating a bar plot representing log2FC from the custom deseq2 function
# Go to the Deseq2 function script and update the metadata category of interest
# Lead the function in
source("scripts/Aim4/DESeq2_function.R")

# Run the function on your own data
res_d_pd_f =  DEseq2_function(abundance_data_filtered_d_pd_f, metadata_d_pd_f, "diabetic_status_and_sex")
res_d_pd_f$feature =rownames(res_d_pd_f)
res_desc_d_pd_f = inner_join(res_d_pd_f,metacyc_daa_annotated_results_df_d_pd_f, by = "feature")
res_desc_d_pd_f = res_desc_d_pd_f[, -c(8:13)]
View(res_desc_d_pd_f)

# Filter to only include significant pathways
sig_res_d_pd_f = res_desc_d_pd_f %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change

sig_res_d_pd_f <- sig_res_d_pd_f[order(sig_res_d_pd_f$log2FoldChange),]
Log2FC_d_pd_f <- ggplot(data = sig_res_d_pd_f, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# Save all plots to data folder 
ggsave(filename="data/Aim4/pathway_heatmap_d_pd_f.png", pathway_heatmap_d_pd_f, width = 15, height = 10)
ggsave(filename="data/Aim4/Log2FC_d_pd_f.png", Log2FC_d_pd_f, width = 15, height = 10)







#### 3. pre-diabetic female vs non-diabetic female
metadata_pd_nd_f <- filter(metadata, diabetic_status_and_sex == "Pre-diabetic_female" | diabetic_status_and_sex == "Non-diabetic_female")

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names_pd_nd_f = metadata_pd_nd_f$'sample-id'
sample_names_pd_nd_f = append(sample_names_pd_nd_f, "pathway")
abundance_data_filtered_pd_nd_f = abundance_data[, colnames(abundance_data) %in% sample_names_pd_nd_f] 

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered_pd_nd_f = abundance_data_filtered_pd_nd_f[, colSums(abundance_data_filtered_pd_nd_f != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered_pd_nd_f) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples_pd_nd_f = rownames(t(abundance_data_filtered_pd_nd_f[,-1])) #Getting a list of the sample names in the newly filtered abundance data
# metadata_d_nd_f = metadata[metadata$`sample-id` %in% abun_samples_d_nd_f,] #making sure the filtered metadata only includes these samples

### TROUBLESHOOTING
# For some reason there's one more in abundance than in metadata
# Get the sample IDs from both data sources
abundance_sample_ids <- colnames(abundance_data_filtered_pd_nd_f)
metadata_sample_ids <- metadata_pd_nd_f$`sample-id`
# Find the extra sample(s) in abundance data
extra_samples <- setdiff(abundance_sample_ids, metadata_sample_ids)
print(extra_samples) # should only have 'pathway' as an extra sample
### commented out a line of code because it was deleting sample ...383 from metadata_d_nd_f even though it was still in the abun_samples_d_nd_f

#### DESEq ####

### TROUBLESHOOTING - pathway DAA using DESEQ2 says that there are non-numeric values
abundance_data_filtered_pd_nd_f <- abundance_data_filtered_pd_nd_f %>%
  mutate(across(-pathway, as.numeric))
### Converted the abundance data to numeric values, but leave pathway column as non-numeric

#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df_pd_nd_f <- pathway_daa(abundance = abundance_data_filtered_pd_nd_f %>% column_to_rownames("pathway"), 
                                               metadata = metadata_pd_nd_f, group = "diabetic_status_and_sex", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df_pd_nd_f <- pathway_annotation(pathway = "MetaCyc", 
                                                              daa_results_df = abundance_daa_results_df_pd_nd_f, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05_pd_nd_f <- abundance_daa_results_df_pd_nd_f %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc_pd_nd_f = inner_join(feature_with_p_0.05_pd_nd_f,metacyc_daa_annotated_results_df_pd_nd_f, by = "feature")
feature_desc_pd_nd_f$feature = feature_desc_pd_nd_f$description
feature_desc_pd_nd_f = feature_desc_pd_nd_f[,c(1:7)]
colnames(feature_desc_pd_nd_f) = colnames(feature_with_p_0.05_pd_nd_f)

#Changing the pathway column to description for the abundance table
abundance_pd_nd_f = abundance_data_filtered_pd_nd_f %>% filter(pathway %in% feature_with_p_0.05_pd_nd_f$feature)
colnames(abundance_pd_nd_f)[1] = "feature"
abundance_desc_pd_nd_f = inner_join(abundance_pd_nd_f,metacyc_daa_annotated_results_df_pd_nd_f, by = "feature")
abundance_desc_pd_nd_f$feature = abundance_desc_pd_nd_f$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc_pd_nd_f = abundance_desc_pd_nd_f[,-c(34:ncol(abundance_desc_pd_nd_f))] 

# Generate a heatmap
pathway_heatmap_pd_nd_f <- pathway_heatmap(abundance = abundance_desc_pd_nd_f %>% column_to_rownames("feature"), metadata = metadata_pd_nd_f, group = "diabetic_status_and_sex")



### TROUBLESHOOTING - metadata must contain a 'sample_name' column for pathway PCA plot
colnames(metadata_d_pd_f)[colnames(metadata_d_pd_f) == "sample-id"] <- "sample_name"
### new problem - there are constant/zero variance columns that PCA can't work with so need to remove?
# Generate pathway PCA plot
pathway_pca_d_pd_f <-pathway_pca(abundance = abundance_data_filtered_d_pd_f %>% column_to_rownames("pathway"), metadata = metadata_d_pd_f, group = "diabetic_status_and_sex")



# Generating a bar plot representing log2FC from the custom deseq2 function
# Go to the Deseq2 function script and update the metadata category of interest
# Lead the function in
source("scripts/Aim4/DESeq2_function.R")

# Run the function on your own data
res_pd_nd_f =  DEseq2_function(abundance_data_filtered_pd_nd_f, metadata_pd_nd_f, "diabetic_status_and_sex")
res_pd_nd_f$feature =rownames(res_pd_nd_f)
res_desc_pd_nd_f = inner_join(res_pd_nd_f,metacyc_daa_annotated_results_df_pd_nd_f, by = "feature")
res_desc_pd_nd_f = res_desc_pd_nd_f[, -c(8:13)]
View(res_desc_pd_nd_f)

# Filter to only include significant pathways
sig_res_pd_nd_f = res_desc_pd_nd_f %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change

sig_res_pd_nd_f <- sig_res_pd_nd_f[order(sig_res_pd_nd_f$log2FoldChange),]
Log2FC_pd_nd_f <- ggplot(data = sig_res_pd_nd_f, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# Save all plots to data folder 
ggsave(filename="data/Aim4/pathway_heatmap_pd_nd_f.png", pathway_heatmap_pd_nd_f, width = 15, height = 10)
ggsave(filename="data/Aim4/Log2FC_pd_nd_f.png", Log2FC_pd_nd_f, width = 15, height = 10)




#### 4. diabetic male vs non-diabetic male
metadata_d_nd_m <- filter(metadata, diabetic_status_and_sex == "Diabetic_male" | diabetic_status_and_sex == "Non-diabetic_female")

#### 5. diabetic male vs pre-diabetic male
metadata_d_pd_m <- filter(metadata, diabetic_status_and_sex == "Diabetic_male" | diabetic_status_and_sex == "Pre-diabetic_female")

#### 6. pre-diabetic male vs non-diabetic male
metadata_pd_nd_m <- filter(metadata, diabetic_status_and_sex == "Pre-diabetic_male" | diabetic_status_and_sex == "Non-diabetic_male")




#### 7. diabetic male vs diabetic female
metadata_d_f_m <- filter(metadata, diabetic_status_and_sex == "Diabetic_female" | diabetic_status_and_sex == "Diabetic_male")

#### 8. pre-diabetic male vs pre-diabetic female
metadata_pd_f_m <- filter(metadata, diabetic_status_and_sex == "Pre-diabetic_female" | diabetic_status_and_sex == "Pre-diabetic_male")

#### 9. non-diabetic male vs non-diabetic female
metadata_nd_f_m <- filter(metadata, diabetic_status_and_sex == "Non-diabetic_female" | diabetic_status_and_sex == "Non-diabetic_male")

