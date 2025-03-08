library(tidyverse)
library(dplyr)

# Load in the datasets

col_meta_fp <- "./data/colombia_metadata.csv"
col_meta <- read.csv(file = col_meta_fp)

mex_meta_fp <- "./data/diabetes_mexico_metadata.csv"
mex_meta <- read.csv(file = mex_meta_fp)



# Annotate diabetic status of individuals

col_meta_diabetes <- col_meta %>%
  mutate(Group = ifelse(glucose > 125 | Hemoglobin_a1c > 6.4, "Diabetic", 
                        ifelse(glucose < 100 & Hemoglobin_a1c < 5.7, "Non-diabetic", "Pre-diabetic")))

# Select columns for merging and check Group counts

col_selected <- col_meta_diabetes %>%
  select(c("X.SampleID", "age_years", "BMI", "country", "sex", "Group"))

Group_counts <- col_selected%>%
  count(Group)

##March 8th 2025 - This section was adjusted to meet requirments for qiime2. Namely changed Sample_ID column name to sample-id and generated a tsv instead of a csv
# Adjust Weight classification and column names

col_merg <- col_selected %>%
  mutate(BMI_class = ifelse(BMI <= 25, "lean", 
                            ifelse(BMI >25 & BMI <= 30, "overweight", 
                                   ifelse(BMI >30 & BMI <=35, "obese", "severe"))))
colnames(col_merg) <- c("sample-id", "Age", "BMI", "Country", "Sex", "Group", "BMI_class")

#Add country to Mexico dataset
mex_meta$Country <- "Mexico"

#Rename Mexico columns
mex_meta<-rename(mex_meta, Sample=Sample_ID,Sample_ID=sample.id, BMI_class=Obesity)

#Make sex lowercase to match with Colombia dataset
mex_meta$Sex <- tolower(mex_meta$Sex)


#Select columns and merge datasets
mex_selected <- mex_meta %>%
  select(c("sample-id", "Age", "BMI", "Country", "Sex", "Group", "BMI_class"))

diabetes_merged <- bind_rows(col_merg, mex_selected)

##old code: write.csv(diabetes_merged, "./data/merged_diabetes_metadata.csv", row.names = FALSE, col.names = TRUE)
##new code:
write.table(diabetes_merged, file = "./data/new_merged_diabetes_metadata.csv", sep = "\t", row.names = FALSE, col.names = TRUE)
