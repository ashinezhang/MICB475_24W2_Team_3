library(tidyverse)

# Load in the dataset

col_meta_fp <- "../data/colombia_metadata.csv"
col_meta <- read.csv(file = col_meta_fp)

# Annotate diabetic status of individuals

col_meta_diabetes <- col_meta %>%
  mutate(Group = ifelse(glucose > 125 | Hemoglobin_a1c > 6.4, "Diabetic", 
                          ifelse(glucose < 100 & Hemoglobin_a1c < 5.7, "Non-diabetic", "Pre-diabetic")))

# Select columns for merging and check Group counts

col_selected <- col_meta_diabetes %>%
  select(c("X.SampleID", "age_years", "BMI", "country", "Group"))

Group_counts <- diabetes_selected %>%
  count(Group)

# Adjust Weight classification and column names

col_merg <- col_selected %>%
  mutate(BMI_class = ifelse(BMI <= 25, "Lean", 
                           ifelse(BMI >25 & BMI <= 30, "Overweight", 
                                  ifelse(BMI >30 & BMI <=35, "Obese", "Severe"))))
colnames(col_merg) <- c("Sample_ID", "Age", "BMI", "Country", "Group", "BMI_class")

write.csv(col_merg, "../data/clean_columbia_metadata.csv", row.names = FALSE, col.names = TRUE)