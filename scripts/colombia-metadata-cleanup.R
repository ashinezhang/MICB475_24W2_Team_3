library(tidyverse)
library(dplyr)

# Load in the datasets

col_meta_fp <- "./data/colombia_metadata.csv"
col_meta <- read.csv(file = col_meta_fp)

# Annotate diabetic status of individuals

col_meta_diabetes <- col_meta %>%
  mutate(Group = ifelse(glucose > 125 | Hemoglobin_a1c > 6.4, "Diabetic", 
                        ifelse(glucose < 100 & Hemoglobin_a1c < 5.7, "Non-diabetic", "Pre-diabetic")))

# Annotate diabetic and sex status of individuals
col_meta_diabetes_sex <- col_meta_diabetes %>%
  mutate(Group_sex = case_when(
    Group == "Diabetic" & sex == "male" ~ "Diabetic_male",
    Group == "Diabetic" & sex == "female" ~ "Diabetic_female",
    Group == "Pre-diabetic" & sex == "male" ~ "Pre-diabetic_male",
    Group == "Pre-diabetic" & sex == "female" ~ "Pre-diabetic_female",
    Group == "Non-diabetic" & sex == "male" ~ "Non-diabetic_male",
    Group == "Non-diabetic" & sex == "female" ~ "Non-diabetic_female"
  ))


# Select columns for merging and check Group counts

col_selected <- col_meta_diabetes_sex %>%
  select(c("X.SampleID", "age_years", "BMI", "country", "sex", "Group", "Group_sex"))

Group_counts <- col_selected%>%
  count(Group)

Group_sex_counts <- col_selected%>%
  count(Group_sex)

# Adjust Weight classification and column names

col_clean <- col_selected %>%
  mutate(BMI_class = ifelse(BMI <= 25, "lean", 
                            ifelse(BMI >25 & BMI <= 30, "overweight", 
                                   ifelse(BMI >30 & BMI <=35, "obese", "severe"))))
colnames(col_BMI) <- c("sample-id", "Age", "BMI", "Country", "Sex", "Group", "Group and Sex", "BMI_class")

# Make table 
write.table(col_clean, file = "./data/new_colombia_diabetes_metadata.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
                 