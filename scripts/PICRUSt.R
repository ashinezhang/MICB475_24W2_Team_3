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
# this can take very long to load
# when asked if you want to update, type "n" for none to make it faster

# After installing all of its above dependencies, install ggpicrust2
install.packages("ggpicrust2")

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

