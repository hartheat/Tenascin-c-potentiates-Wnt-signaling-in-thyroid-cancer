# Author: Matthew Aaron Loberg
# Date: October 22, 2024
# Script: 24-1022_TCGA_TNC_Expression_BRAF_RAS.R

# Goal: Look at log2(TNC) expression in TCGA by BRAF-like or RAS-like designation

##### Load packages #####
library(tidyverse)

##### Load TCGA RSEM data and format for merge with meta data #####
# read in file 
TCGA_RSEM <- read.table(file = 'data_in_use/data_RNA_Seq_v2_expression_median.txt', sep = '\t', header = TRUE)
Gene_Symbols <- TCGA_RSEM$Hugo_Symbol # save the gene names columns as the rownames
TCGA_RSEM_matrix <- as.matrix(TCGA_RSEM[,3:ncol(TCGA_RSEM)]) # make into matrix with just the numerical data
TCGA_Transposed <- t(TCGA_RSEM_matrix)
SAMPLE_ID <- rownames(TCGA_Transposed)
TCGA_Transposed <- as_tibble(TCGA_Transposed)
colnames(TCGA_Transposed) <- Gene_Symbols # Make Gene Symbols the column names
TCGA_Transposed$SAMPLE_ID <- SAMPLE_ID # Set sample IDs
TCGA_Transposed$SAMPLE_ID <- gsub("\\.", "-", TCGA_Transposed$SAMPLE_ID)

######### Read in TCGA meta data  and format ##########
# Load Data
TCGA_Patient_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE) # Load patient data file
TCGA_Patient_Data <- as_tibble(TCGA_Patient_Data) # Make patient data into tibble for easier data cleaning

TCGA_Sample_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE) # Load sample data file
TCGA_Sample_Data <- as_tibble(TCGA_Sample_Data) # Make sample data into tibble for easier data cleaning


# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)
# TCGA_Sample_Data_Restricted$PATIENT_ID <- substr(TCGA_Sample_Data_Restricted$PATIENT_ID, 1, nchar(TCGA_Sample_Data_Restricted$PATIENT_ID)-3) # Old line of code from when I had a restricted file
# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

####### Merge RSEM and patient sample combined ##############
TCGA_merged <- merge(Patient_Sample_Combined, TCGA_Transposed)
rm(TCGA_RSEM, TCGA_RSEM_matrix, TCGA_Transposed, SAMPLE_ID, Gene_Symbols, Patient_Sample_Combined, TCGA_Patient_Data, TCGA_Sample_Data)

##### Plot log2 of TNC by BRAF-like RAS-like designation #########

# TNC - compare BRAF-like vs RAS-like
max(log2(TCGA_merged$TNC+1)) # shows max of 15.64
min(log2(TCGA_merged$TNC+1)) # shows min of 3.7
plot <- ggplot(TCGA_merged %>% subset(BRAFV600E_RAS != ""), aes(BRAFV600E_RAS, log2(TNC+1))) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRAF-RAS Phenotype", y = "log2(TNC+1)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15),
                     limits = c(0, 15.8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1022_TCGA_TNC_by_BRAF_RAS_Phenotype/24-1022_TNC_BRAF_RAS.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

stats_cohort <- TCGA_merged %>% subset(BRAFV600E_RAS != "")
kruskal.test(log2(TNC+1) ~ BRAFV600E_RAS, 
             data = stats_cohort) # p-value < 2.2e-16
pairwise.wilcox.test(stats_cohort$TNC, 
                     stats_cohort$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # <2e-16

## Can test other genes of interest as desired 

##### CLEANING UP
rm(list = ls())
