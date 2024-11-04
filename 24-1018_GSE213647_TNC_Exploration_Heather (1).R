# Author: Matthew Aaron Loberg
# Date: October 18, 2024
# Script: 24-1018_GSE213647_TNC_Exploration_Heather.R

# 24-1018 Update
# Simplified version w/ just plots of interest for Heather

# Goal: Look at TNC expression within GSE213647 from Lee et al. 2024 Nature Communications

# This bulk RNA sequencing data from Lee et al. was TPM normalized and subset to protein coding genes on 24-0910 by Hua-Chang
# See email from Hua-Chang below: 

# Hua-Change email:
# Hi Matt,
# 
# Sure, I have created a protein-coding subset, which is saved at the following location: /data/h_vivian_weiss/20240905_GSE213647_thyroid_cancer/GSE213647_thyroid_cancer.TPM.protein_coding.tsv. This file contains 19,161 protein-coding genes.
# 
# Additionally, I have updated the full TPM table by adding a new "gene_type" column. The gene type information was sourced from the latest HGNC database.
# 
# Please let me know if you need any further details.
# 
# Regards,
# HuaChang

# I then downloaded Hua-Chang's TPM file 
# It is within the 2024_Bulk_RNA_Deconvolution_Analysis "data_in_use" folder saved as follows:
# "GSE213647_thyroid_cancer.TPM.protein_coding.tsv"

##### Load packages #####
library(tidyverse)

##### Load TPM data #####
# read in file + format with genes as columns, including a GEO_ID column as well for sample identification
tpm <- read.table(file = 'data_in_use/GSE213647_thyroid_cancer.TPM.protein_coding.tsv', sep = '\t', header = TRUE)
tpm_t <- t(tpm)
colnames(tpm_t) <- tpm_t[2,]
tpm_t <- tpm_t[5:nrow(tpm_t),]
tpm_t <- as.data.frame(tpm_t) # matrix format was causing issues with adding SampleID variable so I am making it into a data frame
tpm_t$GEO_ID <- rownames(tpm_t)
tpm_t$GEO_ID <- substr(tpm_t$GEO_ID, 1, 10)

# Read in meta data
GSE213647_meta <- read.csv(file = "data_in_use/GSE213647_Meta_Data.csv")

# Merge meta and GSVA
GSE213647_Merged <- merge(GSE213647_meta, tpm_t)

# Cleaning up
rm(tpm, tpm_t, GSE213647_meta)

# Check # of diagnoses 
table(GSE213647_Merged$Diagnosis)


### Plot TNC by Diagnosis

# TNC - all diagnoses
GSE213647_Merged$TNC <- as.numeric(GSE213647_Merged$TNC)
max(log2(GSE213647_Merged$TNC+1)) 
min(log2(GSE213647_Merged$TNC+1)) 
plot <- ggplot(GSE213647_Merged, aes(Diagnosis, log2(TNC+1))) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  #scale_fill_manual(values = c("grey", "lightblue", "red", "purple")) + # INPUT DESIRED COLORS HERE
  labs (x = "Diagnosis", y = "log2(TNC TPM)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     limits = c(0, 8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1003_GSE213647_TNC_Exploration/24-1003_All_Diagnoses_log2_TNC_TPM.png", 
       width = 5, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(TNC ~ Diagnosis, 
             data = GSE213647_Merged) 
pairwise.wilcox.test(GSE213647_Merged$TNC, 
                     GSE213647_Merged$Diagnosis, 
                     p.adjust.method = "bonferroni") 


### TNC PDTC Excluded
plot <- ggplot(GSE213647_Merged %>% subset(Diagnosis != "PDTC"), aes(Diagnosis, log2(TNC+1))) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  #scale_fill_manual(values = c("grey", "lightblue", "red")) + # INPUT DESIRED COLORS HERE
  labs (x = "Diagnosis", y = "log2(TNC TPM)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "ATC")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     limits = c(0, 8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1003_GSE213647_TNC_Exploration/24-1003_All_Diagnoses_log2_TNC_TPM_PDTC_Excluded.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)


### Check STATS w/ non-parametric test
temp <- GSE213647_Merged %>% subset(Diagnosis != "PDTC")
kruskal.test(TNC ~ Diagnosis, 
             data = temp) 
pairwise.wilcox.test(temp$TNC, 
                     temp$Diagnosis, 
                     p.adjust.method = "bonferroni") 
rm(temp)













