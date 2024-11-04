### Author: Matthew Aaron Loberg
### Date: June 5, 2024
### Script: 24-0601_CGRNA_Heather_Metastasis_Exploration.R

# Update 24-0601
# Inputting new clinical file (VUMC.cohort.MAL_5-31-24.csv)
# New clinical file contains data on whether a sample should be used for PFS
# For information on this new column (column title: PFS_Sample), please see the README file on the lab SharePoint
# README file location: Matt > Matt_Sequencing_Cohort_Clinical > README.docx
# Note: moved clinical_cohort_data files to a new "data_in_use" subfolder "clinical_cohort_data"

# Update 24-0605
# The update between this and 24-0601 is that I just did not finish working on this on 24-0601 (got distracted by ENDO2024 conference)
# I will be starting over today and cleaning the code up a bit, saving plots, etc.

# Note on cohort: 
# the subset of malignant samples that are eligible for PFS is in total 136; George had 138
# I need to compare and see which samples are different between us
# 24-0605 - STILL NEED TO DO THIS

# Update 24-1024:
# Restricting out "localdisease" that are not primary tumors 

### Load required packages
library(ggplot2)
library(survival)

library(tidyverse)
library(cowplot)

### Set WD (or open script as a project within defined directory)
# Set your working directory if needed (I do not as working within 2024_Bulk_RNA_Deconvolution_Analysis R Project)
# setwd("_") # Input in the "_" the working directory that you wish to assign

### Load Data
# Note, while doing this, I decided to switch to the most recent version of the clinical file
# My concern was that previous versions of the clinical file may be outdated with regard to some of the met categories
# I switched from "VUMC.cohort.GX_9-8-22_v2.csv" to the most recent in George's folder, which is
# I have now updated to my own most recent: 5-31-24
ClinicalData <- as_tibble(read_csv("VUMCcohort.csv")) # All Data

# Read in normalized RNA seq counts
RNA <- as_tibble(read.table(file = "21-1214_DESeq2_Normalized_RNA_Counts.txt", sep = '\t', header = TRUE))
RNA$RNA.ID <- sub('.', '', RNA$RNA.ID) # Remove the "X" from the start of "RNA.ID" so that it is compatible

# Merge normalized counts and clinical data
cohort <- merge(x = ClinicalData, y = RNA, by.x = "RNA.ID", by.y = "RNA.ID")
rm(ClinicalData, RNA)

# Change "Distant.met.patient" to "YES" or "NO"
for(i in 1:nrow(cohort)){
  if(!is.na(cohort$Distant.met.patient[i])){
    if(cohort$Distant.met.patient[i] == 0){
      cohort$Distant.met.patient[i] <- "NO"
    }
    else if(cohort$Distant.met.patient[i] == 1){
      cohort$Distant.met.patient[i] <- "YES"
    }
  }
  if(!is.na(cohort$Local.met.patient[i])){
    if(cohort$Local.met.patient[i] == "0"){
      cohort$Local.met.patient[i] <- "NO"
    }
    else if(cohort$Local.met.patient[i] == "1"){
      cohort$Local.met.patient[i] <- "YES"
    }
  }
  
}

####### 24-1024 UPDATE #######
# Removing local disease samples that are NOT actually primary
cohort <- cohort %>% subset(Location.type.detailed != "Localdisease" | PFS_Sample != "NO")


# Make a function for all of the plots that I want to make
MetPlots <- function(cohort, output_folder){
    
    # Distant met yes/no for PRIMARIES only
    plot <- ggplot(cohort %>% subset(Location.type == "Primary" & !is.na(Distant.met.patient)), aes(Distant.met.patient, log2(TNC+1))) +
      geom_boxplot(outlier.size = -1, 
                   aes(fill = Distant.met.patient),
                   alpha = 0.4, 
                   show.legend = FALSE) + 
      geom_jitter(aes(),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 1.5, 
                  alpha = 0.7,
                  show.legend = FALSE) +
      scale_fill_manual(values = c("#df3e7d", "#673ab7")) +
      labs (x ="Met patient", y = "log2(TNC)") + 
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 14.5), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 16)) +
      scale_x_discrete(name ="Distant Met Patient", limits = c("NO", "YES")) +
      scale_y_continuous(breaks = c(0,5,10,15),
                         limits = c(0,15.1))
    ggsave(filename = file.path(output_folder, "Distant_Met_Patient_LocalDisease_Only.png"),
           plot, height = 3, width = 3, dpi = 600)
    
    test <- cohort %>% subset(Location.type == "Primary" & !is.na(Distant.met.patient))
    kruskal <- kruskal.test(log2(TNC+1) ~ Distant.met.patient, data = test)
    pairwise <- pairwise.wilcox.test(log2(test$TNC+1), test$Distant.met.patient, 
                                     p.adjust.method = "bonferroni")
    rm(test)
    
    # Local met yes/no for PRIMARIES only
    plot <- ggplot(cohort %>% subset(Location.type == "Primary"), aes(Local.met.patient, log2(TNC+1))) +
      geom_boxplot(outlier.size = -1, 
                   aes(fill = Local.met.patient),
                   alpha = 0.4, 
                   show.legend = FALSE) + 
      geom_jitter(aes(),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 1.5, 
                  alpha = 0.7,
                  show.legend = FALSE) +
      scale_fill_manual(values = c("#f886bb", "#370140")) +
      labs (x ="Met patient", y = "log2(TNC)") + 
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 14.5), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 16)) +
      scale_x_discrete(name ="Local Met Patient", limits = c("NO", "YES")) +
      scale_y_continuous(breaks = c(0,5,10,15),
                         limits = c(0,15.1))
    ggsave(filename = file.path(output_folder, "Local_Met_Patient_LocalDisease_Only.png"),
           plot, height = 3, width = 3, dpi = 600)
    test <- cohort %>% subset(Location.type == "Primary")
    kruskal <- kruskal %>% append(kruskal.test(log2(TNC+1) ~ Local.met.patient, data = test)) 
    pairwise <- pairwise %>% append(pairwise.wilcox.test(log2(test$TNC+1), test$Local.met.patient, 
                                                         p.adjust.method = "bonferroni"))
    rm(test)
    
    # Samples by Location.type (for primary and LN met only)
    plot <- ggplot(cohort %>% subset(Location.type == "Primary" | Location.type == "Localmet"), aes(Location.type, log2(TNC+1))) +
                     geom_boxplot(outlier.size = -1, 
                                  aes(fill = Location.type),
                                  alpha = 0.4, 
                                  show.legend = FALSE) + 
                     geom_jitter(aes(),
                                 position = position_jitter(width = 0.1, height = 0),
                                 size = 1.5, 
                                 alpha = 0.7,
                                 show.legend = FALSE) +
                     scale_fill_manual(values = c("#a84c97", "#ff474c")) + # can modify here to define colors
                     labs (x ="Location Simplified", y = "log2(TNC)") + 
                     theme_classic() + 
                     theme(
                       plot.title = element_text(face = "bold", hjust = 0.5), 
                       axis.title.x = element_text(face = "bold", size = 14.5), 
                       axis.title.y = element_text(face = "bold", size = 12),
                       axis.text.x = element_text(face = "bold", size = 10.5),
                       axis.text.y = element_text(face = "bold", size = 16)) +
                     scale_x_discrete(name ="Location Simplified ", limits = c("Primary", "Localmet")) +
                     scale_y_continuous(breaks = c(0,5,10,15),
                                        limits = c(0,15.1))
    ggsave(filename = file.path(output_folder, "Location_Simplified.png"),
           plot, height = 3, width = 3, dpi = 600)
    test <- cohort %>% subset(Location.type == "Primary" | Location.type == "Localmet")
    kruskal <- kruskal %>% append(kruskal.test(log2(TNC+1) ~ Location.type, data = test)) 
    pairwise <- pairwise %>% append(pairwise.wilcox.test(log2(test$TNC+1), test$Location.type, 
                                                         p.adjust.method = "bonferroni"))
    
    # Return(NULL)
    return(list(kruskal, pairwise))
}

# Create an output directory: 
dir.create("outputs/24-1024_CGRNA_Met_Plots")

# MetPlots malignant cohort (NIFTP Excluded)
output_folder <- "outputs/24-1024_CGRNA_Met_Plots/Malignant_NIFTP_Excluded/"
dir.create(output_folder)
stats <- MetPlots(cohort = cohort %>% subset(Diagnosis != "MNG" & Diagnosis != "FA" &
                                             Diagnosis != "HA" & Diagnosis != "HT" & Diagnosis != "NIFTP"),
                            output_folder = output_folder)
saveRDS(stats, file.path(output_folder, "24-1024_Malignant_Samples_NIFTP_Excluded_Stats_List.RDS"))


# MetPlots WDTC Cohort (NIFTP Excluded)
output_folder <- "outputs/24-1024_CGRNA_Met_Plots/WDTC_NIFTP_Excluded/"
dir.create(output_folder)
stats <- MetPlots(cohort = cohort %>% subset(Diagnosis != "MNG" & Diagnosis != "FA" &
                                             Diagnosis != "HA" & Diagnosis != "HT" & 
                                             Diagnosis != "PDTC" & Diagnosis != "ATC" & Diagnosis != "NIFTP"),
                  output_folder = output_folder)
saveRDS(stats, file.path(output_folder, "24-1024_WDTC_NIFTP_Excluded_Stats_List.RDS"))

##### Cleaning UP ######
rm(list = ls())