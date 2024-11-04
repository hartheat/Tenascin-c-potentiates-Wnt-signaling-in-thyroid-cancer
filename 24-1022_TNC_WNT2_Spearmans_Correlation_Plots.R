# Author: Matthew Aaron Loberg
# Date: 24-1022
# Script: 24-1022_TNC_WNT2_Spearmans_Correlation_Plots.R

# Goal: Make Spearmans correlation plots between WNT2 and TNC
# I will be making correlation plots between TNC and WNT2 across three cohorts: all malignant; WDTC; ATC
# For all of the cohorts, I will be EXCLUDING NIFTP

# Load necessary packages:
library(tidyverse)
library(cowplot)
library(RColorBrewer)

# Load cleaned merged data
cohort <- readRDS(file = "data_in_use/22-0606_CleanedMergedData_DESeq2NormalizedReads.rds")

#### Establish a column for labeling as Follicular, Papillary, or Transformed ####
cohort$Category <- "Benign"
# Anything NOT follicular, papillary, or transformed will have a category of benign

# Follicular lesions
for(i in 1:nrow(cohort)){
  if(cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "NIFTP" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "FC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "HC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "EFVPTC"){
        cohort$Category[i] <- "Follicular"
  }
}

# Papillary lesions
for(i in 1:nrow(cohort)){
  if(cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "IFVPTC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PTC"){
        cohort$Category[i] <- "Papillary"
  }
}

# Transformed lesion
for(i in 1:nrow(cohort)){
  if(cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PDTC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "ATC"){
        cohort$Category[i] <- "Transformed"
  }
}

#### Create a new variable that will be the following groups: Benign, RAS-like WDTC, BRAF-like WDTC, PDTC, ATC
cohort$Diagnosis_Simplified <- cohort$Diagnosis
cohort$BRS <- as.numeric(cohort$BRS)
for(i in 1:nrow(cohort)){
  if(cohort$Category[i] == "Papillary" | cohort$Category[i] == "Follicular"){
    if(cohort$BRS[i] >= 0){
      cohort$Diagnosis_Simplified[i] <- "RAS-Like WDTC"
    }
    else if(cohort$BRS[i] < 0){
      cohort$Diagnosis_Simplified[i] <- "BRAF-Like WDTC"
    }
  }
  if(cohort$Category[i] == "Benign"){
    cohort$Diagnosis_Simplified[i] <- "Benign"
  }
}
table(cohort$Diagnosis_Simplified)
table(cohort$Diagnosis)

########## Plots ###############
# I will be making correlation plots between TNC and WNT2 across three cohorts: all malignant; WDTC; ATC
# For all of the cohorts, I will be EXCLUDING NIFTP

# Subset to malignant cohort
malignant_cohort <- cohort %>% subset(Diagnosis_Simplified != "Benign" & Diagnosis != "NIFTP")
table(malignant_cohort$Diagnosis)

# malignant cohort TNC WNT2 correlation plot
malignant_correlation <- ggplot(malignant_cohort, 
       aes(log2(TNC+1), log2(WNT2+1))) + 
  geom_point(alpha = 0.7) + 
  theme_classic()
ggsave(file = "outputs/24-1022_WNT2_TNC_Correlation/24-1022_WNT2_TNC_Spearmans_Correlation_Malignant_Cohort.png", 
       malignant_correlation, dpi = 600, height = 4.5, width = 4.5)

# Stats
correlation <- cor(log2(malignant_cohort$TNC+1), log2(malignant_cohort$WNT2+1), method = "spearman")
print(correlation)
cor_test_result <- cor.test(x = log2(malignant_cohort$TNC+1), y = log2(malignant_cohort$WNT2+1), method = "spearman")
print(cor_test_result)





# Subset to ATCs
ATC_cohort <- malignant_cohort %>% subset(Diagnosis == "ATC")

# ATC cohort TNC WNT2 correlation plot
ATC_correlation <- ggplot(ATC_cohort, 
                                aes(log2(TNC+1), log2(WNT2+1))) + 
  geom_point(alpha = 0.7) + 
  theme_classic()
ggsave(file = "outputs/24-1022_WNT2_TNC_Correlation/24-1022_WNT2_TNC_Spearmans_Correlation_ATC_Cohort.png", 
       ATC_correlation, dpi = 600, height = 4.5, width = 4.5)

# Stats
correlation <- cor(log2(ATC_cohort$TNC+1), log2(ATC_cohort$WNT2+1), method = "spearman")
print(correlation)
cor_test_result <- cor.test(x = log2(ATC_cohort$TNC+1), y = log2(ATC_cohort$WNT2+1), method = "spearman")
print(cor_test_result)





# Subset to WDTCs
WDTC_cohort <- malignant_cohort %>% subset(Diagnosis != "PDTC" & Diagnosis != "ATC")

# WDTC cohort TNC WNT2 correlation plot
WDTC_correlation <- ggplot(WDTC_cohort, 
                          aes(log2(TNC+1), log2(WNT2+1))) + 
  geom_point(alpha = 0.7) + 
  theme_classic()
ggsave(file = "outputs/24-1022_WNT2_TNC_Correlation/24-1022_WNT2_TNC_Spearmans_Correlation_WDTC_Cohort.png", 
       WDTC_correlation, dpi = 600, height = 4.5, width = 4.5)

# Stats
correlation <- cor(log2(WDTC_cohort$TNC+1), log2(WDTC_cohort$WNT2+1), method = "spearman")
print(correlation)
cor_test_result <- cor.test(x = log2(WDTC_cohort$TNC+1), y = log2(WDTC_cohort$WNT2+1), method = "spearman")
print(cor_test_result)



##### CLEANING UP 
rm(list = ls())




