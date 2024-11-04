# Author: Matthew Loberg
# Date: 24-1022
# Script: 24-1022_CGRNA_Heather_TNC_Expression.R

# Goal: Look at TNC expression by BRAF-like vs RAS-like within the well-differentiated cohort
# THIS TIME I WILL EXCLUDE NIFTP (different vs prior)

# Load necessary packages:
library(tidyverse)
library(cowplot)
library(RColorBrewer)

####### Loading and Setting Up Data ###################
ClinicalData <- as_tibble(read_csv("VUMCcohort.csv")) # All Data

# Read in normalized RNA seq counts
RNA <- as_tibble(read.table(file = "21-1214_DESeq2_Normalized_RNA_Counts.txt", sep = '\t', header = TRUE))
RNA$RNA.ID <- sub('.', '', RNA$RNA.ID) # Remove the "X" from the start of "RNA.ID" so that it is compatible

# Merge normalized counts and clinical data
cohort <- merge(x = ClinicalData, y = RNA, by.x = "RNA.ID", by.y = "RNA.ID")
rm(ClinicalData, RNA)

# Make a BRAF-RAS designation based on the BRAF-RAS score
cohort$BRS <- as.numeric(cohort$BRS)
cohort$BRS_Designation <- "NONE"
for(i in 1:nrow(cohort)){
  if(cohort$BRS[i] > 0){
    cohort$BRS_Designation[i] <- "RAS-LIKE"
  }
  else if(cohort$BRS[i] < 0){
    cohort$BRS_Designation[i] <- "BRAF-LIKE"
  }
}
table(cohort$BRS_Designation)

########## Plots ###############

### Today plots
# WDTC TNC by BRAF-vs RAS

# First subset to WDTC
WDTC_Local_cohort <- cohort %>% subset(Location.type == "Primary" & 
                                       Diagnosis != "MNG" & 
                                       Diagnosis != "HA" &
                                       Diagnosis != "FA" & 
                                       Diagnosis != "HT" & 
                                       Diagnosis != "NIFTP" & 
                                       Diagnosis != "ATC" & 
                                       Diagnosis != "PDTC")

# now find max and min TNC values
max(log2(WDTC_Local_cohort$TNC+1)) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) (13.2)
min(log2(WDTC_Local_cohort$TNC+1)) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  (6.4)


# WDTC BRAF-like vs RAS-like
plot <- ggplot(WDTC_Local_cohort, aes(BRS_Designation, log2(TNC+1))) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRS_Designation),
               alpha = 0.4, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("lightblue", "#0099ff")) +
  labs (x = "Diagnosis", y = "log2(TNC)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS Category", limits = c("RAS-LIKE", "BRAF-LIKE")) +
  scale_y_continuous(breaks = c(6, 8, 10, 12, 14),
                     limits = c(6, 14)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1022_CGRNA_Met_Plots/24-1022_Thyroid_WDTC_BRS_Category_TNC_Boxplot.png", 
       width = 3, height = 3, 
       plot, dpi = 600)
kruskal.test(log2(TNC+1) ~ BRS_Designation, data = WDTC_Local_cohort) 
pairwise.wilcox.test(log2(WDTC_Local_cohort$TNC+1), WDTC_Local_cohort$BRS_Designation, 
                                                     p.adjust.method = "bonferroni")


##### CLEANING UP #####
rm(list = ls())

