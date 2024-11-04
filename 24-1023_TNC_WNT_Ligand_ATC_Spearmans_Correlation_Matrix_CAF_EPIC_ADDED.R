# Author: Matthew Aaron Loberg
# Date: 24-1022
# Script: 24-1022_TNC_WNT_Ligand_ATC_Spearmans_Correlation_Matrix.R

# goal: look at spearmans correlation of Wnt ligands with TNC and plot in a correlation matrix
# Do this for ligands that Diana found were up in ATC in her paper (WNT1, WNT2, WNT5B, WNT6, WNT7A, WNT10A, WNT10B)

# 24-1023 Update 
# Adding CAF EPIC

# Load necessary packages:
library(tidyverse)
library(RColorBrewer)
library("corrplot")


##### Load cohort #####
cohort <- readRDS(file = "data_in_use/22-0606_CleanedMergedData_DESeq2NormalizedReads.rds")

# Make CAF EPIC a numeric
cohort$CAFs <- as.numeric(cohort$`Cancer associated fibroblast_EPIC`)

# Make Matrix of Wnt ligand correlations with TNC using Spearman's correlation
ATC <- cohort %>% subset(Diagnosis == "ATC")
correlation_data <- ATC[c("TNC", "WNT1", "WNT2", "WNT5B", "WNT6", "WNT7A", "WNT10A", "WNT10B", "CAFs")] # WNT ligands increased in ATC + TNC + CAF EPIC
corr_mat=cor(correlation_data,method="s") #create Spearman correlation matrix
#col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444")) # alt color option
#col <- colorRampPalette(c("red", "white", "blue")) # alt color option

png(file = "outputs/24-1023_WNT2_TNC_Correlation/24-1023_ATC_Wnt_Ligand_CorrPlot.png", res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black",
         tl.col = "black")

dev.off()

# alt colors
png(file = "outputs/24-1023_WNT2_TNC_Correlation/24-1023_ATC_Wnt_Ligand_CorrPlot_alt_Colors.png", res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black",
         tl.col = "black",
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it

dev.off()
