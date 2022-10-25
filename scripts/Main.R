source("scripts/Packages.R")



# Load and process raw data ####

## Load traits and include measurements made during my Ph.D. ####
source("scripts/Data_traits.R") 

## Traits: Compute species mean trait values in the G+F and GUs conditions ####
source("scripts/2_1_Process_data_mean_trait.R")

## Abundance: Export abundance data for the focal years ####
source("scripts/2_2_Process_data_abundance.R")
# Rq: sensitivity analysis of the results on the year kept in diachro to be performed.
# (cf. also 2022_10_24_critiques et réponses.docx)

## CWM: Compute CWM ofCSR scores ####
source("scripts/2_3_Process_data_CWM_CSR.R")


# Figures and analyses ####

## Environment and Richness ####
# "3_cover_annuals.R"

## PCA ####
# " 3_PCA.R"


## Trait coverage ####
