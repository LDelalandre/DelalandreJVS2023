source("scripts/Packages.R")

# NB: processed data is available on InDoRES.
# Consequently, jump directly to scripts performing analyses and producing figures,
# which start with the number 3


# Load and process raw data ####
rawdata <- F
if (rawdata == T){
  ## Load traits and include measurements made during my Ph.D. ####
  source("scripts/Data_traits.R") 
  
  ## Traits: Compute species mean trait values in the G+F and GUs conditions ####
  source("scripts/2_1_Process_data_mean_trait.R")
  
  ## Abundance: Export abundance data for the focal years ####
  source("scripts/2_2_Process_data_abundance.R")
  # Rq: sensitivity analysis of the results on the year kept in diachro to be performed.
  # (cf. also 2022_10_24_critiques et réponses.docx)
}



# Figures and analyses ####
# Analyses performed on the data available on InDORES

## Environment and Richness ####
source("scripts/3_Cover_annuals.R") # figure 2

## Trait coverage
source("scripts/3_Trait_coverage.R") # supplementary

## Trait differences
source("scripts/3_Trait_differences.R") # figure 3 and supplementary
source("scripts/3_PCA_final.R") # figure 4

## Intraspecific differences 
source("scripts/3_table_intrasp_differences.R") # table 2

## Comparison ITV and turnover
source("scripts/3_sp_replacement_ITV.R")

