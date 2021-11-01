## Specify the packages of interest ####

packages = c("tidyverse", 
             "FactoMineR", "factoextra", # PCA computation
             "funrar", # Functional distinctiveness computation
             # "cowplot", "gridExtra", "ggsignif", "egg", "ggpubr","ggcorrplot", # Packages used for plots
             # "kableExtra", # Represent tables
             "openxlsx", # to open excel files
             "TAM" # to compute weighted CWV, skewness, and kurtosis
)

## Load, or install and load, packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    library(x, character.only = TRUE)
  }
)


# Complete folder structure ####

# Create folders to export processed data, plots, tables, and figures, if they do not exist #
if(file.exists(file.path("outputs"))==F){
  dir.create(file.path("outputs"))  
}

if(file.exists(file.path("outputs/data"))==F){
  dir.create(file.path("outputs/data"))  
}

if(file.exists(file.path("outputs/figures"))==F){
  dir.create(file.path("outputs/figures"))  
}

if(file.exists(file.path("outputs/plots"))==F){
  dir.create(file.path("outputs/plots"))  
}

if(file.exists(file.path("outputs/tables"))==F){
  dir.create(file.path("outputs/tables"))  
}
