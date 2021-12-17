library(tidyverse)

decinc <- read.table("data/abundance/Decreasers increasers Garnier et al 2018/Supplementary Garnier 2018.txt",sep="\t",header=T)
andecinc <- decinc %>% 
  filter(Life.cycle == "Annual" )
