library(tidyverse)

disturbance <- read.table("data/Disturbance_DivHerbe.txt",header=T,sep="\t",dec=",")
disturbance

ggplot(disturbance, aes(x=Trtmt,y=Tx_CalcPic))+
  geom_boxplot()

