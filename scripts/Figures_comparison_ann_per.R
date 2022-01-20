
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  # merge(name_LH,by="Code_Sp") %>% 
  # relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())





MEAN_CSR2 <- MEAN_CSR %>% 
  select(LifeHistory,C,S,R)




ggplot(,aes(x=))