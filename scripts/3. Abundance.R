source("scripts/1. Packages.R")

# Gets abundance and adds trait data

# Fertilized treatment ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))

Diachro_fer_2005 <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year == 2005) %>%  # /!\ No measure available after 2005! Should I keep 2006, or several years?
  filter(treatment == "Fer")

# Unfertilized treatment #### 
# (Maud's relevés)
soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))

ab_Maud_traits <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", sheet = "abondances par parcelle", 
                            startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(Species,abundance,-plot) %>%
  mutate(Species=str_replace(Species,"_"," ")) %>% 
  full_join(MEAN,by="Species") %>% 
  dplyr::rename(species=Species,code_sp=Code_Sp) %>% 
  filter(!is.na(abundance)) %>% 
  mutate(depth = str_sub(plot,start = 1L,end=1L)) %>% 
  mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
  select(-plot) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(Trtmt == "Nat") %>% 
  filter(abundance >0 )
write.csv2(ab_Maud_traits,"outputs/data/ab_Maud_traits.csv",row.names = F)

ab_Maud_traits %>% 
  select(code_sp,SLA,LDMC,L_Area) %>% 
  unique() %>% 
  filter(!(is.na(SLA))) %>%
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% 
  write.csv2("outputs/data/Pierce CSR/Traits_Maud.csv",row.names=F)