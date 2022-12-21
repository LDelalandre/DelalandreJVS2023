source("scripts/Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)


MEAN_no_subset <- read.csv2("outputs/data/mean_attribute_per_treatment.csv",encoding = "latin1") %>%
  filter(!is.na(SLA)) %>%
  filter(!(LifeForm1 %in% c("DPh","EPh")))%>% 
  filter(!(species== "Geranium dissectum - pétiole")) 
# keep only traits measured in the Nat_Sab
# = compare trait values in Nat_Sab and in fertile
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_completed.csv")%>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole"))

code_sp_lifeform <- MEAN_no_subset %>% 
  select(code_sp,LifeForm1) %>% 
  unique() %>% 
  rename(Code_Sp = code_sp)

# Many species in the data do not appear in the PCA.
code_sp_lifeform %>% group_by(LifeForm1) %>% summarize(n=n())


ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")

# CWM at transect level ####
ab_traits_fer <- ab_fer %>% 
  left_join(MEAN %>% filter(treatment == "Fer"),
            by = c("species","code_sp","LifeForm1","treatment")) 

# I keep it with species level I case I want to compute distance to CWM for some species
CWM_sp_fer <- ab_traits_fer %>% 
  ungroup() %>% 
  group_by(id_transect_quadrat) %>% 
  # filter((LifeForm1 == "The")) %>%
  mutate(relat_ab2 = relat_ab/sum(relat_ab)) %>%  # useful when subsample just annuals or perennials
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab2,na.rm=T),
                         CWV = ~ weighted_var(.,relat_ab2))) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  rename_at( vars( contains( "_CWV") ), list( ~paste("CWV", gsub("_CWV", "", .), sep = "_") ) ) %>% 
  unique()


CWM_fer <- CWM_sp_fer %>% 
  select(paddock, id_transect_quadrat,starts_with("CWM"),starts_with("CWV")) %>% 
  unique()

# Il faudrait regarder la variance, la skewness, la kurtosis...

# " leaf Area"
CWM_fer %>% 
  ggplot(aes(x = CWM_L_Area)) +
  geom_histogram(binwidth = 50) +
  geom_vline(xintercept=300,color = "red") +
  ggtitle("annuals")
  xlim(c(0,300))

MEAN %>% 
  ggplot(aes(x = L_Area,color = LifeForm1)) +
  geom_histogram(binwidth = 50)+
  xlim(c(0,1000))

# SeedMass
CWM_fer %>% 
  ggplot(aes(x = CWM_SeedMass)) +
  geom_histogram() +
  geom_vline(xintercept=1.9,color = "red") +
  ggtitle("all")

