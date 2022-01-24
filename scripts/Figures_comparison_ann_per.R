source("scripts/1. Packages.R")

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  # merge(name_LH,by="Code_Sp") %>% 
  # relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())

# MEAN_CSR2 <- MEAN_CSR %>% 
#   select(LifeHistory,C,S,R)


# Avoir autant de fois l'espèce qu'elle apparait dans les relevés (mauvaise option)
ab_traits_fer_CSR <- ab_fer %>% 
  left_join(MEAN_CSR %>% filter(treatment == "Fer"),
            by = c("species","code_sp","LifeForm1","treatment")) %>% 
    select(code_sp,LifeHistory,C,S,R) %>% 
  filter(!is.na(C))

ab_traits_nat_CSR <- ab_nat %>% 
  left_join(MEAN_CSR %>% filter(treatment == "Nat"),
            by = c("species","code_sp","LifeHistory")) %>% 
  select(depth,code_sp,LifeHistory,C,S,R)%>% 
  filter(!is.na(C))

# Option 2 : prendre les traits des espèces mesurés dans l'un ou l'autre traitement
MEAN_CSR %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  geom_boxplot() +
  facet_wrap(~treatment)


# Option 3 : sélectionner les espèces qui apparaissent dans les relevés de maud superficiel et diachro

# dans le fertile
species_fer <- ab_fer %>% 
  pull(code_sp) %>% 
  unique()

MEAN_CSR %>% 
  filter(code_sp %in% species_fer) %>% 
  select(code_sp,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory)) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  geom_boxplot() +
  ggtitle("Fer") 

# Dans le natif superficiel
species_S <- ab_nat %>% 
  filter(depth == "S") %>% 
  pull(code_sp) %>% 
  unique()

MEAN_CSR %>% 
  filter(code_sp %in% species_S) %>% 
  select(code_sp,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory)) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
    geom_boxplot() +
  ggtitle("Nat superficiel") 


# Dans le natif
species_nat <- ab_nat %>% 
  pull(code_sp) %>% 
  unique()

MEAN_CSR %>% 
  filter(code_sp %in% species_nat) %>% 
  select(code_sp,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory)) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  geom_boxplot() +
  ggtitle("Nat") 



