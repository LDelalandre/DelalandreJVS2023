source("scripts/1. Packages.R")
source("scripts/2. Import files.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer

MEAN_treatment <- MEAN %>% filter(Trtmt == "Nat")

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))

name_LH <- LeafMorpho %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  select(Species,Code_Sp,LifeHistory) %>% 
  unique()

CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  merge(name_LH %>% select(-Code_Sp),by = "Species") %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())


trait <- "SeedMass"
ggplot(MEAN_treatment,aes_string(x=trait))+
  geom_density() +
  geom_point(data = MEAN %>% filter(LifeHistory=="annual"), aes_string(x=trait,y=0),color = "red")


# With abundance ####
ab_maud <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", sheet = "abondances par parcelle", 
                     startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(Species,abundance,-plot) %>%
  mutate(Species=str_replace(Species,"_"," ")) %>% 
  mutate(depth = str_sub(plot,start = 1L,end=1L)) %>% 
  mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
  select(-plot) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(abundance >0 ) %>% 
  merge(name_LH, by="Species") %>% 
  filter(!(Code_Sp == "STRIFSCAB"))

Depth = "S"
ab_maud2 <- ab_maud %>% 
  select(-paddock) %>% 
  group_by(Species,Code_Sp,depth,LifeHistory) %>% 
  summarize(abundance = sum(abundance)) %>% 
  # filter(depth == Depth) %>% 
  arrange(-abundance) 

traits_ab_maud <- merge(ab_maud2 %>% select(-c(Code_Sp,LifeHistory)),MEAN_treatment,"Species") %>% 
  filter(!Code_Sp.x == "STRIFSCAB")

# distribution traits 
traits_ab_maud %>% 
  filter(depth == "S") %>% 
  filter(LifeHistory=="annual") %>%
  ggplot(aes(x=SLA)) +
  geom_density()

# Link to temporal variation Adeline ####
ann_nat_sup <- ab_maud2 %>% 
  filter(LifeHistory == "annual" ) %>%
  filter(depth == "S") %>% 
  pull(Species)

tempo_var_adeline <- read.csv2("data/abundance/Decreasers increasers Adeline/temporal_variation_adeline.CSV") %>% 
  rename(Species = species)

tempo_nat_sup <- tempo_var_adeline %>% 
  filter(Species %in% ann_nat_sup)

tempo_traits <- left_join(tempo_nat_sup,
                          traits_ab_maud,
                          by = "Species")
ggplot(tempo_traits %>% filter(depth == "S"),aes(x=change,y=abundance,label = Species))+
  # ggrepel::geom_label_repel()
  geom_label()



# Distributions of trait values ####
trait <- "LNC"
col <- "LifeForm1.x"
# col <- "LifeHistory"

# ggplot(traits_ab_maud%>% 
#          filter(LifeHistory == "perennial"),
#        aes_string(x=trait))+
#   geom_density(color="blue") +
#   geom_density(data = traits_ab_maud %>% 
#                  filter(LifeHistory == "annual"),
#                aes_string(x=trait),color = "red") 

ggplot(traits_ab_maud, aes_string(x=trait,color = col))+
  geom_density()

traits_ab_maud %>% 
  filter(LifeForm1.x == "Hem" & SLA > 23)


# Which traits are missing
traits_ab_maud %>% 
  filter(LifeHistory=="annual") %>% 
  filter(is.na(SeedMass)) %>% 
  pull(Species)




# score CSR at species level
CSR_nat_Maud <- CSR %>% filter(Trtmt == "Nat") %>% 
  select(Species,C,S,R) %>% 
  merge(traits_ab_maud, by  = "Species")
  
score <- "C"
ggplot(CSR_nat_Maud%>% 
         filter(LifeHistory == "perennial"),
       aes_string(x=score)) +
  geom_density(color="blue") +
  geom_density(data = CSR_nat_Maud%>% 
                 filter(LifeHistory == "annual"),
               aes_string(x=score),color = "red") 
  

ggplot(traits_ab_maud,
       aes_string(x=trait))+
  geom_density() +
  geom_point(data = traits_ab_maud %>% filter(LifeHistory=="annual"), aes_string(x=trait,y=0),color = "red")


ggplot(traits_ab_maud %>% 
         filter(LifeHistory == "perennial"),
       aes_string(x=trait))+
  geom_density() +
  geom_point(data = traits_ab_maud %>% filter(LifeHistory=="annual"), aes_string(x=trait,y=0),color = "red")


# Are some annuals more uncoupled than others?
traits_ab_maud %>% 
  filter(LifeHistory == "perennial")


# Iris fertile ####
ab_iris <- read.xlsx("data/abundance/iris_Relevés bota.xlsx",sheet="Tout compilé") %>% # NB : je peux aussi importer ses données du natif (et du)
  dplyr::rename(abundance = 'Abondance.(tot)',species = Espece) %>% 
  mutate(plot = paste(Parc,Cage,sep="")) %>% 
  mutate(treatment = case_when(Traitement == "P+F+" ~ "Fer",
                               Traitement == "P+F-" ~ "Nat",
                               Traitement == "P-F-" ~ "Tem")) %>% 
  filter(treatment == "Fer") %>% 
  mutate(depth = Cage) %>% 
  select(-c('Abondance.(1)', 'Abondance.(2)','Présence.(1)', 'Présence.(2)', 'Présence.(tot)'))

ab_iris2 <- ab_iris %>% 
  group_by(species) %>% 
  summarize(abundance = sum(abundance)) %>% 
  arrange(species) %>% 
  dplyr::rename(Species = species) %>% 
  merge(names_LH,by = "Species") %>% 
  arrange(-abundance) 

traits_ab_iris <- merge(ab_iris2 %>% select(-c(Code_Sp,LifeHistory)),MEAN %>% filter(Trtmt == "Fer"),"Species")

ggplot(traits_ab_iris %>% 
         filter(LifeHistory == "perennial"),
       aes_string(x=trait))+
  geom_density(color="blue") +
  geom_density(data = traits_ab_iris%>% 
                 filter(LifeHistory == "annual"),
               aes_string(x=trait),color = "red")

traits_ab_iris %>% 
  filter(LifeHistory=="annual")

ggplot(MEAN %>% 
         filter(LifeHistory == "perennial"),
       aes_string(x=trait))+
  geom_density(color="blue") +
  geom_density(data = MEAN%>% 
                 filter(LifeHistory == "annual"),
               aes_string(x=trait),color = "red")



# score CSR at species level
CSR_fer_Iris <- CSR %>% filter(Trtmt == "Fer") %>% 
  select(Species,C,S,R) %>% 
  merge(traits_ab_iris, by  = "Species")

score <- "R"
ggplot(CSR_fer_Iris%>% 
         filter(LifeHistory == "perennial"),
       aes_string(x=score)) +
  geom_density(color="blue") +
  geom_density(data = CSR_nat_Maud%>% 
                 filter(LifeHistory == "annual"),
               aes_string(x=score),color = "red") 

