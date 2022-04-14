source("scripts/1. Packages.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
Annuals <- MEAN %>% 
  filter(LifeForm1=="The") %>% 
  pull(species) %>% 
  unique()


TRY.fage <- readRDS("data/traits/TRY.fage.Rds")


TRY.fage.ellenberg <- TRY.fage %>% 
  filter(grepl("Species environmental indicator value" ,TraitName))

TRY.fage.ellenberg[c(26,31),]
# Ici il y a un problème: c'est la même espèce (AccSpeciesName, nom consolidé, mais pas 
# le même nom d'espèce initial: Aphanes arvensis contre Aphanes arvensis agg., 
# et l'indice n'a pas la même valeur. C'est un agrégat = 
# a species group, i. e., closely related species that are often difficult to distinguish)
# Je retire, en première approx, les agrégats.
TRY.fage.ellenberg2 <- TRY.fage.ellenberg %>% 
  select(AccSpeciesName,SpeciesName,AccSpeciesID,TraitName,OrigValueStr,Dataset) %>% 
  filter(!grepl("agg." ,SpeciesName))

TRY.fage.ellenberg3 <- TRY.fage.ellenberg2 %>% 
  mutate(Trait = str_sub(TraitName, start = 63L, end = -1L)) %>% 
  select(-TraitName) %>% 
  unique() %>% 
  mutate(OrigValueStr = as.numeric(OrigValueStr)) %>% 
  spread(key = Trait, value = OrigValueStr) %>% 
  select(-c(SpeciesName, AccSpeciesID)) 
# On a plusieurs jeux de données d'origine, avec des valeurs d'indices différentes 
# par espèce. Je moyenne tout ça en première approximation.

ellenberg_fage <- TRY.fage.ellenberg3 %>% 
  group_by(AccSpeciesName) %>%
  summarize_at(c("light","moisture","nitrogen","pH","temperature"),mean,na.rm=T) 
ellenberg_annuals <- ellenberg_fage %>% 
  filter(AccSpeciesName %in% Annuals)

write.csv2(ellenberg_fage,"outputs/data/ellenberg_la_fage.csv",row.names = F)
write.csv2(ellenberg_annuals,"outputs/data/ellenberg_annuals.csv",row.names = F)

# A comparer avec les abondances dans le nat sup chez Maud
ab_maud_ann <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(LifeHistory=="annual") %>% 
  group_by(species) %>% 
  summarise(ab = mean(abundance)) %>% 
  arrange(-ab)

# Et avec les increasers et decreasers...

# Braun-Blanquetia ####
read.csv2("")
