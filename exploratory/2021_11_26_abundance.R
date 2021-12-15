source("scripts/1. Packages.R")
source("scripts/2. Import files.R")



MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?"))

names_LH <- MEAN %>% 
  select(Code_Sp,Species,LifeHistory) %>% 
  unique()


# Maud ####
soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))


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
  merge(names_LH, by="Species")


annD <- ab_maud %>% 
  filter(LifeHistory=="annual") %>% 
  filter(depth == "D")%>% 
  pull(Species) %>% 
  unique()

annI <- ab_maud %>% 
  filter(LifeHistory== "annual") %>% 
  filter(depth == "I") %>% 
  pull(Species) %>% 
  unique()

annS <- ab_maud %>% 
  filter(LifeHistory=="annual") %>% 
  filter(depth == "S")
  pull(Species) %>% 
  unique()
  
  

setdiff(annS,annD)
intersect(annS,annD)
intersect(annS,annI)



richness_per_guild <- ab_maud %>% 
  count(depth,paddock,LifeHistory) %>% 
  merge(soil_Maud,by=c("depth","paddock"))


ggplot(richness_per_guild,aes(x=PC1score,y=n,fill=LifeHistory))+
  geom_bar(stat="identity",position = "dodge", width = 0.1)


richness_per_guild %>% 
  group_by(PC1score) %>% 
  # mutate(frac = n/sum(n)) %>% 
  spread(LifeHistory,n) %>% 
  replace(is.na(.), 0) 
  # summarise(var_annual = sd(annual)/mean(annual), var_perennial = sd(perennial)/mean(perennial))



richness_per_guild %>% 
  group_by(PC1score) %>% 
  mutate(frac = n/sum(n)) %>% 
  spread(LifeHistory,frac) %>% 
  replace(is.na(.), 0) %>% 
  summarise()
  # summarise(var_annual = sd(annual)/mean(annual), var_perennial = sd(perennial)/mean(perennial))

ggplot(ab_maud,aes(x=PC1score,y=abundance,color = LifeHistory))+
  geom_point()

ggplot(ab_maud %>% filter(depth == "S"),aes(x=abundance,color = LifeHistory))+
  geom_density()

# Rank-abundance diagram
Depth <- "S"
ab_maud2 <- ab_maud %>% 
  select(-paddock) %>% 
  group_by(Species,Code_Sp,depth,LifeHistory) %>% 
  summarize(abundance = sum(abundance)) %>% 
  filter(depth == Depth) %>% 
  arrange(-abundance) 
ab_maud2$rank <- seq.int(nrow(ab_maud2))

ggplot(ab_maud2,aes(x=rank,y=abundance,fill=LifeHistory))+
  geom_bar(stat="identity") +
  ggtitle(paste("Depth = ",Depth))+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text=element_text(size=18),
    legend.title = element_text(size=18))



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
ab_iris2$rank <- seq.int(nrow(ab_iris2))

ggplot(ab_iris2,aes(x=rank,y=abundance,fill=LifeHistory))+
  geom_bar(stat="identity") +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text=element_text(size=18),
    legend.title = element_text(size=18))


# Alexandre ####
# Ungrazed dans le fertile ! Pas la même chose qu'Ungrazed dans le natif, certes... 
# mais on peut s'attendre à ce qu'il y ait moins d'annuelles ?
ab_alex <- readODS::read_ods("data/abundance/Alexandre/leo_export.ods") %>% 
  dplyr::rename(gradient_level = 'gradient_level (Bas/Intermediaire/Haut)') %>% 
  merge(names_LH %>% rename(species = Species),by = "species")

ab_alex %>% 
  count(grazed,LifeHistory,gradient_level)

aaa <- ab_alex %>% 
  filter(grazed == "Grazed") %>% 
  filter((location %in% c("P6","P9","P8"))) %>% 
  group_by(location, gradient_level) %>% 
  count(LifeHistory)



# Generalism analysis ####
# site-by-species matrix
site_sp_nat <- ab_alex %>% 
  filter(grazed == "Grazed") %>% 
  select(location,gradient_level,Code_Sp,total_cover) %>% 
  filter(location %in% c("P6","P8","P9")) %>% 
  spread(key=Code_Sp,value = total_cover)  %>% 
  mutate(plot_depth = paste(location,gradient_level,sep="_")) %>% 
  column_to_rownames("plot_depth") %>% 
  select(-c("location","gradient_level")) %>% 
  ceiling()

rich <- c()
for (i in c(1:dim(site_sp_nat)[1])){
  rich <- c(rich,sum(site_sp_nat[i,] ,na.rm = T ))
}
cbind(rich,colnames(site_sp_nat))

ab_alex %>% 
  filter(grazed == "Grazed") %>% 
  select(LifeHistory,location,gradient_level,Code_Sp,total_cover) %>% 
  filter(location %in% c("P6","P8","P9"))


  