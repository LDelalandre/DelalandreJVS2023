source("scripts/1. Packages.R")
source("scripts/2. Import files.R")

library("FD")


# Load data ####
bare <- read.csv(file = "data/abundance/Alexandre/data-export/bare.csv")
spab <- read.csv(file = "data/abundance/Alexandre/data-export/spab.csv") %>% 
  group_by(site) %>% 
  mutate(tot_ab = sum(abundance)) %>% 
  mutate(rel_ab = abundance/sum(abundance))

name_LH <- LeafMorpho %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  select(Species,Code_Sp,LifeHistory) %>% 
  unique()

name_LH_attribut <- name_LH %>% 
  mutate(genera_short = str_sub(Species, start = 1L, end = 3L) %>% toupper()) %>%
  mutate(first_space=  str_locate(Species," ")[,1]) %>% 
  mutate(sp_short =  str_sub(Species, start = .$first_space +1 , end =  .$first_space +2) %>% toupper()) %>% 
  mutate(attribut = paste0(genera_short,sp_short)) %>% 
  select(-c(genera_short,sp_short,first_space)) %>% 
  filter(!(Species %in% c("Myosostis ramosissima subsp. ramosissima","Vicia sativa subsp. sativa",
                          "Crepis sancta subsp. sancta","Carex humilis?",
                          "Geranium dissectum - limbe","Geranium dissectum - pétiole"))) # pour éviter les doublons


spab_lh <- merge(spab,name_LH_attribut,by="attribut") # life history added to Alexandre's abundance data

paddocks_nat <- c("P6","P9","P8")
paddocks_fer <- c("C1" ,"C3" ,"C2")
paddocks_tem <-  c("T1")


MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer

# Environment gps
gpsdat <- read.csv("data/abundance/Alexandre/data/gpsdat.csv") %>% 
  mutate(gps = paste(lat,lon,sep = " ")) %>% 
  select(site,gps) %>% 
  unique()

trans <- read.csv("data/abundance/Alexandre/data/trans.csv") 
envir <- read.csv("data/abundance/Alexandre/data/envir.csv") 

meteo <- read.csv2("data/Meteo_LaFage_1973-2006.csv")
mm_year <- meteo %>% group_by(AN) %>% 
  summarize(sumRR = sum(RR,na.rm=T))
hist(mm_year$sumRR)
mean(mm_year$sumRR)


#______________________________________________________________
# Plot ####
lhab_bare_splevel <- merge(spab_lh,bare,by="site") 

# Cover of annuals depending on environmental bare ground
lhab <- spab_lh %>% 
  group_by(site) %>% 
  mutate(ab_tot = sum(abundance), ab_relat = abundance/ab_tot) %>% 
  group_by(site,LifeHistory) %>%
  filter(location %in% c("P6","P9","P8")) %>% 
  summarize(AB = sum(abundance),AB_relat = sum(ab_relat)) 

lhab_bare <- merge(lhab,bare,by="site") %>% 
  filter(LifeHistory == "annual") 
  filter(grazed == "Grazed")
# filter(! (site %in% c("P9IGT2", "P8BGT1") )) # Influent sites

# transect info
ggplot(lhab_bare,aes(x=transect_bg,y = AB_relat,color = LifeHistory))+
  geom_point() +
  facet_wrap(~grazed) +
  geom_smooth(method = "lm") +
  labs(x="fraction of bare ground in the transect", y="relative abundance (fraction of the cover)") +
  ggtitle("Bare ground in the transect")

ggplot(lhab_bare,aes(x=transect_bg,y = AB,color = LifeHistory))+
  geom_point() +
  facet_wrap(~grazed) +
  geom_smooth(method = "lm") +
  labs(x="fraction of bare ground in the transect", y="absolute abundance (cm)") +
  ggtitle("Bare ground in the transect")

# environmental info
ggplot(lhab_bare,aes(x=envir_bg,y = AB_relat,color = LifeHistory))+
  geom_point() +
  facet_wrap(~grazed) +
  geom_smooth(method = "lm")  +
  labs(x="fraction of environmental bare ground", y="relative abundance (fraction of the cover)") +
  ggtitle("Environmental bare ground")

# dung
ggplot(lhab_bare,aes(x=dung,y = AB_relat,color = LifeHistory))+
  geom_point() +
  facet_wrap(~grazed) +
  geom_smooth(method = "lm") +
  labs(x="fraction of dung cover around the transect", y="relative abundance (fraction of the cover)") +
  ggtitle("Dung")

ggplot(lhab_bare,aes(x=litter,y = AB_relat,color = LifeHistory))+
  geom_point() +
  facet_wrap(~grazed) +
  geom_smooth(method = "lm")





#______________________________________________________________
# Identity of annuals ####



grazed_regime <- "Grazed"

spab_lh %>% 
  filter(location=="C1") %>% 
  pull(site) %>% 
  unique()

spab_fer <- spab_lh %>% 
  filter(location %in% paddocks_fer) %>% 
  filter(grazed == grazed_regime) %>% 
  filter(LifeHistory == "annual") %>%
  select(Species) %>% 
  unique()

spab_nat <- spab_lh %>% 
  filter(location %in% paddocks_nat) %>% 
  filter(grazed == grazed_regime) %>% 
  filter(LifeHistory == "annual") %>%
  select(Species) %>% 
  unique()

intersect(spab_fer,spab_nat)

setdiff(spab_fer,spab_nat)
setdiff(spab_nat,spab_fer)



spab_lh %>% 
  pull(location) %>% unique()







# Numbers of species
d1 <- spab_lh %>% 
  filter(grazed == "Grazed") %>%
  filter(location %in% c("P6","P9","P8")) %>% 
  group_by(location,gradient_level,LifeHistory) %>% 
  summarize(n=n()) %>% 
  spread(LifeHistory,n)

ggplot(d1,aes(x=gradient_level,y=perennial))+
  geom_point()+
  geom_point(aes(x=gradient_level,y=annual),color="red")

# abundance of species
d2 <- spab_lh %>% 
  filter(grazed == "Grazed") %>%
  filter(location %in% c("P6","P9","P8")) %>% 
  group_by(location,gradient_level,LifeHistory) %>%
  summarize(AB = sum(abundance)) %>% 
  spread(LifeHistory,AB)

ggplot(d2,aes(x=gradient_level,y=perennial))+
  geom_point()+
  geom_point(aes(x=gradient_level,y=annual),color="red")+
  facet_wrap(~location)


#_______________________________________________________________
# 2. Environment and com level ####

# CWM ####

# Abundance data
comm_nat <- spab_lh %>% 
  filter(location %in% paddocks_nat) %>% 
  select(site,Code_Sp,abundance) %>% 
  pivot_wider(names_from="Code_Sp", values_from = "abundance", values_fill = 0 ) %>% 
  column_to_rownames("site")
comm_nat <- comm_nat[,order(colnames(comm_nat))]
  
traits_nat <- MEAN %>% 
  filter(Trtmt == "Fer") %>% 
  select(-c("Species","Trtmt","Form","LifeHistory")) %>% 
  filter(Code_Sp %in% colnames(comm_nat)) %>% 
  arrange(Code_Sp) %>% 
  column_to_rownames(("Code_Sp"))
  

comm_nat <- comm_nat[,colnames(comm_nat)%in% rownames(traits_nat)]

FD <-dbFD(traits_nat, comm_nat, w.abun = TRUE,
          stand.x = TRUE, asym.bin = NULL,corr = "none",
          calc.FRic = FALSE, m = "max", stand.FRic = FALSE,
          scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
          km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100,
          calc.CWM = TRUE, calc.FDiv = FALSE, dist.bin = 2,print.pco = FALSE,
          messages = TRUE)

envir_per_transect <- envir %>% 
  select(-c(notes,X)) %>% 
  group_by(fertilized,site,location,gradient_level,grazed,transect) %>% 
  summarize_at(c("pebble","cobble","bedrock","dung","moss","litter","veg","bare","hmax"),mean)

CWM_nat <- FD$CWM %>% 
  rownames_to_column("site") %>% 
  left_join(envir_per_transect,by="site")

# CSR
CWM_nat_CSR <- CWM_nat %>% 
  select(site,L_Area,LDMC,SLA) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(CWM_nat_CSR,"outputs/data/Pierce CSR/Alex_CWM_natif.csv" ,row.names=F)

CWM_nat_CSR_completed <- read.csv2("outputs/data/Pierce CSR/Alex_CWM_natif_completed.csv")

# add info to CWM data
CWM_nat_CSR <- CWM_nat %>% 
  merge(CWM_nat_CSR_completed %>% select(site,C,S,R),.,by="site")

# abundance of annuals in the transects
lhab <- spab_lh %>% 
  group_by(site) %>% 
  mutate(ab_tot = sum(abundance), ab_relat = abundance/ab_tot) %>% 
  group_by(site,LifeHistory) %>%
  filter(location %in% c("P6","P9","P8")) %>% 
  summarize(AB = sum(abundance),AB_relat = sum(ab_relat))

CWM_nat_CSR_ab <- lhab %>% filter(LifeHistory == "annual") %>% 
  select(site,AB_relat) %>% 
  rename(AB_relat_annuals = AB_relat) %>% 
  merge(CWM_nat_CSR,by="site")



# PCA envir ####

# Links environment - CWM ####
# Trouver comment résumer les paramètres environnementaux (pour voir si on a un gradient chez Alex aussi comme chez Maud)
# --> ACP sur paramètres environnementaux ?

ggplot(CWM_nat,aes(x=bare+pebble,y=SLA,label = site ))+
  geom_point()


# Links CWM - cover of annuals ####
ggplot(CWM_nat_CSR_ab,aes(x=dung,y=AB_relat_annuals,label=site))+
  geom_label()

ggplot(CWM_nat_CSR_ab,aes(x=S,y=AB_relat_annuals,label=site))+
  geom_label()

# si j'enleve là où beaucoup de crottes (il va surtout falloir corriger pour ça et faire des modèles sympa)
low_dung <- CWM_nat_CSR_ab %>% 
  filter(dung < 1)

ggplot(low_dung,aes(x=C,y=AB_relat_annuals,label=site,color = dung))+
  geom_point()+
  scale_color_gradient(low="blue", high="red" )


# NB : essayer avec les scores CSR

# Links environment - cover of annuals ####

# Cover of annuals ####
ggplot(CWM_nat_CSR_ab,aes(x=AB_relat_annuals))+
  geom_density()

ggplot(CWM_nat_CSR_ab,aes(x =dung, y=AB_relat_annuals, label = site))+
  geom_text()

ggplot(CWM_nat_CSR_ab,aes(x=reorder(site,-AB_relat_annuals),y=AB_relat_annuals,fill=grazed))+
  geom_histogram(stat="identity") +
  theme(axis.text.x = element_text(angle = 45))  
