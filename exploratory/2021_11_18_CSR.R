source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
library('Ternary') # https://mran.microsoft.com/web/packages/Ternary/vignettes/Ternary.html


name_LH <- LeafMorpho %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  select(Code_Sp,LifeHistory) %>% 
  unique()
  
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?"))

# On species (per treatment) ####
MEAN_goodunit <- MEAN %>% 
  select(Species,Code_Sp,Trtmt,L_Area,LDMC,SLA) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(MEAN_goodunit,"outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt.csv" ,row.names=F)


MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  merge(name_LH,by="Code_Sp") %>% 
  relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())


make_list <- function(xy.df){
  # xy.list <- split(xy.df, seq(nrow(xy.df)))
  # xy.list <- setNames(split(xy.df, seq(nrow(xy.df))), rownames(xy.df))
  xy.list <- as.list(as.data.frame(t(xy.df)))
  xy.list
}
# NB : not necessary to make a list for plotting with the Ternary package


to_boxplot <- MEAN_CSR %>% 
  select(-c(LDMC,SLA,L_Area,Red,Green,Blue,Species)) %>% 
  gather(Dimension, Score,-c(Code_Sp,LifeHistory,Trtmt)) %>% 
  filter(Trtmt %in% c("Fer","Nat"))
to_boxplot$Score = as.numeric(to_boxplot$Score)
to_boxplot$Dimension <- as.factor(to_boxplot$Dimension)

ggplot(to_boxplot,aes(x=as.factor(Dimension),y=Score,color=LifeHistory))+
  geom_boxplot() +
  facet_wrap(~Trtmt) 
  # ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                        # map_signif_level = TRUE,vjust = 0.5,col="black")

# mod <- lm(Score ~ LifeHistory*Dimension*Trtmt,data=to_boxplot)
# anova(mod)

# Ternary diagrams
an_fer.df <- MEAN_CSR %>% 
  filter(Trtmt == "Fer" & LifeHistory == "annual") %>% 
  column_to_rownames("Code_Sp") %>%
  select(C,S,R)
an_fer.list <- make_list(an_fer.df)

an_nat.df <- MEAN_CSR %>% 
  filter(Trtmt == "Nat" & LifeHistory == "annual") %>% 
  column_to_rownames("Code_Sp") %>%
  select(C,S,R)
an_nat.list <- make_list(an_nat.df)

per_fer.df <- MEAN_CSR %>% 
  filter(Trtmt == "Fer" & LifeHistory == "perennial") %>% 
  column_to_rownames("Code_Sp") %>%
  select(C,S,R)
per_fer.list <- make_list(per_fer.df)

per_nat.df <- MEAN_CSR %>% 
  filter(Trtmt == "Nat" & LifeHistory == "perennial") %>% 
  column_to_rownames("Code_Sp") %>%
  select(C,S,R)
per_nat.list <- make_list(per_nat.df)


# Species level ####
TernaryPlot(main = "La Fage",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, MEAN_CSR %>% 
               filter(LifeHistory == "perennial") %>% 
               select(C,S,R), col='blue', lty='dotted', lwd=3)

# TernaryPlot(main = "La Fage",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, MEAN_CSR %>% 
               filter(LifeHistory == "annual") %>% 
               select(C,S,R), col='red', lty='dotted', lwd=3)

# Species level per treatment ####
par(mfrow=c(1,1))

# Fertile
TernaryPlot(atip = 'Fertile', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,  per_fer.list, col='blue', lty='dotted', lwd=3)
AddToTernary(points, an_fer.list, col='red', lty='dotted', lwd=3)


TernaryPlot(atip = 'Fertile', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(text, per_fer.list, names(per_fer.list), cex = 0.6, font = 2,col='blue')

TernaryPlot(atip = 'Fertile', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(text, an_fer.list, names(an_fer.list), cex = 0.6, font = 2,col='red')

# Natif
TernaryPlot(atip = 'Natif', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,  per_nat.list, col='blue', lty='dotted', lwd=3)
AddToTernary(points, an_nat.list, col='red', lty='dotted', lwd=3)

TernaryPlot(atip = 'Natif', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(text,  per_nat.list, labels = names(per_nat.list), col='blue', lty='dotted', lwd=3, cex = 0.6)

TernaryPlot(atip = 'Natif', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(text, an_nat.list, labels = names(an_nat.list), col='red', lty='dotted', lwd=3, cex = 0.6)

# Temoin
TernaryPlot(main = "Temoin",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, MEAN_CSR %>% 
               filter(Trtmt == "Tem") %>% 
               filter(LifeHistory == "perennial") %>% 
               select(C,S,R), col='blue', lty='dotted', lwd=3)


# Community Level ####

# Maud ####

# Plot level
CWM_Maud_goodunit <- read.csv2("outputs/data/CWM_Maud.CSV") %>% 
  select(PC1score,depth,CWM_L_Area,CWM_LDMC,CWM_SLA) %>% 
  rename(SLA=CWM_SLA,LDMC=CWM_LDMC,L_Area=CWM_L_Area) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(CWM_Maud_goodunit,"outputs/data/Pierce CSR/Traits_CWM_Maud.csv" ,row.names=F)

CSR_Maud <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_Maud_completed.csv") 

TernaryPlot(main = "Maud, community level",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_Maud %>% 
               filter(depth == "D") %>% 
               select(C,S,R), pch = 1, lwd=2)
AddToTernary(points, CSR_Maud %>% 
               filter(depth == "I") %>% 
               select(C,S,R), pch = 4, lty='dotted', lwd=2)
AddToTernary(points, CSR_Maud %>% 
               filter(depth == "S") %>% 
               select(C,S,R), pch = 3, lty='dotted', lwd=2)
legend("topright", 
       cex = 1,
       bty = "n",
       legend = c('Deep', 'Intermediary', 'Superficial'), 
       pch = c(1,4,3))

#____
# relative abundance of annuals (maud data: pin-point) ####
ab_maud_traits <- read.csv2("outputs/data/ab_Maud_traits.csv")
rel_ab_annuals <- ab_maud_traits %>% 
  group_by(depth,paddock,PC1score) %>% 
  select(PC1score,depth,paddock,species,abundance,LifeHistory) %>% 
  mutate(rel_ab = abundance/sum(abundance)) %>% 
  filter(LifeHistory == "annual") %>% 
  summarize(sum_ab_ann = sum(rel_ab)) %>% 
  ungroup() %>% 
  select(-c("depth"))

CSR_PC1_cover <- left_join(CSR_Maud,rel_ab_annuals,by = "PC1score") %>% 
  replace(is.na(.), 0)
# Rq: cover of annuals is 14% in parc 1, which is the most windy of La Fage

# Environmental info
ggplot(CSR_PC1_cover,aes(x=-PC1score,y=sum_ab_ann,label=depth))+
  geom_point() +
  geom_label()

# CSR info
CSR_PC1_cover_removedinfluent <- CSR_PC1_cover %>% filter(!(paddock=="P1"))

ggplot(CSR_PC1_cover, aes(x=S,y=sum_ab_ann,label = depth))+
  geom_label() 

  
# CSR and environment
ggplot(CSR_PC1_cover,aes(x=-PC1score,y=S))+
  geom_point()
ggplot(CSR_PC1_cover,aes(x=-PC1score,y=C))+
  geom_point()
ggplot(CSR_PC1_cover,aes(x=-PC1score,y=R))+
  geom_point()

  
# en boxplot
ggplot(CSR_PC1_cover,aes(x=depth,y=sum_ab_ann))+
  geom_boxplot()

ggplot(CSR_PC1_cover,aes(x=depth,y=S))+
  geom_boxplot()
ggplot(CSR_PC1_cover,aes(x=depth,y=C))+
  geom_boxplot()
ggplot(CSR_PC1_cover,aes(x=depth,y=R))+
  geom_boxplot()
  


# Ajouter info Adeline fertile
Adeline_ab_tr <- read.csv2("outputs/data/pooled_abundance_and_traits.csv") %>% 
  filter(dataset == "Adeline") %>% 
  group_by(paddock,id_transect_quadrat) # NB choose the level at which to compute moments. group, or  plot...
CSR_Adeline_CWM <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_Adeline_completed.csv")  %>% 
  mutate(treatment = case_when(
    str_detect(paddock, "C") ~ "Fer",
    str_detect(paddock, "N") ~ "Nat",
    str_detect(paddock, "T") ~ "Tem"))

adeline_ab_CSR <- Adeline_ab_tr %>% 
  group_by(id_transect_quadrat) %>% 
  mutate(ab_relat = abundance/sum(abundance)) %>% 
  merge(CSR_Adeline_CWM %>% select(-c("SLA","LDMC","L_Area","paddock","treatment")),by = "id_transect_quadrat")

adeline_ab_CSR_toplot <- adeline_ab_CSR %>% 
  select(soil ,id_transect_quadrat,species,code_sp,treatment,LifeHistory,ab_relat,C,S,R) %>% 
  group_by(soil,C,S,R,id_transect_quadrat,treatment,LifeHistory) %>% 
  summarize(sum_ab_relat = sum(ab_relat)) %>% 
  filter(LifeHistory == "annual" & treatment== 'Fer')

ggplot( adeline_ab_CSR_toplot ,
        aes(x= S  ,y=sum_ab_relat)) +
  geom_point()


  
# Merge with Maud data
maud_tomerge <- CSR_PC1_cover %>% 
  select(PC1score, depth,         C ,       S  ,        R, sum_ab_ann)
adeline_tomerge <- adeline_ab_CSR_toplot %>% 
  mutate(depth = "F", PC1score = "NA") %>% 
  rename(sum_ab_ann = sum_ab_relat) %>% 
  ungroup() %>% 
  select(PC1score, depth,         C ,       S  ,        R, sum_ab_ann)

rbind(maud_tomerge,adeline_tomerge) %>% 
  ggplot(aes(x=reorder(depth,-sum_ab_ann),y=sum_ab_ann))+
  geom_boxplot()

rbind(maud_tomerge,adeline_tomerge) %>% 
  arrange(factor(depth, levels = c("F","D","I","S"))) %>% 
  mutate(depth = factor(depth,levels =  c("F","D","I","S"))) %>% 
  ggplot(aes(x=depth,y=C))+
  geom_boxplot()


# densité abundance annuals
ggplot(CSR_PC1_cover,aes(x=sum_ab_ann))+
  geom_density()

# Caractériser les stratégies autrement que juste CSR ? Plus de traits fonctionnels (par exemple phéno ?)


# Traits des annuelles le long du gradient
ggplot(ab_maud_traits %>% filter(LifeHistory == "annual"),aes(x=PC1score,y=LDMC,color = species))+
  geom_point()

#____

# Compare with species level 
Depth <- "S"

ab_Maud_traits_depth <- read.csv2("outputs/data/ab_Maud_traits.csv") %>% 
  filter(depth==Depth)

CSR_Maud_sp <- read.csv2("outputs/data/Pierce CSR/Traits_Maud_completed.csv") %>% 
  merge(name_LH %>% rename(code_sp=Code_Sp),by="code_sp") %>% 
  filter(code_sp %in% ab_Maud_traits_depth$code_sp)

CSR_Maud_sp %>% 
  filter(LifeHistory=="annual")

# Points
TernaryPlot(main = paste("Maud, species level /",Depth), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,CSR_Maud_sp %>% 
              filter(LifeHistory=="annual") %>% 
              select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
AddToTernary(points,CSR_Maud_sp %>% 
              filter(LifeHistory=="perennial")%>% 
              select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)

# Species names
TernaryPlot(main = "Maud, species level",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
CSR_Maud_sp %>% 
  filter(LifeHistory=="perennial")%>% 
  column_to_rownames("code_sp") %>% 
  select(C,S,R) %>% 
  AddToTernary(text, . , labels = rownames(.),cex=0.6, col='blue',lty='dotted',lwd=3)

TernaryPlot(main = "Maud, species level",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
CSR_Maud_sp %>% 
  filter(LifeHistory=="annual")%>% 
  column_to_rownames("code_sp") %>% 
  select(C,S,R) %>% 
  AddToTernary(text, . , labels = rownames(.),cex=0.6, col='red',lty='dotted',lwd=3)

#_______________________________________________________________________________________
# Iris ####
ab_iris <- read.xlsx("data/abundance/iris_Relevés bota.xlsx",sheet="Tout compilé") %>% # NB : je peux aussi importer ses données du natif (et du)
  dplyr::rename(abundance = 'Abondance.(tot)',species = Espece) %>% 
  mutate(plot = paste(Parc,Cage,sep="")) %>% 
  mutate(treatment = case_when(Traitement == "P+F+" ~ "Fer",
                               Traitement == "P+F-" ~ "Nat",
                               Traitement == "P-F-" ~ "Tem")) 

CWM_iris <- MEAN %>% 
  # merge(name_LH %>% dplyr::rename(code_sp=Code_Sp),by="code_sp") %>% 
  # select(-c(LifeHistory,Code_Sp)) %>% 
  merge(ab_iris %>% dplyr::rename(Species=species,Trtmt = treatment),by=c("Species","Trtmt")) %>% 
  group_by(plot) %>% # NB choose the level at which to compute moments. group, or  plot...
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique() 

Iris_sp <- CWM_iris %>% 
  group_by(Code_Sp,Trtmt,LifeHistory) %>% 
  select(L_Area,LDMC,SLA) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(Iris_sp,"outputs/data/Pierce CSR/Traits_Iris.csv" ,row.names=F)



# Points
CSR_Iris_sp <- read.csv2("outputs/data/Pierce CSR/Traits_Iris_completed.csv")
# Temoin
TernaryPlot(main = paste("Iris, species level / Temoin"), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,CSR_Iris_sp %>%
               filter(treatment == "Tem") %>% 
               filter(LifeHistory=="annual") %>% 
               select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
AddToTernary(points,CSR_Iris_sp %>% 
               filter(LifeHistory=="perennial")%>%
               filter(treatment == "Tem") %>% 
               select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)

# Natif
TernaryPlot(main = paste("Iris, species level / Natif"), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_Iris_sp %>%
               filter(LifeHistory=="annual") %>% 
               filter(treatment == "Nat") %>% 
               select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
AddToTernary(points,CSR_Iris_sp %>% 
               filter(LifeHistory=="perennial")%>%
               filter(treatment == "Nat") %>% 
               select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)

# Fertile
TernaryPlot(main = paste("Iris, species level / Fertile"), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,CSR_Iris_sp %>%
               filter(treatment == "Fer") %>% 
               filter(LifeHistory=="annual") %>% 
               select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
AddToTernary(points,CSR_Iris_sp %>% 
               filter(LifeHistory=="perennial")%>%
               filter(treatment == "Fer") %>% 
               select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)


# Community level #
CWM_iris2 <- CWM_iris %>% 
  group_by(plot,Trtmt) %>% 
  select(starts_with("CWM")) %>% 
  unique() %>% 
  select(CWM_L_Area,CWM_LDMC,CWM_SLA) %>% 
  rename(SLA=CWM_SLA,LDMC=CWM_LDMC,L_Area=CWM_L_Area) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(CWM_iris2,"outputs/data/Pierce CSR/Traits_CWM_Iris.csv" ,row.names=F)

CSR_Iris <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_Iris_completed.csv") 

TernaryPlot(main = "Iris, community level",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_Iris %>% 
               filter(treatment == "Fer") %>% 
               select(C,S,R), col='green', lty='dotted', lwd=3)
AddToTernary(points, CSR_Iris %>% 
               filter(treatment == "Nat") %>% 
               select(C,S,R), col='orange', lty='dotted', lwd=3)
AddToTernary(points, CSR_Iris %>% 
               filter(treatment == "Tem") %>% 
               select(C,S,R), col='black', lty='dotted', lwd=3)


# In the natif 
TernaryPlot(main = "Iris, community level, natif",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
legend("topright", 
       cex = 1,
       bty = "n",
       legend = c('Deep', 'Intermediary', 'Superficial'), 
       pch = c(1,4,3))

AddToTernary(points, CSR_Iris %>% 
               filter(grepl("D",plot)) %>% 
               filter(treatment == "Nat") %>% 
               select(C,S,R), col='black', lty='dotted', lwd=3, pch = 1)
AddToTernary(points, CSR_Iris %>% 
               filter(grepl("I",plot)) %>% 
               filter(treatment == "Nat") %>% 
               select(C,S,R), col='black', lty='dotted', lwd=3, pch = 4)
AddToTernary(points, CSR_Iris %>% 
               filter(grepl("S",plot)) %>% 
               filter(treatment == "Nat") %>% 
               select(C,S,R), col='black', lty='dotted', lwd=3, pch = 3)

#______________________________________________________________________________
# Adeline ####
Adeline_ab_tr <- read.csv2("outputs/data/pooled_abundance_and_traits.csv") %>% 
  filter(dataset == "Adeline") %>% 
  group_by(paddock,id_transect_quadrat) # NB choose the level at which to compute moments. group, or  plot...

Adeline_sp <- Adeline_ab_tr %>% 
  group_by(code_sp,treatment,LifeHistory) %>% 
  select(L_Area,LDMC,SLA) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA)) %>% 
  unique()
write.csv2(Adeline_sp,"outputs/data/Pierce CSR/Traits_Adeline.csv" ,row.names=F)

CWM_adeline0 <- Adeline_ab_tr %>% 
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()

CWM_adeline <- CWM_adeline0 %>% 
  select(CWM_L_Area,CWM_LDMC,CWM_SLA) %>% 
  rename(SLA=CWM_SLA,LDMC=CWM_LDMC,L_Area=CWM_L_Area) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA)) %>% 
  unique()
write.csv2(CWM_adeline,"outputs/data/Pierce CSR/Traits_CWM_Adeline.csv" ,row.names=F)


# Plots

# Points
CSR_Adeline_sp <- read.csv2("outputs/data/Pierce CSR/Traits_Adeline_completed.csv")
# Temoin
TernaryPlot(main = paste("Adeline, species level / Temoin"), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
# AddToTernary(points,CSR_Adeline_sp %>%
#                filter(treatment == "Tem") %>%
#                filter(LifeHistory=="annual") %>%
#                select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
# # NO ANNUAL IN THE TEMOIN
AddToTernary(points,CSR_Adeline_sp %>% 
               filter(LifeHistory=="perennial")%>%
               filter(treatment == "Tem") %>% 
               select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)

# Natif
TernaryPlot(main = paste("Adeline, species level / Natif"), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_Adeline_sp %>%
               filter(LifeHistory=="annual") %>%
               filter(treatment == "Nat") %>%
               select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
AddToTernary(points,CSR_Adeline_sp %>% 
               filter(LifeHistory=="perennial")%>%
               filter(treatment == "Nat") %>% 
               select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)

# Fertile
TernaryPlot(main = paste("Adeline, species level / Fertile"), alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,CSR_Adeline_sp %>%
               filter(treatment == "Fer") %>% 
               filter(LifeHistory=="annual") %>% 
               select(C,S,R), cex=0.8, col='red',lty='dotted',lwd=3)
AddToTernary(points,CSR_Adeline_sp %>% 
               filter(LifeHistory=="perennial")%>%
               filter(treatment == "Fer") %>% 
               select(C,S,R), cex=0.8, col='blue',lty='dotted',lwd=3)


# Community level #

CSR_Adeline_CWM <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_Adeline_completed.csv")  %>% 
  mutate(treatment = case_when(
    str_detect(paddock, "C") ~ "Fer",
    str_detect(paddock, "N") ~ "Nat",
    str_detect(paddock, "T") ~ "Tem"))

TernaryPlot(main = "Adeline, community level",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_Adeline_CWM %>% 
               filter(treatment == "Fer") %>% 
               select(C,S,R), col='green', lty='dotted', lwd=3)
AddToTernary(points, CSR_Adeline_CWM %>% 
               filter(treatment == "Nat") %>% 
               select(C,S,R), col='orange', lty='dotted', lwd=3)
AddToTernary(points, CSR_Adeline_CWM %>% 
               filter(treatment == "Tem") %>% 
               select(C,S,R), col='black', lty='dotted', lwd=3)


# CWM CSR and abundance Adeline ####

adeline_ab_CSR <- Adeline_ab_tr %>% 
  group_by(id_transect_quadrat) %>% 
  mutate(ab_relat = abundance/sum(abundance)) %>% 
  merge(CSR_Adeline_CWM %>% select(-c("SLA","LDMC","L_Area","paddock","treatment")),by = "id_transect_quadrat")

adeline_ab_CSR_toplot <- adeline_ab_CSR %>% 
  select(soil ,id_transect_quadrat,species,code_sp,treatment,LifeHistory,ab_relat,C,S,R) %>% 
  group_by(soil,C,S,R,id_transect_quadrat,treatment,LifeHistory) %>% 
  summarize(sum_ab_relat = sum(ab_relat))


adeline_ab_CSR_toplot %>% 
  filter(LifeHistory == "annual" & treatment== 'Nat') %>% 
  arrange(sum_ab_relat)

ggplot( adeline_ab_CSR_toplot %>% 
         filter(LifeHistory == "annual" & treatment== 'Nat') ,
       aes(x= S  ,y=sum_ab_relat)) +
  geom_point()

#_____________________________________________________________________
# ACP
# Take 
# - the species*treatment level
# - the plot level (on various plots ?)
