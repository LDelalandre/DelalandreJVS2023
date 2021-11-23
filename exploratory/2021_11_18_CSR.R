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
               select(C,S,R), pch = 4, lwd=2)
AddToTernary(points, CSR_Maud %>% 
               filter(depth == "I") %>% 
               select(C,S,R), pch = 3, lty='dotted', lwd=2)
AddToTernary(points, CSR_Maud %>% 
               filter(depth == "S") %>% 
               select(C,S,R), pch = 1, lty='dotted', lwd=2)
legend("topright", 
       cex = 1,
       bty = "n",
       legend = c('Deep', 'Intermediary', 'Superficial'), 
       pch = c(4,3,1))



curve(x^(2) / 2,from = 0,to = 100,col = 'red',type = 'p',pch = 16,n = 20)
curve((1-x^(2))/2 + 5000,from = 0,to = 100,col = 'blue',type = 'p',pch = 15,add = TRUE,n = 20)


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
  group_by(plot,treatment) %>% 
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



#_____________________________________________________________________
# ACP
# Take 
# - the species*treatment level
# - the plot level (on various plots ?)
