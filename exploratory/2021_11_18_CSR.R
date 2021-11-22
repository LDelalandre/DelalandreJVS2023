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

TernaryPlot(atip = 'Fertile', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points,  per_fer.list, col='blue', lty='dotted', lwd=3)
AddToTernary(points, an_fer.list, col='red', lty='dotted', lwd=3)
save

TernaryPlot(atip = 'Fertile', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(text, per_fer.list, names(per_fer.list), cex = 0.6, font = 2,col='blue')

TernaryPlot(atip = 'Fertile', # btip = 'S', ctip = 'R',
            alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(text, an_fer.list, names(an_fer.list), cex = 0.6, font = 2,col='red')


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


# Community Level ####

# Maud ####

# Plot level
CWM_Maud_goodunit <- read.csv2("outputs/data/CWM_Maud.CSV") %>% 
  select(PC1score,CWM_L_Area,CWM_LDMC,CWM_SLA) %>% 
  rename(SLA=CWM_SLA,LDMC=CWM_LDMC,L_Area=CWM_L_Area) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(CWM_Maud_goodunit,"outputs/data/Pierce CSR/Traits_CWM_Maud.csv" ,row.names=F)

CSR_Maud <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_Maud_completed.csv") %>% 
  select(C,S,R)

TernaryPlot(main = "Maud, community level",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_Maud, col='black', lty='dotted', lwd=3)

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
