source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
library('Ternary') # https://mran.microsoft.com/web/packages/Ternary/vignettes/Ternary.html


name_LH <- LeafMorpho %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  select(Code_Sp,LifeHistory)
  
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

# On species
MEAN_goodunit <- MEAN %>% 
  select(Species,Code_Sp,Trtmt,L_Area,LDMC,SLA) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA)) %>% 
  mutate(SLA=SLA*100/1000000) # change from cm²/kg to mm²/mg
write.csv2(MEAN_goodunit,"outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt.csv" ,row.names=F)


MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  full_join(name_LH,by="Code_Sp") %>% 
  relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())

str(MEAN_CSR)


coordss <- list(
  A = c(1, 0, 2),
  B = c(1, 1, 1),
  C = c(1.5, 1.5, 0),
  D = c(0.5, 1.5, 1)
)

coords <- data.frame(
  A = c(1, 0, 0),
  B = c(0, 1, 0)
) %>% 
  as.matrix %>% 
  t() %>% 
  as.data.frame()

an_fer <- MEAN_CSR %>% 
  filter(Trtmt == "Fer" & LifeHistory == "annual") %>% 
  select(C,S,R)

an_nat <- MEAN_CSR %>% 
  filter(Trtmt == "Nat" & LifeHistory == "annual") %>% 
  select(C,S,R)

per_fer <- MEAN_CSR %>% 
  filter(Trtmt == "Fer" & LifeHistory == "perennial") %>% 
  select(C,S,R)

per_nat <- MEAN_CSR %>% 
  filter(Trtmt == "Nat" & LifeHistory == "perennial") %>% 
  select(C,S,R)

par(mfrow=c(1,2))

TernaryPlot(atip = 'Annual', # btip = 'S', ctip = 'R',
            alab = 'C', blab = 'S', clab = 'R')
AddToTernary(points, an_fer, col='green', lty='dotted', lwd=3)
AddToTernary(points, an_nat, col='darkgreen', lty='dotted', lwd=3)

TernaryPlot( atip = 'Perennial',#, btip = 'S', ctip = 'R',
  alab = 'C', blab = 'S', clab = 'R')
AddToTernary(points, per_fer, col='lightblue', lty='dotted', lwd=3)
AddToTernary(points, per_nat, col='darkblue', lty='dotted', lwd=3)

# Community Level
ABUNDANCE_traits <- read.csv2("outputs/data/pooled_abundance_and_traits.csv")


ABUNDANCE_CWM_treatment <- ABUNDANCE_traits %>% 
  group_by(treatment) %>% # NB choose the level at which to compute moments
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )


