source("scripts/1. Packages.R")
source("scripts/2. Import files.R")



#_______________________________________________________________________________
# 2) Add mean species attribute per treatment to abundance data per transect/quadrat ####
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

MEAN_tojoin <- MEAN %>% 
  dplyr::rename(species = Species, code_sp = Code_Sp, treatment = Trtmt)

# Abundance ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  # NB checker comment j'ai construit ce jeu de données.
  # Notamment, est-ce que j'ai bien les mêmes données d'abondance que Maud dans son papier de 2012 ? (ça pourrait expliquer mes résultats différents).
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))

ABUNDANCE_tojoin <- ABUNDANCE %>% 
  filter(!(code_sp %in% c("SOLCAI","SOLTER","SOLLIT"))) %>% 
  select(-c("LifeForm1","LifeForm2"))

ABUNDANCE_traits <- full_join(x= ABUNDANCE_tojoin,y = MEAN_tojoin, by = c("species","treatment","code_sp") )

write.csv2(ABUNDANCE_traits,"outputs/data/pooled_abundance_and_traits.csv",row.names=F)


#_______________________________________________________________________________
# CWM ####
ABUNDANCE_traits <- read.csv2("outputs/data/pooled_abundance_and_traits.csv")

# NB problem for Maud's data !
ABUNDANCE_CWM_treatment <- ABUNDANCE_traits %>% 
  group_by(treatment) %>% # NB choose the level at which to compute moments
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )

# Keep one row per species*treatment, with its trait value and the CWM
attribute_CWM_treatment <- ABUNDANCE_CWM_treatment %>% 
  distinct(species, treatment,.keep_all=T) %>% 
  select(-c("dataset" ,   "method" ,             "year"    ,            "paddock"  ,           "id_transect_quadrat",
            "abundance"    ,       "line_length"    , "soil"       ,         "depth")) # remove variables that


# Maud
soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))

# Maud <- ABUNDANCE_traits %>%
#   filter(dataset=="Maud") %>% 
#   group_by(paddock,depth) %>% 
#   arrange(depth) %>% 
#   filter(!(is.na(LifeHistory))) %>% 
#   full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
#   filter(!is.na(LifeHistory))
# # /!\ problem with raw data (ex: lack of Bupleurum)

ab_Maud_traits <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", sheet = "abondances par parcelle", 
                            startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(Species,abundance,-plot) %>%
  mutate(Species=str_replace(Species,"_"," ")) %>% 
  full_join(MEAN,by="Species") %>% 
  dplyr::rename(species=Species,code_sp=Code_Sp) %>% 
  filter(!is.na(abundance)) %>% 
  mutate(depth = str_sub(plot,start = 1L,end=1L)) %>% 
  mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
  select(-plot) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(Trtmt == "Nat") %>% 
  filter(abundance >0 )
write.csv2(ab_Maud_traits,"outputs/data/ab_Maud_traits.csv",row.names = F)

ab_Maud_traits %>% 
  select(code_sp,SLA,LDMC,L_Area) %>% 
  unique() %>% 
  filter(!(is.na(SLA))) %>%
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% 
  write.csv2("outputs/data/Pierce CSR/Traits_Maud.csv",row.names=F)


# Species level

ggplot(Maud, aes(x=LifeHistory, y=SLA))+
  geom_point() +
  facet_wrap(~depth)

Maud %>% filter(LifeHistory=="annual") %>%
  ggplot(aes(x=PC1score,y=SLA)) +
  geom_point() +
  geom_smooth(method="lm")

Maud %>% 
  group_by(paddock,depth,LifeHistory) %>% 
  summarize(abtot = sum(abundance)) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  spread(LifeHistory,abtot) %>% 
  mutate(relat_ab_annuals = annual/perennial) %>% 
  ggplot(aes(x=PC1score, y=relat_ab_annuals))+
  geom_point()

# Community level
# NB : check that relative values are used to compute the moments!!

vars <- MEAN %>% 
  select(Nb_Lf:Mat) %>% 
  colnames()

CWM_Maud <- ab_Maud_traits %>% 
  group_by(PC1score) %>% 
  mutate_at(vars, # vars: generated earlier in this script
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  # TAM::weighted_skewness for other moments
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  select_at(vars(PC1score, depth, paddock, contains( "CWM") )) %>% 
  unique()
write.csv2(CWM_Maud,"outputs/data/CWM_Maud.CSV",row.names=F)


gather_CWM_Maud<- CWM_Maud %>% 
  gather(key = trait,value=CWM,-c(depth,paddock,PC1score)) 
gather_CWM_Maud$trait <-   str_replace(gather_CWM_Maud$trait,"CWM_","")
gather_CWM_Maud2 <- gather_CWM_Maud %>% 
  arrange(factor(trait,levels=vars)) %>% 
  mutate(trait = factor(trait,levels=vars))

ggplot(gather_CWM_Maud2,aes(x=PC1score,y=CWM))+
  facet_wrap(~trait,scales="free") +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("CWM")


# Distance of species to CWM for various traits ####

CWM_Maud_indiv <- ab_Maud_traits %>% 
  group_by(PC1score) %>% 
  mutate_at(vars, # vars: generated earlier in this script
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )

CWM_Maud_indiv %>% 
  mutate(distance = abs(CWM_LDMC - LDMC)) %>% 
  ggplot(aes(x=LifeHistory,y=distance,color = LifeHistory))+
  geom_point() +
  facet_wrap(~depth)


CWM_Maud_indiv %>% 
  ggplot(aes(x=CWM_LDMC,y=LDMC,color = LifeHistory))+
  geom_point() +
  geom_abline(slope = 1, intercept = 0)
  
  facet_wrap(~depth)


# Fertile : distance of species to CWM for various traits ####

CWM_fer_indiv <- ABUNDANCE_traits %>% 
  filter(dataset=="Diachro" & year ==2004) %>% 
  group_by(id_transect_quadrat) %>% 
  mutate_at(vars, # vars: generated earlier in this script
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )

CWM_fer_indiv %>% 
  filter(!is.na(LifeHistory)) %>% 
  mutate(distance = abs(CWM_SLA - SLA)) %>% 
  ggplot(aes(x=LifeHistory,y=distance,color = LifeHistory))+
  geom_boxplot() +
  facet_wrap(~treatment)

CWM_fer_indiv %>% 
  filter(!is.na(LifeHistory)) %>% 
  mutate(distance = abs(CWM_Mat_Per - Mat_Per)) %>% 
  ggplot(aes(x=LifeHistory,y=Mat_Per,color = LifeHistory))+
  geom_boxplot() +
  facet_wrap(~treatment)

CWM_fer_indiv %>% 
  filter(!is.na(LifeHistory)) %>% 
  ggplot(aes(x=CWM_Mat_Per,y=Mat_Per,color = LifeHistory, shape = treatment))+
  geom_point()+
  facet_wrap(~LifeHistory) +
  geom_abline(slope = 1, intercept = 0)

# Using Maud's code ####


# data spring 2021 ####
data_traits_for_PCA <- LeafMorpho_leo %>% 
  filter(grepl("Nat",Treatment)) %>% 
  select(c(Code_Sp,SLA,LDMC,L_Area,Sple_FM,Sple_DM,Sple_Area)) %>% 
  group_by(Code_Sp) %>% 
  summarize_all(mean) %>% 
  column_to_rownames("Code_Sp")

ACP1<-PCA(data_traits_for_PCA,graph=FALSE)
coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column("code_sp")
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("code_sp")

ggplot(data=coord_ind,aes(x=Dim.1,y=Dim.2,label=code_sp)) +
  geom_point() +
  geom_text(data=coord_var, aes(x=Dim.1*5, Dim.2*5, label=code_sp), size = 5, vjust=1, color="black") +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_label()

