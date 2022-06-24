library(tidyverse)

# 1) List of annuals ####
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(species == "Geranium dissectum - limbe")) %>% 
  filter(!(species == "Geranium dissectum - pétiole")) %>% 
  filter(!(species == "Carex humilis?"))

names_LH <- MEAN %>% 
  select(code_sp,species,LifeHistory) %>% 
  filter(!code_sp=="STRIFSCAB") %>% 
  unique()

annuals <- names_LH %>% 
  filter(LifeHistory=="annual") %>% 
  select(species,code_sp) %>%  
  filter(!grepl("-",code_sp))

# write.csv2(annuals,"outputs/data/ellenberg_annuals_braund-blanquet.csv",row.names=F)
aaa <- read.csv2("outputs/data/ellenberg_annuals_braund-blanquet_completed.csv")

sp_manip <- c("BUPLBALD",
              "ALYSALYS","CAPSBURS","EROPVERN","HORNPETR",
              "ARENSERP","CERAGLOM","CERAPUMI","MINUHYBR",
              "MEDIMINI","TRIFCAMP",
              "GERADISS",
              "BROMHORD","VULPMYUR",
              "SHERARVE",
              "SAXITRIDA",
              "VEROARVE",
              "MYOSRAMOS",
              "FILAPYRAM")

# 2) Temporal evolution abundance ####
diachro <- read.csv2("data/abundance/diachro_releves_tidy2.csv")
# 33 transects, 16 dans le fertile, 17 dans le natif.

# fertilized treatment
diachro2 <- diachro %>% 
  filter(fertilized==T) %>% 
  group_by(line, year) %>% 
  mutate(rel_ab = abundance/sum(abundance)) %>% 
  group_by(species,year) %>% 
  mutate(delta_year = year - 1978, mean_ab = mean(abundance)) %>% 
  filter(LifeForm1=="The")

# Species not in temporal evolution (4 of my species : CERAGLOM, MYOSRAMOS, MINUHYBR, SAXITRIDA)
setdiff(annuals$species,
        diachro2$species %>% unique() )

# measure temporal evolution
diachro3 <- diachro2 %>% 
  ungroup() %>% 
  group_by(species) %>%
  filter(LifeForm1=="The") %>% 
  summarize(Rs = cor(delta_year,mean_ab,method = "spearman"),
            pval = cor.test(delta_year, mean_ab, method = "spearman")$p.value)  %>% 
  merge(names_LH,by="species") # add code_sp 

write.csv2(diachro3,"outputs/data/temporal_evolution_fer.csv",row.names = F)


# 3) Abundances ####

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  filter(LifeForm1=="The") %>% 
  rename(transect = id_transect_quadrat) %>% 
  mutate(code_sp = case_when(species == "Crepis vesicaria ssp. haenseleri" ~ "CREPVESI",
                             species == "Vicia sativa ssp. sativa" ~ "VICISATI",
                             TRUE ~ code_sp))

ab_nat <- read.csv2("outputs/data/abundance_natif.csv")%>% 
  filter(LifeHistory == "annual") %>% 
  mutate(transect=paste(paddock,depth,line,sep="_"))

Nat <- ab_nat %>% 
  group_by(code_sp) %>% 
  summarize(abundance_nat=mean(relat_ab))
Fer <- ab_fer %>% 
  group_by(code_sp) %>% 
  summarize(abundance_fer=mean(relat_ab))
NatFer <- full_join(Nat,Fer,by="code_sp") %>% 
  mutate(abundance_fer = abundance_fer %>%   replace(is.na(.), 0)) %>% 
  mutate(abundance_nat = abundance_nat %>%   replace(is.na(.), 0))
ggplot(NatFer,aes(x=abundance_nat,abundance_fer,label=code_sp)) +
  geom_point() +
  # ggrepel::geom_text_repel() +
  geom_abline(slope = 1,intercept=0) 


# Add info for species whose abundance = 0
ab_spread <- ab_fer %>% 
  select(code_sp,transect,relat_ab) %>%
  rbind( ab_nat %>% 
           select(code_sp,transect,relat_ab) ) %>% 
  spread(key=transect,value=relat_ab,fill=0)

missing_sp <- setdiff(annuals$code_sp,
        ab_spread %>% pull(code_sp) %>% unique() )

for (i in 1:length(missing_sp)){
  ab_spread <- ab_spread %>% 
    add_row(code_sp= missing_sp[i],.before=0)
}

ab_spread <- ab_spread %>% 
  replace(is.na(.),0) %>% 
  arrange(code_sp)

# gather the data frame again
ab_gathered <- ab_spread %>% 
  gather(key = transect,value = relat_ab,-code_sp) 

ab_summarized <- ab_gathered %>% 
  mutate(treatment = case_when(str_detect(transect,"F") ~ "Fer",
                                 TRUE ~ "Nat")) %>% 
  mutate(depth = case_when( str_detect(transect,"F") ~ "Fer",
                            str_detect(transect,"D") ~ "D",
                            str_detect(transect,"I") ~ "I",
                            str_detect(transect,"S") ~ "S" )) %>% 
  group_by(code_sp,depth,treatment) %>% 
  summarize(relat_ab = mean(relat_ab))

# plots
ggplot(ab_summarized,aes(x=code_sp,y=relat_ab, color=treatment, shape = depth)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

# 4) Ellenberg ####
ellenberg <- read.csv2("outputs/data/ellenberg_annuals.csv") %>% 
  rename(species = AccSpeciesName)

# 5) Baseflor ####
baseflor <- read.csv2("outputs/data/baseflor_annuals.csv") %>% 
  rename(species = species_1) %>% 
  rename(caract_ecol = "CARACTERISATION_ECOLOGIQUE_.HABITAT_OPTIMAL.") %>% 
  select(-species_2) %>% 
  unique()

baseflor_manip <- names_LH %>% 
  full_join(baseflor,by = "species") %>% 
  filter(code_sp %in% annuals$code_sp) %>% 
  select(species,code_sp,pollinisation,dissémination,caract_ecol) %>% 
  select(species,dissémination) %>% 
  unique() %>% 
  arrange(dissémination) %>% 
  filter(!dissémination == "")


# II) Group all ####
# Mean abundance per treatment ####
info_species <- ab_summarized %>% 
  select(-treatment) %>% 
  filter(depth %in% c("Fer","S")) %>% 
  spread(key = depth,value=relat_ab) %>%
  merge(names_LH,by="code_sp") %>% 
  full_join(ellenberg,by="species") %>% 
  select(-LifeHistory) %>% 
  full_join(diachro3 %>% select(-c(code_sp,LifeHistory)),by="species") %>% 
  relocate(species,code_sp,Fer,S,Rs,pval) %>% 
  full_join(baseflor_manip,by="species") %>% 
  filter(!species == "Geranium dissectum - pétiole")
  
#___________________________________
# ACP ####
library(FactoMineR)
rownames(info_species) <- NULL
ACP1 <- info_species %>% 
  distinct() %>% 
  column_to_rownames(var = "code_sp") %>% 
  select(nitrogen,moisture,light) %>% 
  PCA(graph = FALSE)

factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained

coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("code_sp")

# species with names
ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,label=code_sp)) +
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point(size=1) +
  ggrepel::geom_text_repel()

# avec axes
ggplot(coord_ind,aes(x=Dim.1,y=Dim.2)) +
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point(size=1) +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  ggrepel::geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 3, vjust=1, color="black") +
  theme(legend.position = "none") +
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) + 
  theme(text=element_text(size=10))

#________________________________________

#
info_species %>% 
  ggplot(aes(x=nitrogen,y=light,label=code_sp))+
  geom_point() +
  geom_jitter() +
  # geom_smooth() 
  ggrepel::geom_text_repel()

# nitrogen in S
info_species %>%
  filter(!is.na(nitrogen)) %>% 
  ggplot(aes(x=nitrogen,y=Rs,label=code_sp))+
  geom_point() +
  geom_smooth(method="lm") 
  ggrepel::geom_text_repel()

mod <- lm(S ~ nitrogen , 
          data = info_species %>% filter(!code_sp=="ARENSERP"))
plot(mod)
anova(mod)
summary(mod)


# Abundance per transect ####
info_species_transect <- ab_gathered %>% 
  filter(!str_detect(transect,"D")) %>% 
  filter(!str_detect(transect,"I")) %>% 
  mutate(treatment = case_when(str_detect(transect,"F") ~ "Fer",
                               TRUE ~ "Nat")) %>% 
  full_join(info_species %>% 
              select(-c(Fer,S)), 
            by = "code_sp")

# plot
info_species_transect %>% 
  filter(treatment == "Fer") %>% 
  ggplot(aes(x=dissémination,y= relat_ab,label=code_sp))+
  geom_boxplot() 


# modèle
mod <- lm(relat_ab ~ nitrogen + light, 
          data = info_species_transect %>% filter(treatment == "Nat"))
plot(mod)
anova(mod)
summary(mod)



info_species_transect %>% 
  filter(treatment == "Nat") %>% 
  ggplot(aes(x=nitrogen,y= relat_ab,label=code_sp))+
  geom_point() +
  geom_smooth(method="lm")



