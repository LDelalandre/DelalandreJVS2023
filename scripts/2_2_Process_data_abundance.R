source("scripts/Packages.R")
names_LH2 <- read.csv2("data/species_names_lifehistory.csv") %>% 
  mutate(species= Species, code_sp = Code_Sp)

# Fertilized treatment ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))
# I use abondance data from the 2005 survey
ab_diachro_2005 <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year == 2005) %>%  # /!\ No measure available after 2005! Should I keep 2006, or several years?
  filter(treatment == "Fer") %>%
  # add relative abundance
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance))

write.csv2(ab_diachro_2005,"outputs/data/abundance_fertile.csv",row.names=F)

ab_diachro_2004 <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year == 2004) %>%  # /!\ No measure available after 2005! Should I keep 2006, or several years?
  filter(treatment == "Fer") %>%
  # add relative abundance
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance))




# Unfertilized treatment #### 
# (Maud's relevés)
soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))
# par ligne
ab_maud <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", 
                     sheet = "abondance par ligne", 
                     startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(species,abundance,-Ligne) %>%
  mutate(species=str_replace(species,"_"," ")) %>% 
  mutate(depth = str_sub(Ligne,start = -2L,end=-2L)) %>% 
  mutate(line = str_sub(Ligne,start = -1L,end=-1L)) %>% 
  mutate(paddock = str_sub(Ligne,start = 1L,end=-3L)) %>%
  # select(-Ligne) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(abundance >0 ) %>% 
  merge(names_LH2, by="species") %>% 
  # add relative abundance
  group_by(line,depth,paddock) %>% 
  filter(!code_sp =="STRIFSCAB")%>% # /!\ à corriger à la source !!
  mutate(relat_ab = abundance/sum(abundance)) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

write.csv2(ab_maud,"outputs/data/abundance_natif.csv",row.names=F)


# Comparison abundances between the two measurements ####
ab_diachro_2004_nat <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year == 2004) %>%  # /!\ No measure available after 2005! Should I keep 2004 instead of 2005 ?
  # voir ce que ça change !
  filter(treatment == "Nat") %>%
  # add relative abundance
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance))

ab_maud_mean <- ab_maud %>% 
  group_by(code_sp,LifeHistory) %>% 
  summarize(mean_relat_ab_maud = mean(relat_ab))

ab_maud_mean_noS <- ab_maud %>% 
  group_by(code_sp,LifeHistory) %>%
  filter(!depth=="S") %>% 
  summarize(mean_relat_ab_maud = mean(relat_ab))

ab_diachro_nat_mean <- ab_diachro_2004_nat %>% 
  group_by(paddock,id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance)) %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  ungroup() %>% 
  group_by(code_sp,LifeHistory) %>% 
  summarize(mean_relat_ab_diachro = mean(relat_ab))

compare_ab <- full_join(ab_maud_mean,ab_diachro_nat_mean,by=c("code_sp","LifeHistory"))  
compare_ab_noS<- full_join(ab_maud_mean_noS,ab_diachro_nat_mean,by=c("code_sp","LifeHistory"))  

ggplot(compare_ab_noS,aes(x=log(mean_relat_ab_maud),y=log(mean_relat_ab_diachro),
                      label=code_sp))+ # ,color=LifeHistory
  geom_point()  +
  geom_smooth(method="lm")
  # ggrepel::geom_text_repel()

mod <- lm(log(mean_relat_ab_maud)~log(mean_relat_ab_diachro),data = compare_ab) 
shapiro.test(residuals(mod)) # normality of residuals
lmtest::bptest(mod) # no homoscedasticity
lmtest::dwtest(mod) # autocorrelation
summary(mod)
anova(mod)

mod$coefficients
sum <- summary(mod)
sum$r.squared
sum$adj.r.squared
sum$coefficients

anov <- anova(mod)
as.data.frame(anov)
rownames(anov) <- c("Species mean abundance (2009)","Residuals")

# export table
table_anova <- anov %>%
  kableExtra::kable( escape = F) %>%
  kableExtra::kable_styling("hover", full_width = F)

cat(table_anova, file = "draft/abundance_2009_2004.doc")



# Comparison abundance annuals in both treatments ####
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
ggplot(NatFer,aes(x= abundance_nat ,y=abundance_fer,label=code_sp)) +
  geom_point() +
  # ggrepel::geom_text_repel() +
  geom_abline(slope = 1,intercept=0) +
  theme_classic() +
  labs(x= "Annuals relative abundance in sandy GU",
       y= "Annuals relative abundance in G+F") 

NatFer2 <- NatFer %>% 
  filter(abundance_fer>0) %>% 
  filter(abundance_nat>0)
cor.test(NatFer2$abundance_nat,NatFer2$abundance_fer)
# correlation negative, but not significant.
