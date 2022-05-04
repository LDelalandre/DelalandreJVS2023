source("scripts/Packages.R")

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  # merge(name_LH,by="Code_Sp") %>% 
  # relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())

# MEAN_CSR2 <- MEAN_CSR %>% 
#   select(LifeHistory,C,S,R)
# ATTENTION on peut aussi utiliser ce MEAN là (recalculer les scores CSR dans ce cas) :
# MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv")%>%
#   filter(!is.na(SLA))

#__________________________________________________
# Option 1 : prendre les traits des espèces mesurés dans l'un ou l'autre traitement ####
# En comparant nat_sab et fer

# i) CSR ####
boxplot_CSR <- MEAN_CSR %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  # filter(score == "R") %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU")) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  theme_classic()+
  ylab("Score (%)")+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  facet_wrap(~zone) 

boxplot_CSR

ggsave("draft/boxplot_CSR.jpg",boxplot_CSR)

data.anovaCSR <- MEAN_CSR %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(treatment%in% c("Nat","Fer"))

options(contrasts=c("contr.treatment","contr.treatment"))
lmCSR <- lm(R ~ treatment * LifeHistory, data = data.anovaCSR)

par(mfrow=c(2,2)) ; plot(lmCSR) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(lmCSR))) # normality_graph
shapiro.test(residuals(lmCSR)) # normality of residuals
lmtest::bptest(lmCSR) # homoscedasticity
lmtest::dwtest(lmCSR)  # non-independence of residuals!
anova(lmCSR)
summary(lmCSR)

# Interprétation facile 


# ii) Autres traits ####

trait <- "SLA"
trait <- "R"
trait <- "Disp"

MEAN_CSR %>% 
  # select(code_sp,treatment,LifeHistory,SeedMass) %>%
  filter(treatment %in% c("Nat","Fer")) %>% 
  ggplot(aes_string(x="LifeHistory", y=trait, color = "LifeHistory"))+
  # geom_point()+
  geom_boxplot() +
  facet_wrap(~treatment) +
  ggsignif::geom_signif()


# Option 2 : MEAN calulé dans Fer et Nat Sup ####

#____________________________________________________________
# Option 3 : sélectionner les espèces qui apparaissent dans les relevés de maud superficiel et diachro ####
# C'est peut-être le plus cohérent, mais on n'a pas assez de points pour
# que les tests soient suffisamment puissants...

# dans le fertile
species_fer <- ab_fer %>% 
  pull(code_sp) %>% 
  unique()

MEAN_CSR %>% 
  filter(code_sp %in% species_fer) %>% 
  filter(treatment=="Fer") %>% 
  select(code_sp,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory)) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  geom_boxplot() +
  ggtitle("Fer") 

# Dans le natif superficiel
species_S <- ab_nat %>% 
  filter(depth == "S") %>% 
  pull(code_sp) %>% 
  unique()

MEAN_CSR %>% 
  filter(code_sp %in% species_S) %>%
  filter(treatment=="Nat") %>% 
  select(code_sp,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory)) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory)) +
  # geom_point() +
  geom_boxplot() +
  ggtitle("Nat superficiel") 


# Dans le natif
species_nat <- ab_nat %>% 
  pull(code_sp) %>% 
  unique()

MEAN_CSR %>% 
  filter(code_sp %in% species_nat) %>% 
  filter(treatment=="Nat") %>% 
  select(code_sp,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory)) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  geom_boxplot() +
  # geom_point()
  ggtitle("Nat") 


# Dans le fertile et le natif superficiel
MEAN_CSR_fer <- MEAN_CSR %>% 
  filter(code_sp %in% species_fer) %>% 
  filter(treatment=="Fer")
MEAN_CSR_natsup <- MEAN_CSR %>% 
  filter(code_sp %in% species_nat) %>% 
  filter(treatment=="Nat")
MEAN_CSR_extremes <- rbind(MEAN_CSR_fer,MEAN_CSR_natsup)

MEAN_CSR_extremes %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  # filter(score == "R") %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  # geom_point()+
  geom_boxplot() +
  facet_wrap(~treatment) 

data.anovaCSR <- MEAN_CSR_extremes %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(treatment%in% c("Nat","Fer"))

options(contrasts=c("contr.treatment","contr.treatment"))
lmCSR <- lm(R ~ treatment * LifeHistory, data = data.anovaCSR)
# plot(lmCSR)
shapiro.test(residuals(lmCSR)) # normality of residuals
lmtest::bptest(lmCSR) # homoscedasticity
lmtest::dwtest(lmCSR)  # non-independence of residuals!
anova(lmCSR)
summary(lmCSR)




