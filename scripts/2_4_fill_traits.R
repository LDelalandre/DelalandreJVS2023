# Take the trait data averaged by species*origin
# Complete missing trait values


library(tidyverse)

# Trait values computed in the G+F and GU(S+I) conditions, i.e. removing GUd for the latter
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int.csv")

# Trait values computed in the G+F and GU conditions
MEAN_no_subset <- read.csv2("outputs/data/mean_attribute_per_treatment.csv",encoding = "latin1") %>%
  filter(!(LifeForm1 %in% c("DPh","EPh")))

# MEAN_Nat_Dol <- read.csv2("outputs/data/mean_attribute_Nat_Dol.csv",encoding = "latin1") %>%
#   filter(!(LifeForm1 %in% c("DPh","EPh")))

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "H_FLORE","Hrepro" , # , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "FLO_FLORE", "Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)

#___________________________________________________
# Seed Mass, height and Ldelta13C ####
# D1 <- MEAN %>% 
#   filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
#   select(code_sp,Hrepro,SeedMass,Ldelta13C)
# 
# D2 <- MEAN_Nat_Dol %>% 
#   filter(LifeHistory == "perennial") %>% 
#   select(code_sp,Hrepro,SeedMass_Dol = SeedMass,Ldelta13C_Dol=Ldelta13C)
# 
# plast_across_nat <- merge(D1,D2)
# # marche bien avec hauteur et avec carbone 13, mais pas avec diamètre.
# 
# plast_across_nat %>% 
#   filter(!(SeedMass == SeedMass_Dol)) %>% 
#   ggplot(aes(x=SeedMass_Dol,y=SeedMass)) + 
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# modH <- lm(Hrepro ~ Hrepro_Dol,plast_across_nat)
# summary(modH)
# # 6 points : pas assez pour faire un modèle. Je prends juste les valeurs dans l'autre traitement.


# Add Seed Mass ####
# part of the dataset for which I have to add height
to_fill_SM <- MEAN %>% 
  filter(LifeHistory == "perennial" & treatment == "Nat" & is.na(SeedMass)) %>% 
  rename(SeedMassNA = SeedMass)
# part of the dataset for which height is ok
ok_SM <- MEAN %>% 
  filter(!(LifeHistory == "perennial" & treatment == "Nat" & is.na(SeedMass)))

to_add_SM <- MEAN_no_subset %>% 
  filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
  select(code_sp,SeedMass)

filled_SM <- left_join(to_fill_SM,to_add_SM) %>% 
  select(-SeedMassNA) %>% 
  select(all_of(colnames(MEAN))) %>% 
  unique()
MEAN_SM <- rbind(filled_SM,ok_SM) 


write.csv2(MEAN_SM,"outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv",row.names=F)


## Complete seed mass in Nat from Seed Mass in Fer ####
ann_SM_missing <- MEAN_SM %>% 
  filter(treatment=="Nat" & LifeHistory == "annual") %>% 
  filter((is.na(SeedMass))) %>% 
  pull(code_sp)

SM_to_complete <- MEAN_SM %>% 
  filter(treatment=="Fer") %>% 
  filter(code_sp %in% ann_SM_missing) %>% 
  select(code_sp,SeedMass) %>% rename(SeedMass_fer = SeedMass) %>% 
  mutate(treatment="Nat")

MEAN_SM2 <- left_join(MEAN_SM,SM_to_complete) %>%
  mutate(SeedMass = case_when(!(is.na(SeedMass_fer))~SeedMass_fer,
                              TRUE ~ SeedMass)) %>%
  select(-SeedMass_fer)

ann_SM_missing2 <- MEAN_SM2 %>% 
  filter(treatment=="Nat" & LifeHistory == "annual") %>% 
  filter((is.na(SeedMass))) %>% 
  pull(code_sp)

ann_SM_missing %>% sort()
ann_SM_missing2 %>% sort()



# Add height ####
# I complete with data from GUD
to_fill_H <- MEAN_SM2 %>% 
  filter(LifeHistory == "perennial" & treatment == "Nat" & is.na(Hrepro)) %>% 
  rename(HreproNA = Hrepro)
# part of the dataset for which height is ok
ok_H <- MEAN_SM2 %>% 
  filter(!(LifeHistory == "perennial" & treatment == "Nat" & is.na(Hrepro)))

to_add_H <- MEAN_no_subset %>% 
  filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
  select(code_sp,Hrepro)

filled_H <- left_join(to_fill_H,to_add_H) %>% 
  select(-HreproNA) %>% 
  select(all_of(colnames(MEAN))) %>% 
  unique()
MEAN_SM_H <- rbind(filled_H,ok_H) 


# Add Ldelta13C ####
# I complete with data from GUD
to_fill_13C <- MEAN_SM_H %>% 
  filter(LifeHistory == "perennial" & treatment == "Nat" & is.na(Ldelta13C)) %>% 
  rename(Ldelta13CNA = Ldelta13C)
# part of the dataset for which height is ok
ok_13C <- MEAN_SM_H %>% 
  filter(!(LifeHistory == "perennial" & treatment == "Nat" & is.na(Ldelta13C)))

to_add_13C <- MEAN_no_subset %>% 
  filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
  select(code_sp,Ldelta13C)

filled_13C <- left_join(to_fill_13C,to_add_13C) %>% 
  select(-Ldelta13CNA) %>% 
  select(all_of(colnames(MEAN))) %>% 
  unique()
MEAN_SM_H_13C <- rbind(filled_13C,ok_13C) 



write.csv2(MEAN_SM_H_13C,"outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv",row.names=F)



# how many species did we add at each step? ####
dim(MEAN %>% filter(treatment == "Nat" & !is.na(SeedMass)))[1]
dim(MEAN_SM2 %>% filter(treatment == "Nat" & !is.na(SeedMass)))[1]

dim(MEAN_SM2 %>% filter(treatment == "Nat" & !is.na(Hrepro)))[1]
dim(MEAN_SM_H %>% filter(treatment == "Nat" & !is.na(Hrepro)))[1]

dim(MEAN_SM_H %>% filter(treatment == "Nat" & !is.na(Ldelta13C)))[1]
dim(MEAN_SM_H_13C %>% filter(treatment == "Nat" & !is.na(Ldelta13C)))[1]

## From other treatment in La Fage ####
# 
# 
# ## flora ####
# traits_flore <- read.table("data/traits/flores/TRAIT_ESP_FLORE.txt",header=T) %>% 
#   select(CODE_ESP,FLO_FLORE,FRU_FLORE,H_FLORE) %>% 
#   rename(code_sp = CODE_ESP)
# 
# MEAN_completed2_flore <- left_join(MEAN_completed2,traits_flore)
# 
# write.csv2(MEAN_completed2_flore %>% unique(),
#            "outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv",
#            row.names=F)
# 
# 
# MEAN_completed2_flore %>% 
#   select(species,code_sp,treatment,LifeHistory,all_of(traits)) %>% 
#   arrange(LifeHistory,treatment,code_sp) %>%
#   select(-c(H_FLORE,FLO_FLORE)) %>% 
#   write.csv2("outputs/tables/mean-trait_values_completed_seedmass_flora.csv",row.names=F)
# 
# # check accuracy of estimation by data from flora
# MEAN_completed2_flore %>% 
#   ggplot(aes(x=Hrepro,y=H_FLORE)) + 
#   geom_point()
# 
# MEAN_completed2_flore %>% 
#   ggplot(aes(x=disp,y=FRU_FLORE)) + 
#   geom_point()
# 
# 
# # Predict phenology from those in the Fertile treatment?
# dispPF <- MEAN_completed2 %>% 
#   filter(LifeHistory == "perennial" & treatment == "Fer") %>% 
#   select(species, code_sp, Disp) %>% 
#   rename(Disp_Fer = Disp)
# 
# dispPN <- MEAN_completed2 %>% 
#   filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
#   select(species, code_sp, Disp) %>% 
#   rename(Disp_Nat = Disp)
# 
# merge(dispPF,dispPN) %>% 
#   ggplot(aes(x=Disp_Nat,y=Disp_Fer)) +
#   geom_point()
# 
# # Predict height from those in the Fertile treatment?
# hPF <- MEAN_completed2 %>% 
#   filter(LifeHistory == "perennial" & treatment == "Fer") %>% 
#   select(species, code_sp, Hrepro) %>% 
#   rename(Hrepro_Fer = Hrepro)
# 
# hPN <- MEAN_completed2 %>% 
#   filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
#   select(species, code_sp, Hrepro) %>% 
#   rename(Hrepro_Nat = Hrepro)
# 
# merge(hPF,hPN) %>% 
#   ggplot(aes(x=Hrepro_Fer,y=Hrepro_Nat)) +
#   geom_point()
# 
# # pas assez de points !
# # utiliser les flores ?
# # je fais un essai avant d'aller chercher les infos dans les flores
# predict_Hrepro <- MEAN %>% 
#   filter(LifeHistory=="perennial") %>% 
#   filter(treatment == "Nat")
#   
# predict_Hrepro %>% 
#   ggplot(aes(x=H_FLORE,y=Hrepro)) +
#   geom_point() +
#   facet_wrap(~treatment) 
#   geom_smooth(method = "lm")
# 
# mod <- lm(Hrepro ~ H_FLORE, data = predict_Hrepro)
# anova(mod) 
# summary(mod)
# plot(mod)
# 
# predict_disp <- MEAN %>% 
#   filter(LifeHistory=="perennial") %>% 
#   filter(treatment == "Fer")
# 
# mod <- lm(Disp ~ FRU_FLORE, data = predict_disp)
# anova(mod) 
# summary(mod)
# plot(mod)
# 
# predict_disp %>% 
#   ggplot(aes(x=FRU_FLORE,y=Disp)) +
#   geom_point() +
#   facet_wrap(~treatment) +
#   geom_smooth(method = "lm")
# 
# 
# 
# #___________________________________________________
# # test seed mass low plasticity
# ftrait <- "SeedMass"
# seed_fer_nat <- MEAN_no_subset %>% 
#   select(code_sp,treatment,all_of(ftrait)) %>% 
#   spread(key = treatment,value = ftrait) 
# 
# seed_fer_nat %>% 
#   ggplot(aes(x=log(Fer),y=log(Nat))) +
#   # ggplot(aes(x=Fer,y=Nat))+
#   geom_point() +
#   geom_abline(slope = 1,intercept=0) +
#   geom_smooth(method="lm")
# 
# mod <- lm(log(Nat) ~ log(Fer), data = seed_fer_nat)
# # plot(mod)
# anova(mod)
# sum <- summary(mod)
# 
# mod$coefficients
# sum$adj.r.squared
# # " OU BIEN: regarder la masse des graines comme une fonction des espèces et du traitement
# # pour montrer que la donnée espèce explique une part de variance bien plus grande
# 
# 
# 
# 
