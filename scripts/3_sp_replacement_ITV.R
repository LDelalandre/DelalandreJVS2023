library(tidyverse)
library(vegan)


# 1) Jaccard distance ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  mutate(code_sp = case_when(species == "Vicia sativa ssp. sativa"~ "VICISATI-SAT",
                             species == "Crepis vesicaria ssp. haenseleri" ~"CREPVESI-HAE",
                             species == "Taraxacum laevigatum" ~ "TARALAEV",
                             species == "Cirsium acaulon" ~ "CIRSACAU",
                             species == "Carthamus mitissimus" ~ "CARTMITI",
                             TRUE~ code_sp)) %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial"))
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial"))



## build presence matrix ####
pres_fer <- ab_fer %>% 
  mutate(treatment = "Fer") %>% 
  rename(line = id_transect_quadrat) %>%  
  mutate(presence = if_else(relat_ab >0,1,0)) %>% 
  select(code_sp,LifeHistory,treatment,line,presence)

pres_nat <- ab_nat %>% 
  mutate(treatment = "Nat") %>% 
  filter(depth == "S") %>% 
  mutate(presence = if_else(relat_ab >0,1,0)) %>%
  select(code_sp,LifeHistory,treatment,Ligne,presence) %>% 
  rename(line = Ligne)

presence <- rbind(pres_fer,pres_nat)

presence_matrix <- presence %>% 
  unique() %>% 
  spread(key = code_sp,value = presence) %>% 
  replace(is.na(.),0)

presence_matrix %>% 
  filter(treatment=="Fer") %>%
  pull(line) %>% 
  length()

presence_matrix %>% 
  filter(treatment=="Nat") %>% 
  pull(line) %>% 
  unique() %>%
  length()


## Species turnover between pairs of transects ####
JAC <- NULL
LH <- NULL
for (lh in c("annual","perennial")){
  fpres_nat <- presence_matrix %>%
    filter(treatment=="Nat") %>% 
    filter(LifeHistory == lh)
  fpres_fer <- presence_matrix %>%
    filter(treatment=="Fer") %>% 
    filter(LifeHistory == lh)
  species <- presence %>% 
    filter(LifeHistory == lh) %>% 
    pull(code_sp) %>% 
    unique()
  for (i in c(1:dim(fpres_nat[1]))){
    for(j in c(1:dim(fpres_fer[1]))){
      A <- fpres_nat[i,] %>% select(-c(LifeHistory,treatment,line))
      B <- fpres_fer[j,] %>% select(-c(LifeHistory,treatment,line))
      jac  <- vegan::vegdist(x = rbind(A,B) %>% 
                               select(all_of(species)),method="jaccard") 
      JAC <- c(JAC,jac)
      LH <- c(LH,lh)
    }
  }
}

jaccard <- data.frame(JAC,LH)
plot_jaccard <- jaccard %>% 
  mutate(LH = if_else(LH == "perennial","Perennials","Annuals")) %>% 
  ggplot(aes(x=LH,y=JAC))+
  geom_boxplot() +
  ylab("Jaccard index between pairs of transects") +
  theme_classic() +
  ggsignif::geom_signif(comparisons = list(c("Annuals", "Perennials")), 
                        map_signif_level=TRUE) +
  xlab('')

plot_jaccard

ggsave("draft/plot_jaccard.png",plot = plot_jaccard)

# mod <- lm(JAC ~ LH,data = jaccard)
# # plot(mod)
# anova(mod)
# summary(mod)

wilcox.test(JAC ~ LH,data = jaccard)





# 2) aITV ####
## Load data ####

# trait data
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv") %>% 
  dplyr::rename(LCCm = LCC) %>% 
  dplyr::rename(LNCm = LNC)
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")%>% 
  dplyr::rename(LCCm = LCC) %>% 
  dplyr::rename(LNCm = LNC)
traits <- c("LDMC","SLA","L_Area",
            "LCCm","LNCm","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

# abundance data
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(transect = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(depth == "S") %>% 
  select(-depth) %>% 
  rename(transect = Ligne)

# Join abundance datasets
ab_fer2 <- ab_fer %>% 
  select(any_of(colnames(ab_nat))) %>% 
  mutate(treatment = "Fer")
ab_nat2 <- ab_nat %>% 
  mutate(treatment = "Nat") %>% 
  select(all_of(colnames(ab_fer2)))
ab <- rbind(ab_fer2,ab_nat2) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
  select(species,code_sp,LifeForm1,LifeHistory,treatment,paddock,transect,abundance)


## Functions ####

compute_cwm <- function(ab,MEAN,MEAN_site, traits,level_of_trait_averaging){
  # This function computes CWM on all the traits.
  # ab: abundance data
  # MEAN: mean trait values at the treatment level
  # MEAN_site: mean trait values at the site level
  # traits: vector of names of traits
  # level_of-trait_averaging: either "treatment" or "site"
  
  if (level_of_trait_averaging == "treatment"){
    fMEAN <- MEAN %>% 
      select("species","code_sp","LifeForm1","LifeHistory","treatment",
             all_of(traits))
  } else if (level_of_trait_averaging == "site") {
    fMEAN <- MEAN_site %>% 
      select("species","code_sp","LifeForm1","LifeHistory",
             all_of(traits))
  }
  
  ab_traits <- ab %>% 
    left_join(fMEAN) %>% 
    group_by(transect,LifeHistory,treatment) %>% 
    mutate(relat_ab = abundance/sum(abundance)) 
  
  data_cwm <- ab_traits %>% 
    mutate_at(vars(all_of(traits)),
              .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
    rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
    unique() %>% 
    select(LifeHistory,treatment,paddock, transect,starts_with("CWM")) %>% 
    unique() %>% 
    rename_at( vars( contains( "CWM_") ), list( ~ gsub("CWM_", "", .) ) )
  
  data_cwm
}


get_cwm_decomposition <- function(ftrait,CWM,CWMfixed){
  # From the data frames of CWM and CWM fixed, generated by compute_cwm, this function returns a data frame 
  # with a row = a life history in a given transect
  # and with info on CWM, CWM fixed, and CWM intra (cf. Siefert 2015, box 1)
  # all computed for the desired trait
  
  # make a data frame with variation at these levels
  fCWM <- CWM %>% 
    select(LifeHistory,treatment,paddock,transect, all_of(ftrait)) %>% 
    mutate(CWM = get(ftrait)) %>% 
    select(-ftrait)
  fCWMfixed <- CWMfixed %>% 
    select(LifeHistory,treatment,paddock,transect, all_of(ftrait)) %>% 
    mutate(CWMfixed = get(ftrait)) %>% 
    select(-ftrait)
  
  # CWM data for the focal trait, at the three levels
  fCWMall <- full_join(fCWM,fCWMfixed) %>% 
    mutate(CWMintra = CWM - CWMfixed)
  fCWMall
}


compute_SS <- function(data,LH){ 
  # data is a data frame with a row = a life history in a given transect
  # and with info on CWM, CWM fixed, and CWM intra (cf. Siefert 2015, box 1)
  # it is the output of get_cwm_decomposition
  data_LH <- data %>% 
    filter(LifeHistory == LH)
  
  modtot <- lm(CWM ~ 1, data = data_LH)
  anotot <- anova(modtot)
  SStot <- anotot$`Sum Sq`
  
  modfixed <- lm(CWMfixed ~ 1, data = data_LH)
  anofixed <- anova(modfixed)
  SSfixed <- anofixed$`Sum Sq`
  
  modintra <- lm(CWMintra ~ 1, data = data_LH)
  anointra <- anova(modintra)
  SSintra <- anointra$`Sum Sq`
  
  SScov <- 100*(SStot - SSintra - SSfixed)/SStot # BIZARRE, CETTE FORMULE NE DONNE PAS somme(SS) = SStot...
  
  data.frame(SStot = SStot,SSfixed = SSfixed,SSintra = SSintra ,SScov = SScov)
}



## Analyses ####

# CWM with both species replacement and intrasp variability  
CWM <- compute_cwm(ab,MEAN,MEAN_site,traits,level_of_trait_averaging = "treatment")
# CWM with species replacement but no intraspecific variability
CWMfixed <- compute_cwm(ab,MEAN,MEAN_site,traits,level_of_trait_averaging = "site")
# CWM intra  = CWM - CWMfixed pour chaque trait
# 

AITV <- NULL
TRAIT <- NULL
LH <- NULL
SUMSQ <- NULL
for (ftrait in traits){
  # compute the 3 CWM on this trait, for each LH in each transect
  fCWMall <- get_cwm_decomposition(ftrait,CWM,CWMfixed) 
  for (fLH in c("annual","perennial")){
    SumSq <- compute_SS(fCWMall,fLH) %>% 
      mutate(trait = ftrait,LifeHistory = fLH)
    # compute aITV for each LH
    aITV <-  SumSq %>% 
      summarize(aITV = log(SSintra/SSfixed)) %>% 
      pull(aITV)
    AITV <- c(AITV,aITV)
    TRAIT <- c(TRAIT,ftrait)
    LH <- c(LH,fLH)
    SUMSQ <- rbind(SUMSQ,SumSq) 
  }
}

## aITV ####
data_aITV <- data.frame(aITV = AITV, trait = TRAIT,LifeHistory = LH) %>% 
  spread(key = LifeHistory,value=aITV) %>% 
  arrange(factor(trait),levels = traits) %>% 
  filter(!(trait == "SeedMass")) 
# filter(!(trait %in% c("Hrepro","SeedMass","Ldelta13C")))
# If aITV inferior to zero: greater contribution of species turnover to total among-community variance


boxplot_aITV <- data_aITV %>% 
  rename(Annuals = annual, Perennials = perennial) %>% 
  gather(key = LifeHistory, value = aITV, - trait) %>% 
  ggplot(aes(x = LifeHistory,y = aITV)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 0,linetype = "dashed") +
  ggsignif::geom_signif(comparisons = list(c("Annuals", "Perennials")), 
                        map_signif_level=TRUE) +
  xlab("") +
  ylab("Contribution of ITV to trait variance (aITV)")




  

# stats sur les boxplot:
# - aITV différent de zéro ?
# - différent entre annuelles et pérennes ?
mod_aITV <- lm(aITV ~ LifeHistory, data = data_aITV %>% 
                 gather(key = LifeHistory, value = aITV, - trait))
anova(mod_aITV) # METTRE LA SORTIE D'ANOVA EN LEGENDE ?
summary(mod_aITV)



# plots_aITV <- gridExtra::grid.arrange(boxplot_aITV,plot_aITV, ncol=2)


library(cowplot)

plot_aITV <- data_aITV %>% 
  mutate(trait = case_when(trait == "L_Area" ~ "LA",
                           trait == "Hrepro" ~ "RPH",
                           trait == "Ldelta13C" ~ "Lδ13C ",
                           TRUE ~ trait)) %>% 
  # filter(!(trait %in% c("L_Area"))) %>%
  ggplot(aes(x=perennial,y=annual,label = trait))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0,linetype = "dashed")+
  geom_vline(xintercept = 0,linetype = "dashed")+
  xlim(c(-9,2)) +
  ylim(c(-9,2)) +
  ggrepel::geom_label_repel() +
  xlab("aITV of perennials") +
  ylab("aITV of annuals") +
  # theme_void() +
  theme_classic() +
  theme(axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  theme(text = element_text(size = 15))

boxplot_aITV_ann <- data_aITV %>% 
  rename(Annuals = annual, Perennials = perennial) %>% 
  gather(key = LifeHistory, value = aITV, - trait) %>% 
  filter(LifeHistory == "Annuals") %>% 
  ggplot(aes(x = LifeHistory,y = aITV)) +
  geom_boxplot() +
  theme_classic()+
  # theme_void() +
  ylim(c(-9,2)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank()) +
  ylab("aITV of annuals")+ 
  ggeasy::easy_remove_y_axis(what = c("ticks", "text", "line"), teach = FALSE) +
  ggeasy::easy_remove_x_axis(what = c("ticks", "text", "line"), teach = FALSE) +
  theme(text = element_text(size = 15))

boxplot_aITV_per <- data_aITV %>% 
  rename(Annuals = annual, Perennials = perennial) %>% 
  gather(key = LifeHistory, value = aITV, - trait) %>% 
  filter(LifeHistory == "Perennials") %>% 
  ggplot(aes(y = LifeHistory,x = aITV)) +
  geom_boxplot() +
  theme_classic()+
  # theme_void() +
  xlim(c(-9,2)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.y=element_blank())+
  xlab("aITV of perennials")+ 
  ggeasy::easy_remove_x_axis(what = c("ticks", "text", "line"), teach = FALSE) +
  ggeasy::easy_remove_y_axis(what = c("ticks", "text", "line"), teach = FALSE) +
  theme(text = element_text(size = 15))



plots_aITV <-ggdraw() +
  draw_plot(plot_aITV, x = 0.2,y=0.2,width = 0.8, height = 0.8) +
  draw_plot(boxplot_aITV_ann, x = 0, y = 0.23, width = .2, height = .8) +
  draw_plot(boxplot_aITV_per, x = 0.237, y = 0, width = .8, height = .2)

plots_aITV

ggsave("draft/plot_aITV.jpg",plots_aITV,width = 8, height = 8)



