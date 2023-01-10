library(tidyverse)
library(ggpubr)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole"))

# I case I want to perform the analyses with species both in the trait and abundance data:
data_fer <- MEAN %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_intersect <- rbind(data_fer,data_nat)


traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "H_FLORE","Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "FLO_FLORE", "Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)

# Chose to take species in trait data, or in both trait and abundance data


#_______________________________________________________________________________
# Table: Jaccard and species number ####
ftrait <- "LDMC"

TABLE <- NULL

for (ftrait in  traits){
  
  # species present for the focal trait
  MEAN_jaccard <- MEAN %>% 
    select(code_sp,LifeForm1,treatment,LifeHistory,ftrait) %>% 
    filter(!(is.na(get(ftrait))))
  
  line_table <- data.frame(Trait=ftrait,Annuals = NA,Perennials = NA)
  
  for (LH in c("annual","perennial")){
    # nb of species of the given life form in both treatments
    nb_sp_fer <-  MEAN_jaccard %>% 
      filter(LifeHistory == LH & treatment == "Fer") %>% 
      pull(code_sp) %>% 
      length()
    
    nb_sp_nat <-  MEAN_jaccard %>% 
      filter(LifeHistory == LH & treatment == "Nat") %>% 
      pull(code_sp) %>% 
      length()
    
    # data in the matrix form (presence-absence)
    MEAN_jaccard_LH <- MEAN_jaccard %>% 
      filter(LifeHistory == LH) %>%
      select(code_sp,treatment) %>% 
      mutate(presence = 1) %>% 
      spread(key = code_sp,value = presence) %>% 
      column_to_rownames("treatment") %>% 
      replace(is.na(.), 0)
    
    Jac <- vegan::vegdist(x = MEAN_jaccard_LH,method="jaccard") %>% 
      round(digits = 2)
    info <- paste0("(",
                   paste(nb_sp_fer,nb_sp_nat,Jac, sep = ";"),
                   ")" )
    
    if(LH == "annual"){
      line_table$Annuals <- info
    } else{
      line_table$Perennials <- info
    }
    
  }
  TABLE <- rbind(TABLE,line_table)
}

TABLE %>% 
  filter(Trait %in% c("LDMC","LCC","Ldelta13C","H_FLORE","FLO_FLORE","SeedMass"))

# AJOUTER JACCARD GLOBAL (= SUR RELEVES BOTANIQUES)

#_______________________________________________________________________________
# Boxplot ####

# example of oneplot
ftrait <- "LDMC"
MEAN %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  ggplot(aes_string(x="LifeHistory", y=ftrait, label = "code_sp",fill = "zone")) +
  theme_classic()+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point(aes(color = LifeHistory),position = position_dodge(width = .75)) +
  scale_fill_manual(values = c("grey", "white")) +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75))
# extraire la légende !


PLOTS <- NULL
i <- 1
for (ftrait in traits){
  
  ## Linear model ####
  data.anovaCSR <- MEAN %>% 
    # select(code_sp,treatment,LifeHistory,C,S,R) %>%
    # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(log_L_Area = log(L_Area))
  
  
  formula <- as.formula(paste0(ftrait, " ~ treatment * LifeHistory", " + (1|code_sp)"))
  formula0 <- as.formula(paste0(ftrait, " ~ 1 + (1|code_sp)"))
  
  if ( length( data.anovaCSR %>% pull(treatment) %>% unique() ) == 2 ){ # if we have data from nat and fer
    mmod <- lme4::lmer( formula , data = data.anovaCSR) # /!\ choose fdata (includes sp just measured in on treatment)
    # or fdata2 (sp measured in both treatments only)
    
    formula0 <- as.formula(paste0(ftrait, " ~ 1 + (1|code_sp)"))
    mmod0 <- lme4::lmer( formula0 , data = data.anovaCSR)
    # 
    # formula1 <- as.formula(paste0(trait, " ~ treatment + (1|code_sp)"))
    # mmod0 <- lme4::lmer( formula0 , data = data.anovaCSR)
    # anova <- anova(mmod0,mmod)
    anov <- anova(mmod,mmod0)
    if( anov$`Pr(>Chisq)`[2]<0.05 ){
      posthoc <- multcomp::cld(emmeans::emmeans(mmod, specs = c("treatment","LifeHistory"),  type = "response",
                                                adjust = "tuckey"),
                               Letters = "abcdefghi", details = T)
      
      comp <- as.data.frame(posthoc$emmeans) %>% 
        mutate(treatment = factor(treatment,levels = c("Fer","Nat"))) %>% 
        mutate(LifeHistory = factor(LifeHistory, levels = c("annual","perennial"))) %>% 
        arrange(LifeHistory,treatment)
    }
    
  }
  
  ## plot ####
  miny <- min( MEAN %>% pull(sym(ftrait)) , na.rm = T)
  maxydata <- max(MEAN %>% pull(sym(ftrait)), na.rm = T)
  
  maxy <- maxydata + (maxydata - miny)/10
  
  # nb of species
  nb1 <- MEAN %>% 
    filter(LifeHistory == "annual") %>%
    filter(treatment =="Fer") %>% 
    pull(sym(ftrait)) %>% 
    na.omit() %>% 
    length()
  nb2 <- MEAN %>% 
    filter(LifeHistory == "annual") %>%
    filter(treatment =="Nat") %>% 
    pull(sym(ftrait)) %>% 
    na.omit() %>% 
    length()
  nb3 <- MEAN %>% 
    filter(LifeHistory == "perennial") %>%
    filter(treatment =="Fer") %>% 
    pull(sym(ftrait)) %>% 
    na.omit() %>% 
    length()
  nb4 <- MEAN %>% 
    filter(LifeHistory == "perennial") %>%
    filter(treatment =="Nat") %>% 
    pull(sym(ftrait)) %>% 
    na.omit() %>% 
    length()
  
  A <- MEAN %>% 
    filter(LifeHistory == "annual") %>%
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
    
    ggplot(aes_string(x="zone", y=ftrait, label = "code_sp",fill = "zone")) +
    theme_classic()+
    theme(axis.title.x=element_blank())+
    # geom_boxplot(aes(color = LifeHistory)) +
    geom_boxplot()+
    geom_point(aes(color = LifeHistory),
               shape = 19,size = 2,
               position = position_dodge(width = .75)) +
    scale_fill_manual(values = c("grey", "white")) +
    # scale_fill_manual(values = c("grey30", "grey80")) +
    geom_line(aes(group = code_sp),
              alpha=0.4) + #,color=LifeHistory
    ylim(c(miny,maxy)) +
    theme(legend.position="none") +
    # ggtitle ("annuals")
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() ,
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    annotate("text", x = 1, y=maxy, label = comp[1,]$.group)+ 
    annotate("text", x = 2, y=maxy, label = comp[2,]$.group) + 
    
    annotate("text", x = 1, y=miny, label = nb1)+ 
    annotate("text", x = 2, y=miny, label = nb2)
  # annotate("text", x = 2, y=maxy + maxy/10, label = trait) 
  
  B <- MEAN %>% 
    filter(LifeHistory == "perennial") %>%
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
    
    ggplot(aes_string(x="zone", y=ftrait, label = "code_sp",fill = "zone")) +
    theme_classic()+
    theme(axis.title.x=element_blank())+
    # geom_boxplot(aes(color = LifeHistory)) +
    geom_boxplot() +
    geom_point(aes(color = LifeHistory),
               shape = 17, size = 2,
               position = position_dodge(width = .75)) +
    # scale_fill_manual(values = c("grey30", "grey80")) +
    scale_fill_manual(values = c("grey", "white")) +
    geom_line(aes(group = code_sp),
              alpha = 0.4) + #,color=LifeHistory
    scale_color_manual(values = "#00BFC4") +
    ylim(c(miny,maxy))+
    theme(legend.position="none") +
    ggeasy::easy_remove_y_axis() +
    # ggtitle("perennials") 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() ,
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    )  + 
    annotate("text", x = 1, y=maxy, label = comp[3,]$.group)+ 
    annotate("text", x = 2, y=maxy, label = comp[4,]$.group)+ 
    
    annotate("text", x = 1, y=miny, label = nb3)+ 
    annotate("text", x = 2, y=miny, label = nb4)
  
  
  
  # Create a text grob (for the title = trait name)
  tgrob <- ggpubr::text_grob(ftrait,size = 10)
  # Draw the text
  plot_0 <- ggpubr::as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
  
  plot <- ggpubr::ggarrange(plot_0,NULL,A, B,
                            ncol = 2,nrow = 2,heights = c(1,5))
  PLOTS[[i]] <- plot
  i <- i+1
}




# Extract the legend alone, from the data frame of species removal expe
plot <- MEAN %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  
  ggplot(aes_string(x="zone", y=ftrait, label = "code_sp",fill = "zone",shape = "LifeHistory")) +
  theme_classic()+
  theme(axis.title.x=element_blank())+
  # geom_boxplot(aes(color = LifeHistory)) +
  geom_boxplot() +
  geom_point(aes(color = LifeHistory,shape = LifeHistory),
             size = 2,
             position = position_dodge(width = .75)) +
  # scale_fill_manual(values = c("grey30", "grey80")) +
  scale_fill_manual(values = c("grey", "white")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ,
        axis.title.x = element_blank()
  )

leg <- ggpubr::get_legend(plot)
legend <- ggpubr::as_ggplot(leg)


boxplot_all_traits <- ggpubr::ggarrange(PLOTS[[1]],PLOTS[[2]],PLOTS[[3]],PLOTS[[4]],PLOTS[[5]],PLOTS[[6]],
                                        PLOTS[[7]],PLOTS[[8]],PLOTS[[9]],
                                        PLOTS[[10]],PLOTS[[11]],PLOTS[[12]],legend)

ggsave("draft/boxplot_all_traits.jpg",boxplot_all_traits,width = 14, height = 8)




#_______________________________________________________________________________
# Pairwise comparisons ####
# Subsample to have the same number of species (bootstrap)

# Between lifehistory within treatment or between treatment within lifehistory


compare_mean <- function(ftrait_big,ftrait_small,nb_bootstrap){ # compare two vectors of trait values
  # input variables are vectors of trait values for the two sets (the biggest, and the smallest)
  # (either annual and perennial species within treatment,
  # or nat and fer within lifeform)
  vect_difference <- c()
  size = length(ftrait_small)
  for (i in 1:nb_bootstrap){
    ftrait_subsample <- sample(ftrait_big,size = size)
    diff <- mean(ftrait_subsample) - mean(ftrait_small)
    vect_difference <- c(vect_difference,diff)
  }
  vect_difference
}

get_vectors_of_traits <- function(within,MEAN,fvariable){
  # within is either "lifehistory" or "treatment"
  # MEAN is the data frame of trait values
  # "fvariable" is the value of ftreatment or flifehistory, depending on the value of the "within" variable
  # so it is either "annual" or perennial in one case, or "Nat" or "Fer" in the other
  if(within=="treatment"){
    ftreatment <- fvariable
    
    # trait values in the two categories (here, treatment)
    MEAN_ftrait <- MEAN %>% 
      filter(treatment == ftreatment) %>% 
      select(code_sp,LifeForm1,treatment,LifeHistory,ftrait) %>% 
      spread(key = LifeHistory, value = get(ftrait))
    
    # vectors of trait values for annuals and perennials
    ftrait_ann <- MEAN_ftrait %>% 
      filter(!is.na(annual )) %>% 
      pull(annual)
    
    ftrait_per <- MEAN_ftrait %>% 
      filter(!is.na(perennial )) %>% 
      pull(perennial)
    
    # subsample the longest trait vector to have the same number of trait values for annuals and perennials
    size = min(length(ftrait_ann),length(ftrait_per))
    if (size == length(ftrait_ann)){
      ftrait_big <- ftrait_per
      ftrait_small <- ftrait_ann
      order <- c("big = perennial","small = annual")
    } else {
      ftrait_big <- ftrait_ann
      ftrait_small <- ftrait_per
      order <- c("big = annual","small = perennial")
    }
  }else if(within == "lifehistory") {
    flifehistory <- fvariable
    
    # trait values in the two categories (here, lifeforms)
    MEAN_ftrait <- MEAN %>% 
      filter(LifeHistory == flifehistory) %>% 
      select(code_sp,treatment,ftrait) %>% 
      spread(key = treatment, value = get(ftrait))
    
    # vectors of trait values for nat and fer
    ftrait_fer <- MEAN_ftrait %>% 
      filter(!is.na(Fer )) %>% 
      pull(Fer)
    
    ftrait_nat <- MEAN_ftrait %>% 
      filter(!is.na(Nat )) %>% 
      pull(Nat)
    
    # subsample the longest trait vector to have the same number of trait values for annuals and perennials
    size = min(length(ftrait_fer),length(ftrait_nat))
    if (size == length(ftrait_fer)){
      ftrait_big <- ftrait_nat
      ftrait_small <- ftrait_fer
      order <- c("big = Nat","small = Fer")
    } else {
      ftrait_big <- ftrait_fer
      ftrait_small <- ftrait_nat
      order <- c("big = Fer","small = Nat")
    }
  }else{
    "error, within must be either lifehistory or treatment"
  }
  output <- list(fvariable,order,ftrait_big,ftrait_small)
  output # output is the absolute value of difference
}



get_distribution_of_differences <- function(within,MEAN,fvariable){
  # if(within=="treatment"){
  #   fvariable <- ftreatment
  # }else{
  #   fvariable <- flifehistory
  # }
  vectors_of_traits <- get_vectors_of_traits(within,MEAN,fvariable)
  # vectors_of_traits[[1]]
  # vectors_of_traits[[2]]
  ftrait_big <- vectors_of_traits[[3]]
  ftrait_small <- vectors_of_traits[[4]]
  
  # set.seed(100)
  vect_difference <- compare_mean(ftrait_big,ftrait_small,nb_bootstrap)
  vect_difference
}


ftrait <- "LDMC" # choose the trait
nb_bootstrap = 1000






# have the distribution of difference for the two groups of traits with subsampling
WITHIN <- c("lifehistory","treatment")
LIFEHISTORY <-  c("annual","perennial")
TREATMENT <-  c("Fer","Nat")


within <- WITHIN[1]
flifehistory <-  LIFEHISTORY[1] # if I compare treatments
ftreatment <- TREATMENT[2] # if I compare lifehistory

LIST_TRAIT_DIFFERENCES <- list() # list of the four comparisons, for each trait
j <- 0
for (ftrait in traits){
  j <- j+1
  LIST_DIFFERENCES <- list() 
  # four comparisons (within lifeform across treatments, and within treatment across lifeforms)
  i <- 0
  for (within in WITHIN){
    if(within == "lifehistory"){
      # compare:
      # annuals in fer and nat
      # perennials in fer and nat
      for (flifehistory in LIFEHISTORY){
        i <- i+1
        distrib <- get_distribution_of_differences(within,MEAN,flifehistory)
        LIST_DIFFERENCES[[i]] <- distrib
      }
      
    }else if(within == "treatment"){
      for (ftreatment in TREATMENT){
        # compare:
        # annuals and perennials in fer
        # annuals and perennials in nat
        i <- i+1
        distrib <- get_distribution_of_differences(within,MEAN,ftreatment)
        LIST_DIFFERENCES[[i]] <- distrib
      }
    }
    
  }
  LIST_TRAIT_DIFFERENCES[[j]] <- LIST_DIFFERENCES
}



# NB: comparaison du LDMC des pérennes dans le fertile et le natif:
# même nombre d'espèces ! Du coup, le bootstrap ne fonctionne pas...



## plot differences within lifeform ####
PLOT_LIFEFORM <- list()
k <- 0
for (ftrait in traits){
  k <- k+1
  
  j <- which(traits == ftrait)
  LIST_DIFFERENCES <- LIST_TRAIT_DIFFERENCES[[j]] # for the focal trait
  
  diff_within_lifeform <- data.frame(
    diff = c( LIST_DIFFERENCES[[1]], LIST_DIFFERENCES[[2]] ),
    LifeForm =  c(rep("annual",nb_bootstrap), rep("perennial",nb_bootstrap) ))
  
  
  PLOT_LIFEFORM[[k]] <- diff_within_lifeform %>%
    mutate(abs_diff = abs(diff)) %>% 
    ggplot(aes(x=abs_diff,color = LifeForm)) +
    geom_density(size=2) +
    ggtitle(paste(ftrait)) + #,": differences within lifeform, across treatment"
    # scale_color_brewer(palette="Set2") +
    theme(legend.position = "none")
}

LFplot <- diff_within_lifeform %>%
  ggplot(aes(x=diff,color = LifeForm)) +
  geom_density(size=2)
LFleg <- ggpubr::get_legend(LFplot)
LFlegend <- ggpubr::as_ggplot(LFleg)

density_lifeform <- ggpubr::ggarrange(PLOT_LIFEFORM[[1]],PLOT_LIFEFORM[[2]],PLOT_LIFEFORM[[3]],
                                      PLOT_LIFEFORM[[4]],PLOT_LIFEFORM[[5]],PLOT_LIFEFORM[[6]],
                                      PLOT_LIFEFORM[[7]],PLOT_LIFEFORM[[8]],PLOT_LIFEFORM[[9]],
                                      PLOT_LIFEFORM[[10]],PLOT_LIFEFORM[[11]],PLOT_LIFEFORM[[12]],
                                      LFlegend)
final_plot_lifeform <- 
  ggpubr::annotate_figure(density_lifeform,
                          top = text_grob("Trait differences within lifeform, across treatment ", 
                          color = "black", face = "bold", size = 14))

final_plot_lifeform




    



## plot differences within lifeform ####
PLOT_TREATMENT <- list()
k <- 0
for (ftrait in traits){
  k <- k+1
  
  j <- which(traits == ftrait)
  LIST_DIFFERENCES <- LIST_TRAIT_DIFFERENCES[[j]] # for the focal trait
  
  diff_within_treatment <- data.frame( 
    diff = c(LIST_DIFFERENCES[[3]],LIST_DIFFERENCES[[4]]),
    treatment = c(rep("Fer",nb_bootstrap), rep("Nat",nb_bootstrap)))
  
  
  PLOT_TREATMENT[[k]] <- diff_within_treatment %>% 
    mutate(abs_diff = abs(diff)) %>%
    ggplot(aes(x=abs_diff,color = treatment)) +
    geom_density(size = 2) +
    ggtitle(paste(ftrait))  +
    scale_color_brewer(palette="Set2")+
    theme(legend.position = "none")
}

TRplot <- diff_within_treatment %>% 
  ggplot(aes(x=diff,color = treatment)) +
  geom_density(size = 2)  +
  scale_color_brewer(palette="Set2")
TRleg <- ggpubr::get_legend(TRplot)
TRlegend <- ggpubr::as_ggplot(TRleg)

density_treatment <- ggpubr::ggarrange(PLOT_TREATMENT[[1]],PLOT_TREATMENT[[2]],PLOT_TREATMENT[[3]],
                                       PLOT_TREATMENT[[4]],PLOT_TREATMENT[[5]],PLOT_TREATMENT[[6]],
                                       PLOT_TREATMENT[[7]],PLOT_TREATMENT[[8]],PLOT_TREATMENT[[9]],
                                       PLOT_TREATMENT[[10]],PLOT_TREATMENT[[11]],PLOT_TREATMENT[[12]],
                                      TRlegend)

final_plot_treatment <- ggpubr::annotate_figure(density_treatment,
                top = text_grob("Trait differences within treatment, across lifeform ", 
                                color = "black", face = "bold", size = 14))

final_plot_treatment

# diff is the effect size (correct it for small sample sizes, cf. garnier 2019 (Flora)?)




# Trait correlation within lifeform ####