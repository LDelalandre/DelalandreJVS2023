library(tidyverse)
library(ggpubr)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv") %>%
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  mutate(log_LA = log(L_Area)) %>% 
  unique()

# I case I want to perform the analyses with species both in the trait and abundance data:
data_fer <- MEAN %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_intersect <- rbind(data_fer,data_nat)


traits <- c("LDMC","SLA","log_LA",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)#"H_FLORE","FLO_FLORE",

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
  filter(Trait %in% c("LDMC","LCC","Ldelta13C","Hrepro","Disp","SeedMass")) #"H_FLORE","FLO_FLORE",

# AJOUTER JACCARD GLOBAL (= SUR RELEVES BOTANIQUES)

#_______________________________________________________________________________
# Boxplot ####


# changer les noms des traits
# passer LA en log
# changer les noms des zones

PLOTS <- NULL
i <- 1
for (ftrait in traits){
  print(ftrait)
  ## Linear model ####
  data.anovaCSR <- MEAN %>% 
    # select(code_sp,treatment,LifeHistory,C,S,R) %>%
    # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(log_L_Area = log(L_Area))
  
  formula0 <- as.formula(paste0(ftrait, " ~ 1 + (1|code_sp)"))
  formula_treatment <- as.formula(paste0(ftrait, " ~ treatment ", " + (1|code_sp)"))
  formula_LH <- as.formula(paste0(ftrait, " ~ LifeHistory ", " + (1|code_sp)"))
  formula_no_interaction <- as.formula(paste0(ftrait, " ~ treatment + LifeHistory", " + (1|code_sp)"))
  formula <- as.formula(paste0(ftrait, " ~ treatment * LifeHistory", " + (1|code_sp)"))
  
  mmod0 <- lme4::lmer( formula0 , data = data.anovaCSR)
  
  if ( length( data.anovaCSR %>% pull(treatment) %>% unique() ) == 2 ){ # if we have data from nat and fer
    
    mmod_treatment <- lme4::lmer( formula_treatment , data = data.anovaCSR)
    pval_treatment <- car::Anova(mmod_treatment) %>% pull('Pr(>Chisq)')
    
    mmod_LH <- lme4::lmer( formula_LH , data = data.anovaCSR)
    pval_LH <- car::Anova(mmod_LH) %>% pull('Pr(>Chisq)')
    
    if(pval_treatment <= 0.25 & pval_LH <= 0.25){
      mmod_no_interaction <- lme4::lmer( formula_no_interaction , data = data.anovaCSR)
      mmod <- lme4::lmer( formula , data = data.anovaCSR) # /!\ choose fdata (includes sp just measured in on treatment)
      # or fdata2 (sp measured in both treatments only)
      
      anov0 <- anova(mmod_no_interaction,mmod0)
      print(anov0$`Pr(>Chisq)`[2])
      
      anov1 <- anova(mmod,mmod_no_interaction)
      print(anov1$`Pr(>Chisq)`[2])
      
    }else if (pval_treatment <= 0.25 & pval_LH > 0.25){
      anov <-  anova(mmod_treatment,mmod0)
      print("treatment only")
      print(anov0$`Pr(>Chisq)`[2])
    }else if (pval_treatment > 0.25 & pval_LH <= 0.25){
      anov <-  anova(mmod_LH,mmod0)
      print("LH only")
      print(anov0$`Pr(>Chisq)`[2])
    }else{
      print("no effect")
    }
    
    # 
    # formula1 <- as.formula(paste0(trait, " ~ treatment + (1|code_sp)"))
    # mmod0 <- lme4::lmer( formula0 , data = data.anovaCSR)
    # anova <- anova(mmod0,mmod)
    anov <- anova(mmod,mmod0)
    car::Anova(mmod)
    if( anov$`Pr(>Chisq)`[2]<0.05 ){
      posthoc <- multcomp::cld(emmeans::emmeans(mmod, specs = c("treatment","LifeHistory"),  type = "response",
                                                adjust = "tuckey"),
                               Letters = "abcdefghi", details = T)
      
      # contrasts = différences entre les deux traitements pour annuelles et pérennes
      posthoc$comparisons %>% 
        filter(contrast %in% c("Nat annual - Fer annual","Nat perennial - Fer perennial")) %>%
        mutate(contrast = if_else(contrast == "Nat annual - Fer annual", "Annuals","Perennials")) %>% 
        select(contrast,estimate,p.value)
      
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
    # theme(axis.title.x=element_blank())+
    # geom_boxplot(aes(color = LifeHistory)) +
    geom_boxplot()+
    geom_point(#aes(color = LifeHistory),
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
          # axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    annotate("text", x = 1, y=maxy, label = comp[1,]$.group)+ 
    annotate("text", x = 2, y=maxy, label = comp[2,]$.group) + 
    
    annotate("text", x = 1, y=miny, label = nb1)+ 
    annotate("text", x = 2, y=miny, label = nb2) +
    labs(x = "Annuals")
  # annotate("text", x = 2, y=maxy + maxy/10, label = trait) 
  
  B <- MEAN %>% 
    filter(LifeHistory == "perennial") %>%
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
    
    ggplot(aes_string(x="zone", y=ftrait, label = "code_sp",fill = "zone")) +
    theme_classic()+
    xlab("Perennials") +
    # theme(axis.title.x=element_blank())+
    # geom_boxplot(aes(color = LifeHistory)) +
    geom_boxplot() +
    geom_point(#aes(color = LifeHistory),
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
          #axis.title.x = element_blank(),
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
  
  ggplot(aes_string(x="zone", y=ftrait, label = "code_sp",fill = "zone")) + #,shape = "LifeHistory"
  theme_classic()+
  theme(axis.title.x=element_blank())+
  # geom_boxplot(aes(color = LifeHistory)) +
  geom_boxplot() +
  geom_point(#aes(color = LifeHistory,shape = LifeHistory),
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


compare_mean <- function(ftrait_1,ftrait_2,nb_bootstrap){ # compare two vectors of trait values
  # input variables are vectors of trait values for the two sets (the biggest, and the smallest)
  # (either annual and perennial species within treatment,
  # or nat and fer within lifeform)
  
  # always: 1 = ann and 2 = per (to compute per - ann)
  # or 1 = fer and 2 = nat (to compute nat - fer)
  vect_difference <- c()
  length_1 <- length(ftrait_1)
  length_2 <- length(ftrait_2)
  size = min(length_1, length_2)
  if(size == length_1){
    ftrait_small <- ftrait_1
    ftrait_big <- ftrait_2
  }else{
    ftrait_big <- ftrait_1
    ftrait_small <- ftrait_2
  }
  
  for (i in 1:nb_bootstrap){
    ftrait_subsample <- sample(ftrait_big,size = size)
    diff <- mean(ftrait_subsample) - mean(ftrait_small)
    vect_difference <- c(vect_difference,diff)
  }
  
  # I want to return 2 - 1
  if(size == length_1){ # i.e. if 1 is smaller, vect_difference is 2 - 1
    vect_difference_good_sign <- vect_difference
  }else{
    vect_difference_good_sign <- -vect_difference
  }
  
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
    
    ftrait_1 <- ftrait_ann
    ftrait_2 <- ftrait_per
    
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
    
    ftrait_1 <- ftrait_fer
    ftrait_2 <- ftrait_nat
    
  }else{
    "error, within must be either lifehistory or treatment"
  }
  output <- list(fvariable,ftrait_1,ftrait_2)
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
  ftrait_1 <- vectors_of_traits[[2]]
  ftrait_2 <- vectors_of_traits[[3]]
  
  # set.seed(100)
  vect_difference <- compare_mean(ftrait_1,ftrait_2,nb_bootstrap)
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
TABLE_LIFEFORM <- NULL
k <- 0
for (ftrait in traits){
  k <- k+1
  
  j <- which(traits == ftrait)
  LIST_DIFFERENCES <- LIST_TRAIT_DIFFERENCES[[j]] # for the focal trait
  
  diff_within_lifeform <- data.frame(
    diff = c( LIST_DIFFERENCES[[1]], LIST_DIFFERENCES[[2]] ),
    LifeForm =  c(rep("annual",nb_bootstrap), rep("perennial",nb_bootstrap) ))

  # effect size
  mean_ES = diff_within_lifeform %>% 
    group_by(LifeForm) %>% 
    summarize(mean = mean(diff)) %>% 
    mutate(mean = round(mean,digits = 3))
  
  # probability of the real difference being zero
  diff_ann <- diff_within_lifeform %>% filter(LifeForm=="annual") %>% pull(diff) 
  inf_zero_ann <- length(which(diff_ann < 0)) / length(diff_ann)
  sup_zero_ann <- length(which(diff_ann > 0)) / length(diff_ann)
  pb_include_zero_ann = min(inf_zero_ann,sup_zero_ann)
  
  ES_ann <- mean_ES %>% filter(LifeForm=="annual") %>% pull(mean)
  ES_ann_table = paste0(ES_ann," (",
                  pb_include_zero_ann,")")
  
  diff_per <- diff_within_lifeform %>% filter(LifeForm=="perennial") %>% pull(diff)
  inf_zero_per <- length(which(diff_per < 0)) / length(diff_per)
  sup_zero_per <- length(which(diff_per > 0)) / length(diff_per)
  pb_include_zero_per = min(inf_zero_per,sup_zero_per)
  
  ES_per <- mean_ES %>% filter(LifeForm=="perennial") %>% pull(mean)
  ES_per_table = paste0(ES_per ," (",
                  pb_include_zero_per,")")

  # U mann Whitney = Wilcox.test
  mean_comparison <- wilcox.test(diff ~ LifeForm, data = diff_within_lifeform)
  signif <- paste0(
  mean_comparison$statistic, " (",
  mean_comparison$p.value , ")"
  )

  table_lifeform <- data.frame(trait = ftrait,
                               ES_ann = ES_ann_table,
                               ES_per = ES_per_table,
                               # diff_ES = ES_per - ES_ann,
                               signif_comparison_ES = signif 
  )
  TABLE_LIFEFORM <- rbind(TABLE_LIFEFORM,table_lifeform)
  
  PLOT_LIFEFORM[[k]] <- diff_within_lifeform %>%
    ggplot(aes(x=diff,color = LifeForm)) +
    geom_density(size=2) +
    ggtitle(paste(ftrait)) + #,": differences within lifeform, across treatment"
    # scale_color_brewer(palette="Set2") +
    theme(legend.position = "none") +
    geom_vline(xintercept=0,size=2)
}


TABLE_LIFEFORM

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
    ggplot(aes(x=diff,color = treatment)) +
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