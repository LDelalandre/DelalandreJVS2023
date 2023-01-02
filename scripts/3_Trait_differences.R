library(tidyverse)

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
            "H_FLORE",#"Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "FLO_FLORE", #Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)

# Chose to take species in trait data, or in both trait and abundance data



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
  scale_fill_manual(values = c("grey", "white"))
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
    annotate("text", x = 2, y=maxy, label = comp[2,]$.group) 
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
    annotate("text", x = 2, y=maxy, label = comp[4,]$.group)
  
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
  
  ggplot(aes_string(x="zone", y=trait, label = "code_sp",fill = "zone",shape = "LifeHistory")) +
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
                                        PLOTS[[7]],PLOTS[[8]],PLOTS[[9]],legend)

ggsave("draft/boxplot_all_traits.jpg",boxplot_all_traits,width = 14, height = 8)





# Pairwise comparisons ####
# Subsample to have the same number of species (bootstrap)

## Between lifeforms within treatment ####

## Between treatments within lifeform ####



# Trait correlation within lifeform ####