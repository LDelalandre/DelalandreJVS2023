library(tidyverse)
library(ggpubr)
library(kableExtra)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv") %>%
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  mutate(log_LA = log(L_Area)) %>% 
  mutate(log_SeedMass = log(SeedMass)) %>% 
  unique()

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

# I case I want to perform the analyses with species both in the trait and abundance data:
data_fer <- MEAN %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_intersect <- rbind(data_fer,data_nat)


traits <- c("LDMC","SLA","log_LA",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro" ,  #, "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "Disp",#"Mat_Per", #"Mat","Flo",
            "log_SeedMass"
)#"H_FLORE","FLO_FLORE",

traits_names <- c("Leaf Dry Matter Content (mg/g)", "Specific Leaf Area (m²/kg)"," log(Leaf Area (cm²))",
                  "Leaf Carbon Content/mass (mg/g)","Leaf Nitrogen Content/mass (mg/g)", 
                  paste("Leaf δ13C (part per thousand)"),
                  "Reproductive Height (cm)", 
                  "Date of first dispersal (Julian day)",
                  "log(Seed Mass (mg))")


# Chose to take species in trait data, or in both trait and abundance data


#_______________________________________________________________________________
# Table: Jaccard and species number ####

AF <- ab_fer %>% 
  filter(LifeForm1=="The") %>% 
  select(species) %>% 
  unique() %>% 
  mutate(presence_fer=1)
nb_ann_fer <- dim(AF)[1]
AN <- ab_nat %>% 
  filter(LifeForm1=="The") %>% 
  select(species) %>% 
  unique() %>% 
  mutate(presence_nat=1)
nb_ann_nat <- dim(AN)[1]

table_ann <- full_join(AF,AN) %>% 
  replace(is.na(.),0) %>%
  column_to_rownames("species") %>% 
  as.matrix() %>% 
  t()

Jac_ann <- vegan::vegdist(x = table_ann,method="jaccard") %>% 
  round(digits = 2)
info_ann <- paste0(Jac_ann," (",
                   paste(nb_ann_fer,nb_ann_nat, sep = ";"),
                   ")" )

# Perennials in Fertile
PF <- ab_fer %>% 
  filter(!(LifeForm1=="The")) %>% 
  select(species) %>% 
  unique() %>% 
  mutate(presence_fer=1)
nb_per_fer <- dim(PF)[1]
PN <- ab_nat %>% 
  filter(!(LifeForm1=="The")) %>% 
  select(species) %>% 
  unique() %>% 
  mutate(presence_nat=1)
nb_per_nat <- dim(PN)[1]

full_join(PF,PN) %>% 
  replace(is.na(.),0) %>%
  column_to_rownames("species") %>% 
  as.matrix() %>% 
  t()

table_per <- full_join(PF,PN) %>% 
  replace(is.na(.),0) %>%
  column_to_rownames("species") %>% 
  as.matrix() %>% 
  t()



Jac_per <- vegan::vegdist(x = table_per,method="jaccard") %>% 
  round(digits = 2)
info_per <- paste0(Jac_per," (",
                   paste(nb_per_fer,nb_per_nat, sep = ";"),
                   ")" )



TABLE <- data.frame(Trait = "Releves",
                    Annuals = info_ann,
                    Perennials = info_per)
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
    info <- paste0(Jac," (",
                   paste(nb_sp_fer,nb_sp_nat, sep = ";"),
                   ")" )
    
    if(LH == "annual"){
      line_table$Annuals <- info
    } else{
      line_table$Perennials <- info
    }
    
  }
  TABLE <- rbind(TABLE,line_table)
}

table_sp_replacement <- TABLE %>% 
  # filter(Trait %in% c("LDMC","LCC","Ldelta13C","Hrepro","Disp","SeedMass")) %>% #"H_FLORE","FLO_FLORE",
  kableExtra::kable( escape = F,
                     col.names = c("Origin", "Annuals", "Perennials")) %>%
  kableExtra::kable_styling("hover", full_width = F)

cat(table_sp_replacement, file = "draft/table_species_replacement.doc")

# AJOUTER JACCARD GLOBAL (= SUR RELEVES BOTANIQUES)

# Faire plutôt un indice de jaccard par paire de transect

#_______________________________________________________________________________
# Interspecific comparisons ####

species_experiment <- read.csv2("data/annual_species_experiment.csv") %>% 
  pull(code_sp)

# changer les noms des traits
# passer LA en log
# changer les noms des zones

MEAN <- MEAN %>% 
  filter((code_sp %in% species_experiment))

TABLE_PVAL <- NULL
PLOTS <- NULL
i <- 1
for (ftrait in traits){
  ## Linear model ####
  data.anovaCSR <- MEAN %>% 
    # select(code_sp,treatment,LifeHistory,C,S,R) %>%
    # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(log_L_Area = log(L_Area)) %>% 
    mutate(treatment = as.factor(treatment) ) %>% 
    mutate(LifeHistory = as.factor(LifeHistory))
  
  formula0 <- as.formula(paste0(ftrait, " ~ 1 + (1|code_sp)"))
  formula <- as.formula(paste0(ftrait, " ~ treatment * LifeHistory", " + (1|code_sp)"))
  
  mmod0 <- lme4::lmer( formula0 , data = data.anovaCSR)
  
  if ( length( data.anovaCSR %>% pull(treatment) %>% unique() ) == 2 ){ # if we have data from nat and fer
    mmod <- lme4::lmer( formula , data = data.anovaCSR) # /!\ choose fdata (includes sp just measured in on treatment)
    # or fdata2 (sp measured in both treatments only)
    
    # formula1 <- as.formula(paste0(trait, " ~ treatment + (1|code_sp)"))
    # mmod0 <- lme4::lmer( formula0 , data = data.anovaCSR)
    # anova <- anova(mmod0,mmod)
    # anov <- anova(mmod,mmod0)
    
    anov <- car::Anova(mmod)
    pval <- anov %>%
      as.data.frame() %>% 
      rename(pval = 'Pr(>Chisq)') %>%
      pull(pval) %>% 
      format(scientific = TRUE, digits = 2) %>% 
      as.numeric()
    table_pval <-  data.frame(Trait = ftrait,
                              Treatment = pval[1],
                              Life_History = pval[2],
                              Interaction = pval[3]
    )
    
    # part of variance explained
    variance <- MuMIn::r.squaredGLMM(mmod) 
    vfixed <- variance[1] # marginal R squared value associated with fixed effects
    vcomplete <- variance[2] # conditional R2 value associated with fixed effects plus the random effects. 
    part_var_species <- (vcomplete-vfixed/vcomplete)

    # if( table_pval$Interaction < 0.05 ){
    
    # pb version package: la version récente de emmeans demande 
    # remove.packages("emmeans")
    # library(remotes)
    # install_version("emmeans", "1.5.2-1")
    EM <- emmeans::emmeans(mmod, specs = c("treatment","LifeHistory"),  type = "response",
                           adjust = "tuckey")
    # emmeans (version 1.5.2-1)
    posthoc <- multcomp::cld(EM, #emmeans::as.glht(EM)
                             Letters = "abcdefghi", details = T)
    
    # contrasts = différences entre les deux traitements pour annuelles et pérennes
    contrasts <- posthoc$comparisons %>% 
      filter(contrast %in% c("Nat annual - Fer annual","Nat perennial - Fer perennial",
                             "Fer annual - Nat annual","Fer perennial - Nat perennial")) %>%
      mutate(estimate = case_when(contrast %in% c("Fer annual - Nat annual","Fer perennial - Nat perennial") ~ -estimate,
                                  TRUE ~ estimate)) %>% 
      mutate(contrast = case_when(contrast == "Fer annual - Nat annual" ~ "Nat annual - Fer annual",
                                  contrast == "Fer perennial - Nat perennial" ~ "Nat perennial - Fer perennial",
                                  TRUE ~ contrast)) %>% 
      mutate(contrast = if_else(contrast == "Nat annual - Fer annual", "Annuals","Perennials")) %>% 
      select(contrast,estimate,p.value) %>% 
      mutate(estimate = round(estimate,digits =1)) %>% 
      mutate(p.value = format(p.value,scientific = TRUE, digits = 2)) %>% 
      mutate(est_pval = paste0(estimate, " (",p.value,")")) %>%
      select(contrast,est_pval) %>% 
      spread(key = contrast,value = est_pval)
      
      table_pval2 <- cbind(table_pval,contrasts)
      table_pval2$percent_var_code_sp <- part_var_species
      TABLE_PVAL <- rbind(TABLE_PVAL,table_pval2)
      
      comp <- as.data.frame(posthoc$emmeans) %>% 
        mutate(treatment = factor(treatment,levels = c("Fer","Nat"))) %>% 
        mutate(LifeHistory = factor(LifeHistory, levels = c("annual","perennial"))) %>% 
        arrange(LifeHistory,treatment)
    # }else{
    # }
    
  }
}

for (ftrait in traits){
  
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
    # scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    geom_line(aes(group = code_sp),
              alpha=0.4) + #,color=LifeHistory
    ylim(c(miny- 1/10*(maxy-miny),
           maxy)) +
    theme(legend.position="none") +
    # ggtitle ("annuals")
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() ,
          # axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    annotate("text", x = 1, y=maxy, label = comp[1,]$.group)+ 
    annotate("text", x = 2, y=maxy, label = comp[2,]$.group) + 
    
    annotate("text", x = 1, y=miny - 1/10*(maxy-miny), label = nb1)+ 
    annotate("text", x = 2, y=miny - 1/10*(maxy-miny), label = nb2) +
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
    scale_fill_manual(values = c("grey", "white")) +
    # scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    geom_line(aes(group = code_sp),
              alpha = 0.4) + #,color=LifeHistory
    scale_color_manual(values = "#00BFC4") +
    ylim(c(miny - 1/10*(maxy-miny), 
           maxy))+
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
    
    annotate("text", x = 1, y=miny - 1/10*(maxy-miny), label = nb3)+ 
    annotate("text", x = 2, y=miny- 1/10*(maxy-miny), label = nb4)
  
  
  
  # Create a text grob (for the title = trait name)
  tgrob <- ggpubr::text_grob(traits_names[i],size = 10)
  # Draw the text
  plot_0 <- ggpubr::as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
  
  plot <- ggpubr::ggarrange(A, B,
                            ncol = 2,nrow = 1,heights = c(1,5))
  plot2 <- annotate_figure(plot, top = text_grob(traits_names[i], 
                                        color = "black", 
                                        # face = "bold", 
                                        size = 12))
  PLOTS[[i]] <- plot2
  i <- i+1
}

## Boxplot ####

# Extract the legend alone, from the data frame of species removal expe
plot <- MEAN %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>%
  mutate(Management = if_else(treatment=="Fer","Intensive","Extensive")) %>% 
  mutate(Management = factor(Management,levels = c("Intensive","Extensive"))) %>% 
  ggplot(aes_string(x="Management", y=ftrait, label = "code_sp",fill = "Management")) + #,shape = "LifeHistory"
  theme_classic() +
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
                                        ncol = 3,nrow = 3)


boxplot_all_traits_legend <- ggpubr::ggarrange(boxplot_all_traits,
                                               legend,
                                               ncol = 2, nrow = 1,
                                               widths = c(1, 0.2))
# boxplot_all_traits_legend
ggsave("draft/boxplot_all_traits.jpg",boxplot_all_traits_legend,width = 12, height = 8)

## Table ####
table_trait_diff <- TABLE_PVAL %>% 
  mutate(percent_var_code_sp = round(percent_var_code_sp,digits = 2)*100) %>% 
  transmute(
    Trait = Trait, 
    Treatment = ifelse(Treatment =="<0.05",
                       cell_spec(Treatment, bold = T),
                       cell_spec(Treatment, bold=F)),
    Life_History = ifelse(Treatment =="<0.05",
                         cell_spec(Life_History, bold = T),
                         cell_spec(Life_History, bold=F)),
    Interaction = ifelse(Treatment =="<0.05",
                         cell_spec(Interaction, bold = T),
                         cell_spec(Interaction, bold=F)),
    Annuals     =Annuals     ,
    Perennials = Perennials
    # Var_ITV = percent_var_code_sp
  ) %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait", "Regime", "LifeHistory","Interaction",
                                   "Int. - Ext. (annuals)", "Int. - Ext. (perennials)"
                                   # "Var explained by species (%)"
                                   )) %>%
  kableExtra::kable_styling("hover", full_width = F)
  


cat(table_trait_diff, file = "draft/table_mixed_model_boxplot.doc")


