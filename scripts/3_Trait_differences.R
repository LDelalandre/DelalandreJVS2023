library(tidyverse)
library(ggpubr)
library(kableExtra)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv") %>%
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  mutate(log_LA = log(L_Area)) %>% 
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
            "SeedMass"
)#"H_FLORE","FLO_FLORE",

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

TABLE %>% 
  # filter(Trait %in% c("LDMC","LCC","Ldelta13C","Hrepro","Disp","SeedMass")) %>% #"H_FLORE","FLO_FLORE",
  kableExtra::kable( escape = F,
                     col.names = c("Origin", "Annuals", "Perennials")) %>%
  kableExtra::kable_styling("hover", full_width = F)

# AJOUTER JACCARD GLOBAL (= SUR RELEVES BOTANIQUES)

#_______________________________________________________________________________
# Interspecific comparisons ####


# changer les noms des traits
# passer LA en log
# changer les noms des zones

TABLE_PVAL <- NULL
PLOTS <- NULL
i <- 1
for (ftrait in traits){
  ## Linear model ####
  data.anovaCSR <- MEAN %>% 
    # select(code_sp,treatment,LifeHistory,C,S,R) %>%
    # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(log_L_Area = log(L_Area))
  
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

    # if( table_pval$Interaction < 0.05 ){
      posthoc <- multcomp::cld(emmeans::emmeans(mmod, specs = c("treatment","LifeHistory"),  type = "response",
                                                adjust = "tuckey"),
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
        mutate(est_pval = paste0(estimate, "(",p.value,")")) %>%
        select(contrast,est_pval) %>% 
        spread(key = contrast,value = est_pval)
      
      table_pval2 <- cbind(table_pval,contrasts)
      TABLE_PVAL <- rbind(TABLE_PVAL,table_pval2)
      
      comp <- as.data.frame(posthoc$emmeans) %>% 
        mutate(treatment = factor(treatment,levels = c("Fer","Nat"))) %>% 
        mutate(LifeHistory = factor(LifeHistory, levels = c("annual","perennial"))) %>% 
        arrange(LifeHistory,treatment)
    # }else{
    # }
    
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

## Boxplot ####

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
                                        PLOTS[[7]],PLOTS[[8]],PLOTS[[9]],legend)

ggsave("draft/boxplot_all_traits.jpg",boxplot_all_traits,width = 14, height = 8)

## Table ####
table_trait_diff <- TABLE_PVAL %>% 
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
  ) %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait", "Treatment", "LifeHistory","Interaction",
                                   "Fer-Nat (annuals)", "Fer-Nat (perennials)")) %>%
  kableExtra::kable_styling("hover", full_width = F)
  


cat(table_trait_diff, file = "draft/differences_traits.doc")



# Intraspecific comparisons ####


diff_to_random <- function(intrasp_var){
  # create columns randfer and randnat, which randomize whether 
  # the trait value was observed in fer or nat of the given species 
  random_diff <- intrasp_var %>% 
    mutate(rand1 = rbinom(length(Fer),1,0.5),
           rand2 = 1 - rand1) %>% 
    mutate(randfer = if_else(rand1 == 1, Fer,Nat),
           randnat = if_else(rand2==1, Fer,Nat)) %>% 
    mutate(randdiff = randnat - randfer) %>% 
    group_by(LifeHistory) %>% 
    summarize(mean_diff = mean(randdiff)) %>% 
    arrange(LifeHistory)
  # vector with random differences for annuals and perennials
  random_diff$mean_diff
}


diff_ann_per <- function(intrasp_var){
  random_diff <- intrasp_var %>% 
    select(code_sp,LifeHistory,diff) %>% 
    mutate(LifeHistory_random = sample(LifeHistory)) %>% 
    group_by(LifeHistory_random) %>% 
    summarise(mean_diff = mean(diff)) %>% 
    arrange(LifeHistory_random )
  random_diff$mean_diff[2] %>% abs() -  random_diff$mean_diff[1] %>% abs()
}

## test1: existence of intrasp var ####

PVAL_ANN <- NULL
PVAL_PER <- NULL
for (ftrait in traits){
  intrasp_var <- MEAN %>% 
    # filter(LifeHistory=="annual") %>% 
    select(species,code_sp,LifeHistory,treatment,ftrait) %>% 
    spread(key=treatment,value=ftrait) %>% 
    mutate(trait = ftrait) %>% 
    na.omit() %>%  # keep only sp measured in the 2 treatments
    mutate(ratio = Nat / Fer) %>% 
    mutate(diff = Nat-Fer)
  
  # ggplot(intrasp_var,aes(x=LifeHistory,y=diff)) +
  #   geom_boxplot() +
  #   geom_point() +
  #   ggtitle(ftrait)
  
  real_difference <- intrasp_var %>%
    group_by(LifeHistory) %>% 
    summarize(mean_diff = mean(diff)) %>% 
    arrange(LifeHistory) %>% 
    pull(mean_diff)
  
  TEST_DIFF_ZERO <- data.frame(trait = ftrait,
                               randomization = "real",
                               mean_diff_ann = real_difference[1],
                               mean_diff_per = real_difference[2])
  nb_boot <- 50
  for (i in c(1:nb_boot)){
    randomized_difference <- diff_to_random(intrasp_var)
    test_diff_zero <- data.frame(trait = ftrait,
                                 randomization = i,
                                 mean_diff_ann = randomized_difference[1],
                                 mean_diff_per = randomized_difference[2])
    TEST_DIFF_ZERO <- rbind(TEST_DIFF_ZERO,test_diff_zero)
  }
  
  # plot(density(TEST_DIFF_ZERO %>% 
  #                filter(!(randomization == "real")) %>% 
  #                pull(mean_diff_ann)
  #              ))
  # abline(v = TEST_DIFF_ZERO[1,3])
  # 
  # plot(density(TEST_DIFF_ZERO %>% 
  #                filter(!(randomization == "real")) %>% 
  #                pull(mean_diff_per)
  # ))
  # abline(v = TEST_DIFF_ZERO[1,4])
  
  pval_ann <- min(
    length(which(TEST_DIFF_ZERO$mean_diff_ann < TEST_DIFF_ZERO[1,3] )) / dim(TEST_DIFF_ZERO)[1] ,
    length(which(TEST_DIFF_ZERO$mean_diff_ann > TEST_DIFF_ZERO[1,3] )) / dim(TEST_DIFF_ZERO)[1]
  )
  pval_per <- min(
    length(which(TEST_DIFF_ZERO$mean_diff_per < TEST_DIFF_ZERO[1,4] )) / dim(TEST_DIFF_ZERO)[1],
    length(which(TEST_DIFF_ZERO$mean_diff_per > TEST_DIFF_ZERO[1,4] )) / dim(TEST_DIFF_ZERO)[1]
  )
  PVAL_ANN <- c(PVAL_ANN,pval_ann)
  PVAL_PER <- c(PVAL_PER,pval_per)
  
}


test_plasticity <- data.frame(trait = traits,
           pval_ann = PVAL_ANN,
           pval_per = PVAL_PER)

  
## test2: difference in intrasp var between annuals and perennials ####
PVAL <- NULL
for (ftrait in traits){
  intrasp_var <- MEAN %>% 
    # filter(LifeHistory=="annual") %>% 
    select(species,code_sp,LifeHistory,treatment,ftrait) %>% 
    spread(key=treatment,value=ftrait) %>% 
    mutate(trait = ftrait) %>% 
    na.omit() %>%  # keep only sp measured in the 2 treatments
    mutate(ratio = Nat / Fer) %>% 
    mutate(diff = Nat-Fer)
  
  DIFF_LH <- NULL
  for (j in c(1:nb_boot)){
    diff_LH <- diff_ann_per(intrasp_var)
    DIFF_LH <- c(DIFF_LH,diff_LH)
  }
  real_diff <- intrasp_var %>% 
    group_by(LifeHistory) %>% 
    summarize(mean_diff = mean(diff)) %>%
    arrange(LifeHistory) %>% 
    pull(mean_diff)
  real_abs_diff <- abs(real_diff[2]) - abs(real_diff[1])
  
  # plot(density(DIFF_LH))
  # abline(v = real_abs_diff)
  pval <- length(which(DIFF_LH < real_abs_diff)) / length(DIFF_LH)
  
  PVAL <- c(PVAL,pval)
  
}


## table ####
difference_in_plasticity <- data.frame(trait = traits,
           p.value_difference = PVAL)

table_intrasp <- merge(test_plasticity,difference_in_plasticity) %>% 
  arrange(match(trait,traits)) %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait", "Annuals", "Perennials", "Difference")) %>%
  kableExtra::kable_styling("hover", full_width = F) 
  

cat(table_intrasp, file = "draft/differences_variation_intrasp.doc")
