library(tidyverse)
library(ggpubr)
library(kableExtra)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv") %>%
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  mutate(log_LA = log(L_Area)) %>% 
  unique() %>% 
  dplyr::rename(LCCm = LCC) %>% 
  dplyr::rename(LNCm = LNC)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

# I case I want to perform the analyses with species both in the trait and abundance data:
data_fer <- MEAN %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_intersect <- rbind(data_fer,data_nat)


traits <- c("LDMC","SLA","log_LA",
            "LCCm","LNCm","Ldelta13C",#"LPC",
            "Hrepro" ,  #, "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)#"H_FLORE","FLO_FLORE",

traits_names <- c("Leaf Dry Matter Content (mg/g)", "Specific Leaf Area (cm²/kg)"," log(Leaf Area (cm²))",
                  "Leaf Carbon Content (mg/g)","Leaf Nitrogen Content (mg/g)", "Leaf delta 13C (part per thousand)",
                  "Reproductive Height (cm)", 
                  "Date of first dispersal (Julian day)",
                  "Seed Mass (mg)")

# Chose to take species in trait data, or in both trait and abundance data


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

diff_to_random_all <- function(intrasp_var){
  # create columns randfer and randnat, which randomize whether 
  # the trait value was observed in fer or nat of the given species 
  random_diff <- intrasp_var %>% 
    mutate(rand1 = rbinom(length(Fer),1,0.5),
           rand2 = 1 - rand1) %>% 
    mutate(randfer = if_else(rand1 == 1, Fer,Nat),
           randnat = if_else(rand2==1, Fer,Nat)) %>% 
    mutate(randdiff = randnat - randfer) %>% 
    summarize(mean_diff = mean(randdiff))
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
# randomize treatment in which species were measured
# measure difference between trait values in the two treatments
# compute proportion of times when real difference was either greater or smaller than randomized
nb_boot <- 1000

PVAL_ANN <- NULL
PVAL_PER <- NULL
PVAL_ALL <- NULL
for (ftrait in traits){
  # compute trait difference and ratio across the two treatments
  intrasp_var <- MEAN %>% 
    # filter(LifeHistory=="annual") %>% 
    select(species,code_sp,LifeHistory,treatment,ftrait) %>% 
    spread(key=treatment,value=ftrait) %>% 
    mutate(trait = ftrait) %>% 
    na.omit() %>%  # keep only sp measured in the 2 treatments
    mutate(ratio = Nat / Fer) %>% 
    mutate(diff = Nat-Fer) %>% 
    mutate(RDPI = diff/Fer)
  # Relative Distances Plasticity Index (RDPI)
  
  # ggplot(intrasp_var,aes(x=LifeHistory,y=diff)) +
  #   geom_boxplot() +
  #   geom_point() +
  #   ggtitle(ftrait)
  
  # compute difference 
  # real_difference <- intrasp_var %>%
  #   group_by(LifeHistory) %>% 
  #   summarize(mean_diff = mean(RDPI)) %>% 
  #   arrange(LifeHistory) %>% 
  #   pull(mean_diff)
  real_difference_all <- intrasp_var %>%
    summarize(mean_diff = mean(diff)) %>% 
    pull(mean_diff)
  
  
  TEST_DIFF_ZERO <- data.frame(trait = ftrait,
                               randomization = "real",
                               # mean_diff_ann = real_difference[1],
                               # mean_diff_per = real_difference[2],
                               mean_diff_all = real_difference_all)
  
  for (i in c(1:nb_boot)){
    # randomized_difference <- diff_to_random(intrasp_var)
    randomized_difference_all <- diff_to_random_all(intrasp_var)
    test_diff_zero <- data.frame(trait = ftrait,
                                 randomization = i,
                                 # mean_diff_ann = randomized_difference[1],
                                 # mean_diff_per = randomized_difference[2],
                                 mean_diff_all = randomized_difference_all)
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
  
  # pval_ann <- min(
  #   length(which(TEST_DIFF_ZERO$mean_diff_ann <= TEST_DIFF_ZERO[1,3] )) / dim(TEST_DIFF_ZERO)[1] ,
  #   length(which(TEST_DIFF_ZERO$mean_diff_ann >= TEST_DIFF_ZERO[1,3] )) / dim(TEST_DIFF_ZERO)[1]
  # )
  # 
  # pval_per <- min(
  #   length(which(TEST_DIFF_ZERO$mean_diff_per <= TEST_DIFF_ZERO[1,4] )) / dim(TEST_DIFF_ZERO)[1],
  #   length(which(TEST_DIFF_ZERO$mean_diff_per >= TEST_DIFF_ZERO[1,4] )) / dim(TEST_DIFF_ZERO)[1]
  # )
  
  pval_all <- min(
    length(which(TEST_DIFF_ZERO$mean_diff_all <= TEST_DIFF_ZERO[1,3] )) / dim(TEST_DIFF_ZERO)[1],
    length(which(TEST_DIFF_ZERO$mean_diff_all >= TEST_DIFF_ZERO[1,3] )) / dim(TEST_DIFF_ZERO)[1]
  )
  
  # PVAL_ANN <- c(PVAL_ANN,pval_ann)
  # PVAL_PER <- c(PVAL_PER,pval_per)
  PVAL_ALL <- c(PVAL_ALL,pval_all)
  
}

write.csv2(data.frame(trait = traits,
                      pval_all = PVAL_ALL), 
           "outputs/data/test_plasticity_1000_bootstrap.csv",row.names=F)




## test2: difference in intrasp var between annuals and perennials ####

# wilcox.test
PVAL <- NULL
STAT <- NULL
for (ftrait in traits){
  intrasp_var <- MEAN %>% 
    # filter(LifeHistory=="annual") %>% 
    select(species,code_sp,LifeHistory,treatment,ftrait) %>% 
    spread(key=treatment,value=ftrait) %>% 
    mutate(trait = ftrait) %>% 
    na.omit() %>%  # keep only sp measured in the 2 treatments
    mutate(ratio = Nat / Fer) %>% 
    mutate(diff = Nat-Fer)
  # https://perso.ens-lyon.fr/lise.vaudor/test-de-wilcoxon-mann-whitney/
  test <- wilcox.test(diff ~ LifeHistory, data = intrasp_var)
  
  PVAL <- c(PVAL, test$p.value )
  STAT <- c(STAT,test$statistic)
}
difference_in_plasticity <- data.frame(trait = traits,
                                       W = STAT,
                                       p.value = PVAL %>% round(digits = 2))



# randomizations (avoid it: adds to much complexity)
# PVAL <- NULL
# for (ftrait in traits){
#   intrasp_var <- MEAN %>% 
#     # filter(LifeHistory=="annual") %>% 
#     select(species,code_sp,LifeHistory,treatment,ftrait) %>% 
#     spread(key=treatment,value=ftrait) %>% 
#     mutate(trait = ftrait) %>% 
#     na.omit() %>%  # keep only sp measured in the 2 treatments
#     mutate(ratio = Nat / Fer) %>% 
#     mutate(diff = Nat-Fer)
#   
#   DIFF_LH <- NULL
#   for (j in c(1:nb_boot)){
#     diff_LH <- diff_ann_per(intrasp_var)
#     DIFF_LH <- c(DIFF_LH,diff_LH)
#   }
#   real_diff <- intrasp_var %>% 
#     group_by(LifeHistory) %>% 
#     summarize(mean_diff = mean(diff)) %>%
#     arrange(LifeHistory) %>% 
#     pull(mean_diff)
#   real_abs_diff <- abs(real_diff[2]) - abs(real_diff[1])
#   
#   # plot(density(DIFF_LH))
#   # abline(v = real_abs_diff)
#   pval <- length(which(DIFF_LH < real_abs_diff)) / length(DIFF_LH)
#   
#   PVAL <- c(PVAL,pval)
#   
# }

# write.csv2 ("outputs/data/comparison_plastcity.csv")

## table ####
test_plasticity <- read.csv2("outputs/data/test_plasticity_1000_bootstrap.csv")
# ATTENTION NE PAS TESTER PLASTICITE SUR H_repro, Disp, SeedMass = j'ai complété les données


table_intrasp <- merge(test_plasticity %>% select(trait,pval_all),difference_in_plasticity) %>% 
  arrange(match(trait,traits)) %>% 
  mutate(pval_all = round(pval_all,digits =3)) %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait", "p.value", "W","p.value")) %>%
  kableExtra::kable_styling("hover", full_width = F)  %>% 
  kableExtra::add_header_above(c(" "=1,
                                 "Intraspecific variation" = 1, 
                                 "Effect of life history" = 2))

cat(table_intrasp, file = "draft/table_intrasp_change.doc")
