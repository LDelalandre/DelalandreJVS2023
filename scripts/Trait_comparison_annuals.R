library("tidyverse")
source("scripts/Data_traits.R")

# Model traits of annuals as a function of their origin
trait_sheet <- c("LeafMorpho","LeafCN","Leaf13C","Biovolume","Pheno","Seed")

TRAITS <- list(c("LDMC","SLA","L_Area"),
            c("LCC","LNC"),# "LPC",
            c("Ldelta13C"),
            c("Hveg"  ,    "Hrepro"   , "Dmax"  , "Dmin") ,
            c("Disp"),#,"Mat_Per", #"Mat","Flo",
            c("SeedMass")
)
trait_unit <- read.csv2("data/trait_unit.csv")

BOXPLOT <- list() # boxplot of trait comparisons between annuals in fer and in nat sup
P.VAL <- c()
INTERCEPT <- c()
ESTIMATE <- c()
count <- 0 # position in the BOXPLOT list
for (i in 1:length(trait_sheet)){
  fsheet <- trait_sheet[i] # name of the sheet we focus on
  fdata <- get(fsheet) %>% # content of that sheet
    filter(Treatment %in% c("Fer_Clc","Fer_Dlm","Nat_Sab")) %>% 
    mutate(treatment = str_sub(Treatment,1,3)) %>% 
    filter(LifeForm1=="The")  
  
  in_nat <- fdata %>% 
    filter(treatment =="Nat") %>% 
    pull(Code_Sp) %>% 
    unique() 
  in_fer <- fdata %>% 
    filter(treatment =="Fer") %>% 
    pull(Code_Sp) %>% 
    unique()
  sp_both <- intersect(in_nat,in_fer)
  
  fdata2 <- fdata %>% filter(Code_Sp %in% sp_both)
  
      
  ftraits <- TRAITS[[i]] # traits in this sheet
  for ( j in 1:length(ftraits) ){
    trait <- ftraits[j] # choose one trait in the sheet
    
    # plot
    plot <- fdata %>% 
      # filter(measurementDeterminedBy=="Leo Delalandre") %>%
      ggplot(aes_string(x="treatment",y= trait ))+
      geom_point()
    count <- count+1
    BOXPLOT[[count]] <- plot

    # stats
    formula <- as.formula(paste0(trait, " ~ treatment", " + (1|Code_Sp)"))
    formula0 <- as.formula(paste0(trait, " ~ 1 + (1|Code_Sp)"))
    
    if ( length( fdata %>% pull(treatment) %>% unique() ) == 2 ){ # if we have data from nat and fer
      mmod <- lme4::lmer( formula , data = fdata) # /!\ choose fdata (includes sp just measured in on treatment)
      # or fdata2 (sp measured in both treatments only)
      # mmod0 <- lme4::lmer( formula0 , data = fdata)
      # anova <- anova(mmod0,mmod)
      anova <- car::Anova(mmod)
      summary <- summary(mmod)
      
      pval <- anova$`Pr(>Chisq)`
      P.VAL <- c(P.VAL,pval)
      
      estimate <- summary$coefficients[2,1]
      ESTIMATE <- c(ESTIMATE,estimate)
      
      intercept <- summary$coefficients[1,1]
      INTERCEPT <- c(INTERCEPT,intercept)
    }else{
      P.VAL <- c(P.VAL,NA)
      ESTIMATE <- c(ESTIMATE,NA)
      INTERCEPT <- c(INTERCEPT,NA)
    }

    
    # qqnorm(residuals(mmod),main="")
    # plot(fitted(mmod),residuals(mmod),xlab="Fitted",ylab="Residuals") ; abline(h=0)
    
  }
}
BOXPLOT

vec_traits <- unlist(TRAITS, recursive = TRUE, use.names = TRUE)

df_mod <- data.frame(Trait = vec_traits,
                     Intercept = INTERCEPT,
                    Estimate = ESTIMATE,
                    p.value = P.VAL) %>% 
  mutate(Estimate = round(Estimate,digits = 2)) %>% 
  mutate(Estimate = Estimate + Intercept) %>% 
  mutate(p.value = format(p.value, scientific=TRUE, digits=3)) %>% 
  mutate_if(is.numeric, round, digits = 2) 


# export table
table_mod <- df_mod %>%
  merge(trait_unit,by="Trait") %>% 
  select(Trait,Unit,Intercept,Estimate,p.value) %>% 
  arrange(factor(Trait,levels = vec_traits)) %>% 
  kableExtra::kable( escape = F,
         col.names = c("Trait",
                       "Unit",
                       "Value in G+F",
                       "Value in GUs",
                       "p.value")) %>%
  kableExtra::kable_styling("hover", full_width = F)

table_mod

cat(table_mod, file = "draft/mixed_model_intra_annual.doc")

