source("scripts/Packages.R")

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

# MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  
  ## merge(name_LH,by="Code_Sp") %>%
  ## relocate(C,S,R) %>%
  # mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  # mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  # mutate(R=str_replace(R,",",".") %>% as.numeric())

MEAN_CSR_shallow <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_completed_seed_mass.csv",
                              dec=",",
                              encoding="latin1") %>% 
  # merge(name_LH,by="Code_Sp") %>% 
  # relocate(C,S,R) %>% 
  # mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  # mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  # mutate(R=str_replace(R,",",".") %>% as.numeric()) %>% 
  filter(!(species == "Geranium dissectum - pétiole"))
  


# CSR ####

variable_names <- list(
  "G+F" = expression(paste("G"^'+',"F",sep='')),
  "GU-S" = expression(paste("GU"[S],sep='')) )

variable_labeller <- function(variable,value){
  return(variable_names[value])
}

# Natif: only shallow (sandy) soil
# boxplot_CSR_shallow <- MEAN_CSR_shallow %>% 
#   select(code_sp,treatment,LifeHistory,C,S,R) %>%
#   gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
#   # filter(score == "C") %>%
#   filter(treatment%in% c("Nat","Fer")) %>% 
#   mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
#   mutate(score = factor(score,levels = c("C","S","R"))) %>% 
#   ggplot(aes(x=score, y=value, color = LifeHistory))+
#   theme_classic()+
#   ylab("Score (%)")+
#   theme(axis.title.x=element_blank())+
#   geom_boxplot() +
#   # geom_point(position = position_dodge(width = .75)) +
#   facet_wrap(~zone,labeller=variable_labeller) 
#   # geom_line(aes(group = score,code_sp),position = position_dodge(width = .75))
  
boxplot_CSR_shallow <- MEAN_CSR_shallow %>%
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  # filter(score == "C") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  mutate(score = factor(score,levels = c("C","S","R"))) %>% 
  mutate(grouping = interaction(LifeHistory,code_sp)) %>% 
  
  ggplot(aes(x=LifeHistory, y=value, fill = zone))+
  theme_classic()+
  ylab("Score (%)")+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  # geom_point(position = position_dodge(width = .75), aes(color = LifeHistory)) +
  facet_wrap(~score)  +
  scale_fill_manual(values = c("grey", "white")) 
  geom_line(aes(group= grouping))

boxplot_CSR_shallow

ggsave("draft/boxplot_CSR.jpg",boxplot_CSR_shallow)

data.anovaCSR <- MEAN_CSR_shallow %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(treatment%in% c("Nat","Fer"))

options(contrasts=c("contr.treatment","contr.treatment"))
lmCSR <- lm(S ~ treatment * LifeHistory, data = data.anovaCSR )

par(mfrow=c(2,2)) ; plot(lmCSR) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(lmCSR))) # normality_graph
shapiro.test(residuals(lmCSR)) # normality of residuals
lmtest::bptest(lmCSR) # homoscedasticity
lmtest::dwtest(lmCSR)  # non-independence of residuals!
anova(lmCSR)
summary(lmCSR)

# # specifically test interaction
# lmCSR_all <- lm(R ~ treatment * LifeHistory, data = data.anovaCSR )
# lmCSR_no_interaction <- lm(R ~ treatment + LifeHistory, data = data.anovaCSR )
# anova(lmCSR_all,lmCSR_no_interaction)
# summary(lmCSR_all)

posthoc <- multcomp::cld(emmeans::emmeans(lmCSR, specs = c("treatment","LifeHistory"),  type = "response",
                                          adjust = "tuckey"),
                         Letters = "abcdefghi", details = T)

comp <- as.data.frame(posthoc$emmeans)
comp

# CSR between annuals ####
MEAN_CSR_shallow %>% 
  filter(LifeHistory == "annual") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  # filter(score == "R") %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  ggplot(aes(x=score, y=value, color = zone))+
  theme_classic()+
  ylab("Score (%)")+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  ggtitle("annuals")



MEAN_CSR_shallow %>% 
  filter(LifeHistory == "annual") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(score == "C") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  spread(key = score, value = value) %>%
  ggplot(aes(x=zone, y=C, color = zone)) +
  theme_classic()+
  ylab("Score (%)") +
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp)) +
  ggtitle("C score ")


MEAN_CSR_shallow %>% 
  filter(LifeHistory == "annual") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(score == "S") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  spread(key = score, value = value) %>%
  ggplot(aes(x=zone, y=S, color = zone,label = code_sp)) +
  theme_classic()+
  ylab("Score (%)") +
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp)) +
  ggtitle("S score ") 
# ggrepel::geom_text_repel()
# CERAGLOM, GERADISS et SHERARVE réagissent à l'envers en score S

MEAN_CSR_shallow %>% 
  filter(LifeHistory == "annual") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(score == "R") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  spread(key = score, value = value) %>%
  ggplot(aes(x=zone, y=R, color = zone)) +
  theme_classic()+
  ylab("Score (%)") +
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp)) +
  ggtitle("R score ")



# CSR between perennials ####

MEAN_CSR_shallow %>% 
  filter(LifeHistory == "perennial") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  # filter(score == "R") %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  ggplot(aes(x=score, y=value, color = zone))+
  theme_classic()+
  ylab("Score (%)")+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  ggtitle("perennials")



MEAN_CSR_shallow %>% 
  filter(LifeHistory == "perennial") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(score == "C") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  spread(key = score, value = value) %>%
  ggplot(aes(x=zone, y=C, color = zone)) +
  theme_classic()+
  ylab("Score (%)") +
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp)) +
  ggtitle("C score ")

MEAN_CSR_shallow %>% 
  filter(LifeHistory == "perennial") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(score == "S") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  spread(key = score, value = value) %>%
  ggplot(aes(x=zone, y=S, color = zone,label = code_sp)) +
  theme_classic()+
  ylab("Score (%)") +
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp)) +
  ggtitle("S score ") 



MEAN_CSR_shallow %>% 
  filter(LifeHistory == "perennial") %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(score == "R") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  spread(key = score, value = value) %>%
  ggplot(aes(x=zone, y=R, color = zone)) +
  theme_classic()+
  ylab("Score (%)") +
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp)) +
  ggtitle("R score ")

# compute the % of sp in common, in both guilds





# On any trait ####
TRAITS <- list(c("LDMC","SLA","L_Area"),
               c("LCC","LNC"),# "LPC",
               c("Ldelta13C"),
               c("Hveg"  ,    "Hrepro"   , "Dmax"  , "Dmin") ,
               c("Disp"),#,"Mat_Per", #"Mat","Flo",
               c("SeedMass")
)

trait <- "L_Area"

# /!\ log pour leaf area
MEAN_CSR_shallow %>% 
  filter(LifeHistory == "annual") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  # ggplot(aes(x=zone, y=log(L_Area), label = code_sp)) +
  ggplot(aes_string(x="zone", y=trait, label = "code_sp")) +
  theme_classic()+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = code_sp))  
  # ggrepel::geom_text_repel()
  



  
  

## model on everyone, as for CSR ####

options(contrasts=c("contr.treatment","contr.treatment"))
mod <- lm(LDMC ~ treatment * LifeHistory, data = MEAN_CSR_shallow )

par(mfrow=c(2,2)) ; plot(mod) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(mod))) # normality_graph
shapiro.test(residuals(mod)) # normality of residuals
lmtest::bptest(mod) # homoscedasticity
lmtest::dwtest(mod)  # non-independence of residuals!
anova(mod)
summary(mod)

# # specifically test interaction
# lmCSR_all <- lm(R ~ treatment * LifeHistory, data = data.anovaCSR )
# lmCSR_no_interaction <- lm(R ~ treatment + LifeHistory, data = data.anovaCSR )
# anova(lmCSR_all,lmCSR_no_interaction)
# summary(lmCSR_all)

posthoc <- multcomp::cld(emmeans::emmeans(mod, specs = c("treatment","LifeHistory"),  type = "response",
                                          adjust = "tuckey"),
                         Letters = "abcdefghi", details = T)

comp <- as.data.frame(posthoc$emmeans)
comp



## model on annuals, but on traits averaged at sp level ####
trait <- "LDMC"

formula <- as.formula(paste0(trait, " ~ treatment * LifeHistory", " + (1|code_sp)"))
formula0 <- as.formula(paste0(trait, " ~ 1 + (1|code_sp)"))

mmod <- lme4::lmer( formula , data = MEAN_CSR_shallow) # /!\ choose fdata (includes sp just measured in on treatment)
# or fdata2 (sp measured in both treatments only)
# mmod0 <- lme4::lmer( formula0 , data = fdata)
# anova <- anova(mmod0,mmod)
car::Anova(mmod)
summary(mmod)

anova$`Pr(>Chisq)`

summary$coefficients[2,1]




## TESTS PLOTS ####
MEAN_CSR_shallow %>% 
  # filter(LifeHistory == "annual") %>%
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  # ggplot(aes(x=zone, y=log(L_Area), label = code_sp)) +
  ggplot(aes_string(x="LifeHistory", y=trait, label = "code_sp",fill = "zone")) +
  theme_classic()+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  geom_point(aes(color = LifeHistory),position = position_dodge(width = .75)) +
  scale_fill_manual(values = c("grey", "white"))
# extraire la légende !

# Si seulement dans abondance
data_fer <- MEAN_CSR_shallow %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN_CSR_shallow %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_CSR_shallow_in_abundance <- rbind(data_fer,data_nat)


TRAITS <- c("LDMC","SLA","L_Area",
"LCC","LNC","Ldelta13C",
"Hveg"  ,    "Hrepro"   , "Dmax"  , "Dmin" ,
"Disp",
"SeedMass",
"C","S","R"
)
# TRAITS <- c("C","S","R")


PLOTS <- NULL
i <- 1
for (trait in TRAITS){
  
  miny <- min(MEAN_CSR_shallow %>% pull(sym(trait)))
  maxy <- max(MEAN_CSR_shallow %>% pull(sym(trait)))
  
  A <- MEAN_CSR_shallow_in_abundance %>% 
    filter(LifeHistory == "annual") %>%
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
    
    ggplot(aes_string(x="zone", y=trait, label = "code_sp",fill = "zone")) +
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
    ggtitle(trait)
  
  B <- MEAN_CSR_shallow_in_abundance %>% 
    filter(LifeHistory == "perennial") %>%
    filter(treatment%in% c("Nat","Fer")) %>% 
    mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
    
    ggplot(aes_string(x="zone", y=trait, label = "code_sp",fill = "zone")) +
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
    ggeasy::easy_remove_y_axis()+
    # ggtitle("perennials") 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() ,
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) 
  
  plot <- ggpubr::ggarrange(A,B)
  PLOTS[[i]] <- plot
  i <- i+1
}


# Extract the legend alone, from the data frame of species removal expe
plot <- MEAN_CSR_shallow_in_abundance %>% 
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
# Simplifier la légende 
# Juste points de couleur pour LifeHistory
# juste boxplots pour zone
# Faire des symboles pour les noms de zones
# Changer le nom LifeHistory



boxplot_all_traits <- ggpubr::ggarrange(PLOTS[[1]],PLOTS[[2]],PLOTS[[3]],PLOTS[[4]],PLOTS[[5]],PLOTS[[6]],
                  PLOTS[[7]],PLOTS[[8]],PLOTS[[9]],PLOTS[[10]],PLOTS[[11]],PLOTS[[12]],
                  PLOTS[[13]],PLOTS[[14]],PLOTS[[15]],legend)

boxplot_CSR_traits <- ggpubr::ggarrange(PLOTS[[13]],PLOTS[[14]],PLOTS[[15]],legend)

ggsave("draft/boxplot_all_traits.jpg",boxplot_all_traits)
ggsave("draft/boxplot_CSR_traits.jpg",boxplot_CSR_traits)
