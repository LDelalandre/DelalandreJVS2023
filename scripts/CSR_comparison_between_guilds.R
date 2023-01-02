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
  
# Life forms ####
MEAN_CSR_shallow %>% 
  filter( (treatment == "Fer" & code_sp %in% ab_fer$code_sp) | (treatment == "Nat" & code_sp %in% ab_nat$code_sp) ) %>% 
  group_by(treatment,LifeForm1) %>% 
  summarize(n=n())


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




# Species replacement ####
ab_ann_nat <- ab_nat %>% filter(LifeHistory == "annual") %>% pull(code_sp) %>% unique()
ab_per_nat <- ab_nat %>% filter(LifeHistory == "perennial") %>% pull(code_sp) %>% unique()
ab_ann_fer <- ab_fer %>% filter(LifeForm1== "The") %>% pull(code_sp) %>% unique()
ab_per_fer <- ab_fer %>% filter(!(LifeForm1== "The")) %>% pull(code_sp) %>% unique()

length(ab_ann_fer)
length(ab_ann_nat)
length(intersect(ab_ann_nat,ab_ann_fer))

length(ab_per_fer)
length(ab_per_nat)
length(intersect(ab_per_nat,ab_per_fer))


platelay <- data.frame(rown = rep(letters[1:8], 4),
                       coln = rep(1:4,each = 8),
                       colorvar = rnorm(32, 0.3, 0.2))
ggplot(platelay,aes(x=rown,y=colorvar,fill=coln)) + 
  geom_bar(position="dodge",stat="identity")

data.frame(LifeHistory = c("Annual","Annual","Perennial","Perennial"), treatment = c("Fer","Nat","Fer","Nat"), nb_sp = , status # dire si elles sont en commun dans les trtmts ) 


ggplot(SUMM3,aes(x=LifeHistory,y=nb_sp,fill=status))+
  scale_color_manual(values = cols, aesthetics = "fill") +
  geom_histogram(stat="identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 66,hjust=1)) +
  ylim(c(0,2))+
  facet_wrap(~site) +
  geom_text(aes(label=count), vjust=0,position = position_dodge(width = 1)) +
  labs(fill = "Category of species") +
  xlab("Property") +
  ylab("Mean effect of key species") +
  ggsave(paste0("figures_tables/mean_effect_key_sp.png"),height=7,width=7)