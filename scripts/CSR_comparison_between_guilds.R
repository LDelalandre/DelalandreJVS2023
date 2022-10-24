source("scripts/Packages.R")

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  
  # merge(name_LH,by="Code_Sp") %>% 
  # relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())

MEAN_CSR_shallow <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_subset_nat_sab_completed.csv",dec=",",fileEncoding="latin1") %>% 
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
boxplot_CSR_shallow <- MEAN_CSR_shallow %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  # filter(score == "R") %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  ggplot(aes(x=score, y=value, color = LifeHistory))+
  theme_classic()+
  ylab("Score (%)")+
  theme(axis.title.x=element_blank())+
  geom_boxplot() +
  # geom_point() +
  facet_wrap(~zone,labeller=variable_labeller) 
  

boxplot_CSR_shallow

ggsave("draft/boxplot_CSR.jpg",boxplot_CSR_shallow)

data.anovaCSR <- MEAN_CSR_shallow %>% 
  select(code_sp,treatment,LifeHistory,C,S,R) %>%
  # gather(key = score, value = value, -c(code_sp, LifeHistory,treatment)) %>% 
  filter(treatment%in% c("Nat","Fer"))

options(contrasts=c("contr.treatment","contr.treatment"))
lmCSR <- lm(R ~ treatment * LifeHistory, data = data.anovaCSR )

par(mfrow=c(2,2)) ; plot(lmCSR) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(lmCSR))) # normality_graph
shapiro.test(residuals(lmCSR)) # normality of residuals
lmtest::bptest(lmCSR) # homoscedasticity
lmtest::dwtest(lmCSR)  # non-independence of residuals!
anova(lmCSR)
summary(lmCSR)


lmCSR_all <- lm(R ~ treatment * LifeHistory, data = data.anovaCSR )
lmCSR_no_interaction <- lm(R ~ treatment + LifeHistory, data = data.anovaCSR )
anova(lmCSR_all,lmCSR_no_interaction)


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