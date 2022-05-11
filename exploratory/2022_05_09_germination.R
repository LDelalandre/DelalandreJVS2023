library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv")%>%
  filter(!is.na(SLA)) 

# Germination ####
germination <- read.csv2("data/traits/germination.csv") %>% 
  mutate(code_sp = toupper(code_sp)) %>% 
  filter(!(code_sp == "CAPSBURS"))

germination %>% 
  mutate(percent_germ = germ_nb/Nseed_per_pot) %>% 
  ggplot(aes(x=code_sp,y=percent_germ,color = treatment)) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_boxplot() +
  scale_color_manual(values = c("blue","red"))


germination_line <- germination %>% 
  group_by(code_sp, treatment, line_plaque, table, column, line) %>% 
  summarize(nb_seed = sum(Nseed_per_pot),
            nb_seed_germ = sum(germ_nb)) %>% 
  mutate(percent_germ = nb_seed_germ / nb_seed) %>% 
  ungroup() %>% 
  select(code_sp,treatment,percent_germ) %>% 
  filter(!is.na(percent_germ))

germination_line %>% 
  ggplot(aes(x=code_sp,y=percent_germ, color=treatment ))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c("blue","red"))

germination_mean <- germination_line %>% 
  group_by(code_sp,treatment) %>% 
  summarize(mean = mean(percent_germ))

# Différence de taux de germination en fonction du Rs, ou de l'indice d'Ellenberg ?
MEAN %>% 
  merge(germination_mean,by=c("code_sp","treatment")) %>% 
  ggplot(aes(x= LNC ,y=mean,color = treatment)) +
  geom_point()


