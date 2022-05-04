# Charger les données auparavant

# Regarder dans diachro nombre points contacts f(traitement)
ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  # filter(year == 2004) %>%
  group_by(id_transect_quadrat,year) %>% 
  mutate(sumab = sum(abundance)) %>% 
  ggplot(aes(x=treatment,y=sumab))+
  geom_boxplot() +
  facet_wrap(~year) 
# ggsave("outputs/plots/comparison_intercept_nat_fer.png",width = 10,height=10)

# Regarder la météo ####
meteo <- read.csv2("data/environment/Meteo_LaFage_1973-2006.csv")
meteo %>% 
  group_by(AN) %>% 
  summarize(sum_RR = sum(RR), mean_TX = mean(TX), mean_RGC = mean(RGC)) %>% 
  ggplot(aes(x=AN,y=mean_TX))+
  geom_line()
