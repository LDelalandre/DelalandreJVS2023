library(tidyverse)
library(vegan)



# Jaccard distance ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  mutate(code_sp = case_when(species == "Vicia sativa ssp. sativa"~ "VICISATI-SAT",
                             species == "Crepis vesicaria ssp. haenseleri" ~"CREPVESI-HAE",
                             species == "Taraxacum laevigatum" ~ "TARALAEV",
                             species == "Cirsium acaulon" ~ "CIRSACAU",
                             species == "Carthamus mitissimus" ~ "CARTMITI",
                             TRUE~ code_sp)) %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial"))
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial"))



# build presence matrix ####
pres_fer <- ab_fer %>% 
  mutate(treatment = "Fer") %>% 
  rename(line = id_transect_quadrat) %>%  
  mutate(presence = if_else(relat_ab >0,1,0)) %>% 
  select(code_sp,LifeHistory,treatment,line,presence)

pres_nat <- ab_nat %>% 
  mutate(treatment = "Nat") %>% 
  filter(depth == "S") %>% 
  mutate(presence = if_else(relat_ab >0,1,0)) %>%
  select(code_sp,LifeHistory,treatment,Ligne,presence) %>% 
  rename(line = Ligne)

presence <- rbind(pres_fer,pres_nat)

presence_matrix <- presence %>% 
  unique() %>% 
  spread(key = code_sp,value = presence) %>% 
  replace(is.na(.),0)

presence_matrix %>% 
  filter(treatment=="Fer") %>%
  pull(line) %>% 
  length()

presence_matrix %>% 
  filter(treatment=="Nat") %>% 
  pull(line) %>% 
  unique() %>%
  length()


# Species turnover between pairs of transects ####
lh <- "perennial"



JAC <- NULL
LH <- NULL
for (lh in c("annual","perennial")){
  fpres_nat <- presence_matrix %>%
    filter(treatment=="Nat") %>% 
    filter(LifeHistory == lh)
  fpres_fer <- presence_matrix %>%
    filter(treatment=="Fer") %>% 
    filter(LifeHistory == lh)
  species <- presence %>% 
    filter(LifeHistory == lh) %>% 
    pull(code_sp) %>% 
    unique()
  for (i in c(1:dim(fpres_nat[1]))){
    for(j in c(1:dim(fpres_fer[1]))){
      A <- fpres_nat[i,] %>% select(-c(LifeHistory,treatment,line))
      B <- fpres_fer[j,] %>% select(-c(LifeHistory,treatment,line))
      jac  <- vegan::vegdist(x = rbind(A,B) %>% 
                               select(all_of(species)),method="jaccard") 
      JAC <- c(JAC,jac)
      LH <- c(LH,lh)
    }
  }
}

jaccard <- data.frame(JAC,LH)
plot_jaccard <- jaccard %>% 
  mutate(LH = if_else(LH == "perennial","Perennials","Annuals")) %>% 
  ggplot(aes(x=LH,y=JAC))+
  geom_boxplot() +
  ylab("Jaccard index between pairs of transects") +
  theme_classic() +
  ggsignif::geom_signif(comparisons = list(c("Annuals", "Perennials")), 
              map_signif_level=TRUE) +
  xlab('')

plot_jaccard

ggsave("draft/plot_jaccard.png",plot = plot_jaccard)

# mod <- lm(JAC ~ LH,data = jaccard)
# # plot(mod)
# anova(mod)
# summary(mod)

wilcox.test(JAC ~ LH,data = jaccard)
