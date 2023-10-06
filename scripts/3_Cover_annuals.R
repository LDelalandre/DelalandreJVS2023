source("scripts/Packages.R")

# https://r-charts.com/base-r/margins/

# I) Load data ####
# NB: data mainly generated from 2021_11_18_CSR.R

## 1) Vegetation consumption ####
disturbance <- read.table("data/environment/Disturbance_DivHerbe.txt",header=T,sep="\t",dec=",") %>% 
  filter(!(Trtmt == "Tem"))


## 2) INN ####
INraw <- read.xlsx("data/environment/IN La Fage.xlsx", sheet = "resultat", startRow = 1, colNames = TRUE) %>% 
  mutate()
nom <- read.xlsx("data/environment/IN La Fage.xlsx", sheet = "nom", startRow = 1, colNames = TRUE)  %>% 
  filter(!is.na(plot))
IN <- merge(nom,INraw,by="nature.echantillon") %>% 
  filter(!(treatment =="temoin")) %>% 
  mutate(soil = if_else(treatment=="fertilise","Fer",plot)) %>% 
  select(soil,INN,INP) %>% 
  filter(soil %in% c("Fer","Sable"))

# II) Figure_envt ####

#specify path to save PDF to
destination = "draft/fig_envt.pdf"
# destination2 = "outputs/figures/1_abundance_richness_annuals.jpg"

#open PDF
pdf(file=destination,width = 6.5, height = 6)
# png(file=destination2)

#specify to save plots in 2x2 grid
par(mar = c(2, 6.5, 1.3, 1),mfrow = c(2,2),mpg = c(1,1,0))


## Bar plot of INN and INP ####
IN_gathered <- IN %>% 
  gather(key = index, value = value, -soil) %>% 
  na.omit() %>% 
  mutate(soil_index = paste(soil,index,sep="_")) 

IN_gathered$soil_index <- factor(IN_gathered$soil_index, levels = c("Fer_INN","Fer_INP","Sable_INN","Sable_INP"))

boxplot(IN_gathered$value ~ IN_gathered$soil_index,
        xlab = NA,
        ylab = "Nutrition index (%)",
        xaxt = "n",
        medlwd = 1,
        col = gray(c(0.3,0.9)))
legend("topright",
       legend = c("NNI (%)","PNI (%)"),
       pch = 15,
       col = gray(c(0.3,0.9)))

title("A - Nutrient availability",adj = 0,line = 0.5)

## Biomass production ####
biomass <- read.xlsx("data/environment/Biomasses et indices La Fage.xlsx", 
                     sheet = "2009", 
                     startRow = 1, colNames = TRUE, rowNames = F) %>% 
  mutate(Dates = as.Date(Dates- 25569, origin = "1970-01-01"))

biomass_may <- biomass %>% 
  filter(Parcs %in% c("C1","C2","1","6","8","10")) %>% 
  filter(Dates == "2009-05-01") %>% 
  mutate(rdt.T.ha = as.numeric(rdt.T.ha)) %>% 
  mutate(Position = str_sub(Position,1,1)) %>% 
  mutate(Position = case_when(is.na(Position) ~ "Fer",
                              TRUE ~ Position)) 
biomass_may$Position <- factor(biomass_may$Position , levels = c("Fer","D","I","S"))
biomass_may_reduced <- biomass_may %>% 
  filter(Position %in% c("Fer","S"))
biomass_may_reduced$Position <- factor(biomass_may_reduced$Position , levels = c("Fer","S"))

boxplot(biomass_may_reduced$rdt.T.ha ~ biomass_may_reduced$Position,
        xlab = NA,
        ylab = "Maximum standing biomass (t/ha)",
        xaxt = "n",
        medlwd = 1
)
title("B - Biomass",adj = 0,line = 0.5)

## Disturbance ####

labels <- c("Intensive",
            "Extensive")

boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  
        width=c(1,1),
        # col=c("orange" , "seagreen"),
        xlab = "",
        ylab = "Proportion of biomass eaten",
        # xaxt = "n",
        at = c(1,2),
        medlwd = 1,
        names = labels
)

title("C - Disturbance",adj = 0,line = 0.5)

# ggplot(aes(x = Trtmt, y = Tx_CalcPic), data = disturbance)+ 
#   geom_boxplot()+
#   geom_point()




##Richness of annual species ####
#specify path to save PDF to
# destination = "outputs/figures/fig_richness_annuals.pdf"
# # destination2 = "outputs/figures/1_abundance_richness_annuals.jpg"
# 
# #open PDF
# pdf(file=destination,width = 4, height = 4)
# png(file=destination2)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))

richness_per_guild_nat <- ab_nat %>% 
  count(depth,paddock,LifeHistory,line) %>% 
  merge(soil_Maud,by=c("depth","paddock"))

richness_per_guild_fer <- ab_fer %>%
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  mutate(depth = "Fer") %>% 
  count(depth,paddock,LifeHistory,line) %>% 
  mutate(PC1score = NA)

richness_per_guild <- rbind(richness_per_guild_nat,richness_per_guild_fer)
richness_per_guild$depth <- factor(richness_per_guild$depth , levels = c("Fer","D","I","S"))

richness_per_guild_toplot <- richness_per_guild %>% 
  spread(LifeHistory,n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(relative_richness_annual = annual / (annual + perennial)) %>% 
  mutate(zone = case_when(depth == "Fer"~ "G+F",
                          depth == "D" ~ "GU-D",
                          depth == "I" ~ "GU-I",
                          TRUE ~ "GU-S"))
richness_per_guild_toplot_reduced <- richness_per_guild_toplot %>% 
  filter(zone %in% c("G+F","GU-S"))
# relative richness

boxplot(relative_richness_annual *100 ~ zone,
        data = richness_per_guild_toplot_reduced,
        ylim=c(0,90),
        xlab = NA,
        ylab = "Relative richness of annuals (%)",
        # xaxt = "n",
        medlwd = 1,
        # xaxt = "n",
        names = labels)

title("D - Richness",adj = 0,line = 0.3)


dev.off()




## Abundance of annual species ####
png(file="draft/boxplot abundance.png")

ann_fer <- ab_fer %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  group_by(LifeHistory,year,paddock,line) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  filter(LifeHistory =="annual") %>% 
  mutate(depth = "Fer") %>% 
  relocate(LifeHistory,year,depth,paddock,line,tot_relat_ab)

ann_nat <- ab_nat %>% 
  group_by(LifeHistory,depth,paddock,line) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  filter(LifeHistory =="annual") %>%
  mutate(year = 2009) %>% 
  relocate(LifeHistory,year,depth,paddock,line,tot_relat_ab) %>% 
  mutate(line = as.character(line))

cover_annuals <- rbind(ann_fer,ann_nat)
cover_annuals$depth <- factor(cover_annuals$depth , levels = c("Fer","D","I","S"))
cover_annuals2 <- cover_annuals %>% 
  filter(depth %in% c("Fer","S"))
cover_annuals2$depth <- factor(cover_annuals2$depth , levels = c("Fer","S"))


boxplot(tot_relat_ab * 100 ~ regime,
        data = cover_annuals2 %>% 
          mutate(regime = if_else(depth == "Fer","Intensive","Extensive")) %>% 
          mutate(regime = factor(regime,levels = c("Intensive","Extensive"))), # cover_annuals, # 
        # ylim=c(0,90),
        xlab = NA,
        ylab = "Relative abundance of annuals (%)",
        # xaxt = "n",
        medlwd = 1)
# title("Abundance",adj = 0,line = 0.3)

dev.off()




# peut-être faire test non paramétrique (kruskal-wallis)

# glm " mais ça marche pas, je suis sur de la richesse relative, donc proportion, pas poisson !!
# --> A refaire !
mod_binom0 <- glm(cbind(annual,perennial) ~ 1, family = "binomial", data=richness_per_guild_toplot)
mod_binom <- glm(cbind(annual,perennial) ~ depth, family = "binomial", data=richness_per_guild_toplot)
# ou family =  binomial(logit)
anova(mod_binom)
summary(mod_binom)
anovaII <- car::Anova(mod_binom)

# donc je regarde la RS absolue des annuelles, "annual", et pas relat, "relative_richness_annual"
# mod_glm <- glm(relative_richness_annual ~ depth, family = "quasipoisson", data = richness_per_guild_toplot)
# quasipoisson to deal with overdispersion.
emm <- emmeans::emmeans(mod_binom, 
                        specs = pairwise ~ "depth",  
                        type = "response",
                        adjust = "tukey")
emm$emmeans
emm$contrasts
# A sortir en appendix

comp_binom <- as.data.frame(emm$emmeans) %>% 
  arrange(factor(depth, levels = levels(richness_per_guild_toplot$depth))) %>% 
  mutate(.group = c("a","b","b","c")) # /!\ recode manually significant levels

# text(x=c(1:4),y=c(85,40,40,60), labels = comp_binom$.group)


#turn off PDF plotting
dev.off() 



# Figures pour ppt ####
## Bar plot of INN and INP ####

barplot(to_barplot[1,c(1,4)],
        col=c("#F8766D","#00BFC4"),
        cex.axis=3,
        las = 2 ,
        ylim = c(0, 80)) # texte axe y horizontal)

barplot(to_barplot[2,c(1,4)],
        col=c("#F8766D","#00BFC4"),
        cex.axis=3,
        las = 2 ,
        ylim = c(0, 80))

## Prod ####
boxplot(biomass_may_reduced$rdt.T.ha ~ biomass_may_reduced$Position,
        xlab = NA,
        ylab = NA,
        xaxt = "n",
        medlwd = 1,
        # border=c("#F8766D","#00BFC4"),
        col=c("#F8766D","#00BFC4"),
        # boxlwd = 8,
        # whisklwd = 5,
        # staplelwd = 8,
        cex.axis=3,
        las = 2 # texte axe y horizontal
)

## Disturb ####
boxplot(disturbance$Tx_CalcPic * 100 ~ disturbance$Trtmt ,  
        width=c(1,1),
        # col=c("orange" , "seagreen"),
        xlab = NA,
        ylab = NA,
        xaxt = "n",
        at = c(1,2),
        medlwd = 1,
        names = labels,
        col=c("#F8766D","#00BFC4"),
        las = 2 ,
        cex.axis=2.5
)

## Richness ####
boxplot(relative_richness_annual *100 ~ zone,
        data = richness_per_guild_toplot_reduced,
        ylim=c(0,90),
        xlab = NA,
        ylab = NA,
        # xaxt = "n",
        medlwd = 1,
        xaxt = "n",
        names = labels,
        col=c("#F8766D","#00BFC4"),
        las = 2 ,
        cex.axis=2.5)

