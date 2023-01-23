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
  select(soil,INN,INP)

IN2 <- IN %>% 
  group_by(soil) %>% 
  summarise(meanINN = mean(INN,na.rm = T),meanINP = mean(INP,na.rm=T),
            seINN = sd(INN,na.rm=T)/n(), seINP = sd(INP,na.rm=T)/n()) %>% 
  arrange(-meanINN)


to_barplot <- t(as.matrix(IN2[, c("meanINN", "meanINP")]))
to_barplot_se <- t(as.matrix(IN2[, c("seINN", "seINP")]))
#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


## 3) CSR ####
# CWM2_fer <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_fer_completed.csv")
# CWM2_nat <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_nat_completed.csv")

CWM3_fer <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_fer.csv") # not yet completed CSR
CWM3_nat <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_subset_nat_sab_int.csv")

## 4) Abundance ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))







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
BarPlot <- barplot(to_barplot[,c(1,4)],
                   beside = TRUE, names.arg = c("Intensive", "Extensive"),#names.arg = IN2$soil, # yaxt = "n",
                   ylim=c(0, max(c(IN2$meanINN,IN2$meanINP)+3)  ),  
                   col = c("black","grey"),
                   ylab = "Nutrition index (%)",
                   xaxt = "n"
)
title("A - Nutrient availability",adj = 0,line = 0.5)

legend("topright",
       legend = c("NNI (%)","PNI (%)"),
       pch = 15,
       col = c("black","grey"))
error.bar(BarPlot,to_barplot[,c(1,4)],to_barplot_se[,c(1,4)])




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
        ylab = "Productivity (T/ha)",
        xaxt = "n",
        medlwd = 1
)
title("B - Productivity",adj = 0,line = 0.5)

## Disturbance ####

labels <- c("Intensive",
            "Extensive")

boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  
        width=c(1,4),
        # col=c("orange" , "seagreen"),
        xlab = "",
        ylab = "Proportion of biomass eaten",
        # xaxt = "n",
        at = c(1,2),
        medlwd = 1,
        names = labels
)

title("C - Disturbance",adj = 0,line = 0.5)






##Richness of annual species ####
#specify path to save PDF to
# destination = "outputs/figures/fig_richness_annuals.pdf"
# # destination2 = "outputs/figures/1_abundance_richness_annuals.jpg"
# 
# #open PDF
# pdf(file=destination,width = 4, height = 4)
# png(file=destination2)

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
        ylab = "Relative richness \n of annuals (%)",
        # xaxt = "n",
        medlwd = 1,
        # xaxt = "n",
        names = labels)

title("D - Richness",adj = 0,line = 0.3)

dev.off()






























# legend("topleft", legend="d", bty='n')


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

text(x=c(1:4),y=c(85,40,40,60), labels = comp_binom$.group)


#turn off PDF plotting
dev.off() 





# disturbance %>% 
#   group_by(Trtmt) %>% 
#   summarize(mean = mean(Tx_CalcPic))


#________________________________________________________________________________
# V) Optional graphs ####

## only C score ####
boxplot(
  value ~ score * depth , data = CSR_toplot %>% filter(score == "CWM_C"), xaxt = "n",
  xlab = "", ylab = "Score (%)",
  col = "black",
  medlwd = 1
)
legend(
  "topleft", title = "Score",
  legend = c("C"), fill = c("black"),
  cex = 0.7
)


## CSR Triangle ####
CSR_toplot2 <- CSR_toplot %>%
  spread(key = score,value = value)

library(Ternary)
TernaryPlot(main = "La Fage",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
AddToTernary(points, CSR_toplot2 %>%
               filter(depth == "Fer") %>%
               select(CWM_C,CWM_S,CWM_R), col='green', lty='dotted', lwd=3, pch = 0)
AddToTernary(points, CSR_toplot2 %>%
               filter(depth == "D") %>%
               select(CWM_C,CWM_S,CWM_R), col='black', lty='dotted', lwd=3, pch = 1)
AddToTernary(points, CSR_toplot2 %>%
               filter(depth == "I") %>%
               select(CWM_C,CWM_S,CWM_R), col='black', lty='dotted', lwd=3, pch = 4)
AddToTernary(points, CSR_toplot2 %>%
               filter(depth == "S") %>%
               select(CWM_C,CWM_S,CWM_R), col='black', lty='dotted', lwd=3, pch = 3)
legend("topright",
       cex = 1,
       bty = "n",
       legend = c('Deep', 'Intermediary', 'Superficial'),
       pch = c(1,4,3))

## Annual cover ####
jpeg("draft/cover_nat.jpeg")
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


boxplot(tot_relat_ab * 100 ~ depth,
        data = ann_nat, # cover_annuals, # 
        # ylim=c(0,90),
        xlab = NA,
        ylab = "Relative abundance \n of annuals (%)",
        xaxt = "n",
        medlwd = 1)

ann_nat %>%
  ungroup() %>% 
  filter(depth == "S") %>% 
  summarize(mean_ab = mean(tot_relat_ab))

mod <- lm(tot_relat_ab ~ depth, data =  cover_annuals%>% filter(!(depth=="Fer"))) # ou ann_nat
anova(mod)
summary(mod)
shapiro.test(residuals(mod)) # not normal --> kruskal-wallis ?
lmtest::bptest(mod)
lmtest::dwtest(mod)

posthoc <- multcomp::cld(emmeans::emmeans(mod, specs = "depth",  type = "response",
                                          adjust = "tuckey"),
                         Letters = "abcdefghi", details = T)

comp <- as.data.frame(posthoc$emmeans) %>% 
  arrange(factor(depth, levels = levels(cover_annuals$depth)))

# text(x=c(1:4),y=c(82,30,40,50), labels = comp$.group)
text(x=c(1:3),y=c(20,20,20), labels = comp$.group)

dev.off()

# plutôt faire kruskal test
# En gros : si je bosse sur le natif seul, il y a un effet du gradient de prof de sol.
# Mais cet effet est noyé par la grosse différence avec le fertile, qui 
# embarque la majorité de la variance.

# C'est aussi visible avec un test non paramétrique (kruskal-wallis)
cover_annuals %>% 
  ggplot(aes(x=depth,y=tot_relat_ab))+
  geom_point()
kruskal.test(tot_relat_ab  ~ depth,data = cover_annuals )
FSA::dunnTest(tot_relat_ab ~ depth,data = cover_annuals )

kruskal.test(tot_relat_ab  ~ depth,data = cover_annuals %>% filter(!(depth=="Fer")) )
FSA::dunnTest(tot_relat_ab ~ depth,data = cover_annuals %>% filter(!(depth=="Fer")) )

cover_annuals %>% 
  ggplot(aes(x=tot_relat_ab))+
  geom_histogram(sts="identity") +
  facet_wrap(~depth)
# On n'a pas la même forme de distribution dans les trois cas de figure --> pas approprié d'utiliser Kruskal pour fertile...


cover_annuals %>% 
  group_by(depth) %>% 
  summarize(mean_relat_ab = mean(tot_relat_ab))
         