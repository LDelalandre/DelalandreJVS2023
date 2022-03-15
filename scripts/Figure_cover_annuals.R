source("scripts/1. Packages.R")

# https://r-charts.com/base-r/margins/

# I) Load data ####
# NB: data mainly generated from 2021_11_18_CSR.R

# 1) Vegetation consumption ####
disturbance <- read.table("data/environment/Disturbance_DivHerbe.txt",header=T,sep="\t",dec=",") %>% 
  filter(!(Trtmt == "Tem"))


# 2) INN ####
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


# 3) CSR ####
CWM2_fer <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_fer_completed.csv")
CWM2_nat <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_nat_completed.csv")

CWM3_fer <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_fer.csv") # not yet completed CSR
CWM3_nat <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_nat.csv")

# 4) Abundance ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(!code_sp =="STRIFSCAB") # /!\ à corriger à la source !!

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))



# II) Figure ####
# NB: Version pré-finale, à modifier sous inkscape
#specify path to save PDF to
destination = "outputs/figures/Fig_envir_cover/main figure.pdf"

#open PDF
pdf(file=destination,width = 3, height = 7)

#specify to save plots in 2x2 grid
par(mar = c(2, 5.5, 0.1, 1),mfrow = c(6,1),mpg = c(1,1,0))

# 0) Number of annual species ####
# NB : must be done at the level of plot, I guess.
richness_per_guild_nat <- ab_nat %>% 
  count(depth,paddock,LifeHistory) %>% 
  merge(soil_Maud,by=c("depth","paddock"))

richness_per_guild_fer <- ab_fer %>%
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  mutate(depth = "Fer") %>% 
  count(depth,paddock,LifeHistory) %>% 
  mutate(PC1score = NA)

richness_per_guild <- rbind(richness_per_guild_nat,richness_per_guild_fer)
richness_per_guild$depth <- factor(richness_per_guild$depth , levels = c("Fer","D","I","S"))

richness_per_guild_toplot <- richness_per_guild %>% 
  spread(LifeHistory,n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(relative_richness_annual = annual / (annual + perennial))

boxplot(relative_richness_annual *100 ~ depth,
        data = richness_per_guild_toplot,
        # ylim=c(0,100),
        xlab = NA,
        ylab = "Relative richness \n of annuals (%)",
        xaxt = "n",
        medlwd = 1)

# 1) Annual cover ####
ann_fer <- ab_fer %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  group_by(LifeHistory,year,paddock,id_transect_quadrat) %>% 
  dplyr::rename(line = id_transect_quadrat) %>% 
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

# (Test modèles) ####
# abondance absolue
ab_nat %>% 
  filter(LifeHistory=="annual") %>% 
  group_by(depth,paddock,PC1score,line) %>% 
  summarize(abundance = mean(abundance)) %>% 
  ggplot(aes(x=PC1score,y=abundance,label=paddock))+
    geom_point()
  # geom_smooth(method = "lm")

#abondance relative
ann_nat_all <- ab_nat %>% 
  group_by(LifeHistory,depth,paddock,line) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  # filter(LifeHistory =="annual") %>% 
  mutate(year = 2009) %>% 
  relocate(LifeHistory,year,depth,paddock,line,tot_relat_ab) %>% 
  mutate(line = as.character(line))

ann_nat_all %>% 
  group_by(depth,paddock,line) %>% 
  summarize(sum = sum(tot_relat_ab))

ann_nat_soil <- merge(ann_nat,soil_Maud,by = c("depth","paddock"))
ggplot(ann_nat_soil,aes(x=PC1score, y=tot_relat_ab))+
  geom_point() 
  # geom_smooth(method="lm")
mod <- lm(data = ann_nat_soil, tot_relat_ab ~ PC1score)
plot(mod)
anova(mod)
summary(mod)

ggplot(ann_nat_soil,aes(x=depth, y=tot_relat_ab))+
  geom_boxplot() 
mod <- lm(data = ann_nat_soil, tot_relat_ab ~ depth)
# plot(mod)
anova(mod)
summary(mod)

# graphe avec les pérennes aussi
ann_nat_soil_all <- merge(ann_nat_all,soil_Maud,by = c("depth","paddock"))
ggplot(ann_nat_soil_all,aes(x=depth, y=tot_relat_ab,color= LifeHistory))+
  geom_boxplot() 

# richesse absolue
# NB refaire : j'ai mois de points pour la richesse (par plot) 
# que pour l'abondance (par transect) --> plus de puissance!
# evolution richesse absolue comm gt
richness_per_guild_nat %>% 
  filter(LifeHistory=="annual") %>% 
  ggplot(aes(x=PC1score,y=n))+
  geom_point() 
  # geom_smooth(method="lm")

#annuelles
richness_per_guild_nat %>% 
  filter(LifeHistory=="annual") %>% 
  ggplot(aes(x=depth,y=n))+
  geom_boxplot() 
mod <- lm(data = richness_per_guild_nat %>% filter(LifeHistory=="annual"), 
          n ~ depth)

# Dans Fer pour comparer richesse absolue
richness_per_guild_fer %>% 
  ggplot(aes(x=depth,y=n,color=LifeHistory))+
  geom_boxplot() +
  ylim(c(0,85))
richness_per_guild_nat %>% 
  ggplot(aes(x=depth,y=n,color=LifeHistory))+
  geom_boxplot() +
  ylim(c(0,85))


mod <- lm(data = richness_per_guild_nat %>% filter(LifeHistory=="annual"), 
          n ~ PC1score)
# plot(mod)
anova(mod)
summary(mod)
# NB refaire par transect et pas plot pour la puissance stat.

# " pérennes"  
richness_per_guild_nat %>% 
  filter(LifeHistory=="perennial") %>% 
  ggplot(aes(x=PC1score,y=n))+
  geom_point()+
  geom_smooth(method="lm")

# toutes
richness_per_guild_nat %>% 
  group_by(PC1score,depth) %>%
  summarize(n=sum(n)) %>%
  ggplot(aes(x=PC1score,y=n))+
  geom_point()+
  geom_smooth(method="lm")

# plot annuals rich relat
ggplot(richness_per_guild_toplot %>% filter(!(depth == "Fer")),
       aes(x=PC1score,y=relative_richness_annual)) +
  geom_point()
mod <- lm(data = richness_per_guild_toplot %>% filter(!(depth == "Fer")), 
          relative_richness_annual ~ PC1score)
# plot(mod)
anova(mod)
summary(mod)

ggplot(richness_per_guild_toplot %>% filter(!(depth == "Fer")),
       aes(x=depth,y=relative_richness_annual)) +
  geom_boxplot()
mod <- lm(data = richness_per_guild_toplot %>% filter(!(depth == "Fer")), 
          relative_richness_annual ~ depth)
# plot(mod)
anova(mod)
summary(mod)
# fin ds tests 

cover_annuals <- rbind(ann_fer,ann_nat)
cover_annuals$depth <- factor(cover_annuals$depth , levels = c("Fer","D","I","S"))


boxplot(tot_relat_ab * 100 ~ depth,
        data = cover_annuals, # %>% filter(!(depth == "Fer")),
        ylim=c(0,10),
        xlab = NA,
        ylab = "Relative abundance \n of annuals (%)",
        xaxt = "n",
        medlwd = 1)
# Est-ce qu'il faut que je bosse au niveau de la ligne, ou du parc? Au niveau du parc, j'ai deux points pour le fertile,
# sauf si je prends plusieurs années.

# 2) CWM CSR scores ####
CSR_toplot_nat <- CWM2_nat %>% 
  mutate(id_com = paste(depth,paddock,line,sep = "_")) %>% 
  select(id_com,depth,CWM_C,CWM_S,CWM_R) %>% 
  gather(key = "score", value = "value", -c(depth,id_com) ) 

CSR_toplot_fer <- CWM2_fer %>% 
  mutate(depth = "Fer") %>% 
  select(depth,id_transect_quadrat,CWM_C,CWM_S,CWM_R) %>% 
  gather(key = "score", value = "value", -c(depth,id_transect_quadrat) ) %>% 
  rename(id_com = id_transect_quadrat)

CSR_toplot <- rbind(CSR_toplot_nat,CSR_toplot_fer)
CSR_toplot$score <- factor(CSR_toplot$score , levels = c("CWM_C","CWM_S","CWM_R"))
CSR_toplot$depth <- factor(CSR_toplot$depth , levels = c("Fer","D","I","S"))


boxplot(
  value ~ score * depth , data = CSR_toplot, 
  xaxt = "n",
  xlab = "", 
  ylab = "Score (%)",
  col = c("black", "grey","white"),
  medlwd = 1
)
legend(
  "topleft", title = "Score",
  legend = c("C", "S", "R"), fill = c("black", "grey","white"),
  cex = 0.7
)

# only C score
# boxplot(
#   value ~ score * depth , data = CSR_toplot %>% filter(score == "CWM_C"), xaxt = "n",
#   xlab = "", ylab = "Score (%)",
#   col = "black",
#   medlwd = 1
# )
# legend(
#   "topleft", title = "Score",
#   legend = c("C"), fill = c("black"),
#   cex = 0.7
# )


# Triangle
# CSR_toplot2 <- CSR_toplot %>% 
#   spread(key = score,value = value)
# 
# TernaryPlot(main = "La Fage",alab = "C \u2192", blab = "S \u2192", clab = "\u2190 R")
# AddToTernary(points, CSR_toplot2 %>% 
#                filter(depth == "Fer") %>% 
#                select(CWM_C,CWM_S,CWM_R), col='green', lty='dotted', lwd=3, pch = 0)
# AddToTernary(points, CSR_toplot2 %>% 
#                filter(depth == "D") %>% 
#                select(CWM_C,CWM_S,CWM_R), col='black', lty='dotted', lwd=3, pch = 1)
# AddToTernary(points, CSR_toplot2 %>% 
#                filter(depth == "I") %>% 
#                select(CWM_C,CWM_S,CWM_R), col='black', lty='dotted', lwd=3, pch = 4)
# AddToTernary(points, CSR_toplot2 %>% 
#                filter(depth == "S") %>% 
#                select(CWM_C,CWM_S,CWM_R), col='black', lty='dotted', lwd=3, pch = 3)
# legend("topright", 
#        cex = 1,
#        bty = "n",
#        legend = c('Deep', 'Intermediary', 'Superficial'), 
#        pch = c(1,4,3))


# Link score C in fertile and abundance of annuals
C_ab_ann_fer <- merge(CSR_toplot_fer %>% rename(line = id_com),ann_fer,by=c("depth","line"))
ggplot(C_ab_ann_fer %>% filter(score=="CWM_C"),aes(x=value,y=tot_relat_ab))+
  geom_point()
# Il faudrait faire un vrai modèle stat, avec plusieurs variables explicatives...

C_rich_ann_fer <- merge(CSR_toplot_fer %>% rename(line = id_com),ann_fer,by=c("depth","line"))
ggplot(C_ab_ann_fer %>% filter(score=="CWM_C"),aes(x=value,y=tot_relat_ab))+
  geom_point()


# 2bis) Ellenberg ####
ellenberg_toplot_nat <- CWM3_nat %>% 
  mutate(id_com = paste(depth,paddock,line,sep = "_")) %>% 
  select(id_com,depth,CWM_light,CWM_moisture,CWM_nitrogen,CWM_pH) %>% 
  gather(key = "score", value = "value", -c(depth,id_com) ) 

ellenberg_toplot_fer <- CWM3_fer %>% 
  mutate(depth = "Fer") %>% 
  select(id_transect_quadrat,depth,CWM_light,CWM_moisture,CWM_nitrogen,CWM_pH) %>% 
  gather(key = "score", value = "value", -c(depth,id_transect_quadrat) ) %>% 
  rename(id_com = id_transect_quadrat)

ellenberg_toplot <- rbind(ellenberg_toplot_nat,ellenberg_toplot_fer)
ellenberg_toplot$score <- factor(ellenberg_toplot$score , levels = c("CWM_light","CWM_nitrogen","CWM_moisture"))
ellenberg_toplot$depth <- factor(ellenberg_toplot$depth , levels = c("Fer","D","I","S"))


boxplot(
  value ~ score * depth , data = ellenberg_toplot, 
  # xaxt = "n",
  xlab = "", 
  ylab = "Score (%)",
  col = c("black", "grey","white"),
  medlwd = 1
)
legend(
  "topleft", title = "Score",
  legend = c("light", "nitrogen", "moisture"), fill = c("black", "grey","white"),
  cex = 0.7
)

# 3) Vegetation cover ####
tot_ab_fer <- ab_fer %>% 
  group_by(year,paddock,id_transect_quadrat) %>% 
  dplyr::rename(line = id_transect_quadrat) %>% 
  summarise(tot_ab = sum(abundance)) %>% 
  mutate(depth = "Fer") %>% 
  relocate(year,depth,paddock,line,tot_ab)

tot_ab_nat <- ab_nat %>% 
  group_by(depth,paddock,line) %>% 
  summarise(tot_ab = sum(abundance)) %>% 
  mutate(year = 2009) %>% 
  relocate(year,depth,paddock,line,tot_ab) %>% 
  mutate(line = as.character(line))

cover <- rbind(tot_ab_fer,tot_ab_nat)
cover$depth <- factor(cover$depth , levels = c("Fer","D","I","S"))

boxplot(tot_ab ~depth, data = cover ,
        xlab = NA,
        ylab = "Number of point \n intersept",
        medlwd = 1,
        xaxt = "n")

# 4) Biomass consumption ####
boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  
        width=c(1,4),
        # col=c("orange" , "seagreen"),
        xlab = NA,
        ylab = "Proportion of \n biomass eaten",
        xaxt = "n",
        # log = "x",
        at = c(1,2),
        medlwd = 1
)

# 5) Biomass production ####
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

boxplot(biomass_may$rdt.T.ha ~ biomass_may$Position,
        xlab = NA,
        ylab = "Productivity (T/ha)",
        xaxt = "n",
        medlwd = 1
)

# 6) Bar plot of INN and INP ####
BarPlot <- barplot(to_barplot,
                   beside = TRUE, names.arg = c("Fer", "D","I","S"),#names.arg = IN2$soil, # yaxt = "n",
                   ylim=c(0, max(c(IN2$meanINN,IN2$meanINP)+3)  ),  
                   col = c("black","grey"),
                   ylab = "Nutrition index (%)"
                   )
legend("topright",
       legend = c("INN (%)","INP (%)"),
       pch = 15,
       col = c("black","grey"))
error.bar(BarPlot,to_barplot,to_barplot_se)

#turn off PDF plotting
dev.off() 


# III) Stats ####

# 0) Annual richness ####
lm_richness <- lm(relative_richness_annual ~ depth, data = richness_per_guild_toplot )
anova(lm_richness)
summary(lm_richness)
multcomp::cld(emmeans::emmeans(lm_richness, specs = "depth", type = "response",
                               adjust = "tuckey"),
              Letters = "abcdefghi", details = T)



# 1) Annual cover ####
cover_annuals_nat <- cover_annuals %>% 
  filter(!(depth == "Fer"))

lm_cover <- lm(tot_relat_ab ~ depth, data = cover_annuals_nat, 
               contrasts = list(depth = "contr.SAS"))
plot(lm_cover)

shapiro.test(residuals(lm_cover))
lmtest::bptest(lm_cover)
lmtest::dwtest(lm_cover)

anova(lm_cover)
summary(lm_cover)

multcomp::cld(emmeans::emmeans(lm_cover, specs = "depth", type = "response",
                               adjust = "tuckey"),
              Letters = "abcdefghi", details = T)
