source("scripts/Packages.R")

# https://r-charts.com/base-r/margins/

# I) Load data ####
# NB: data mainly generated from 2021_11_18_CSR.R

## Environmental data ####
env_data <- read.csv("outputs/data/env_data.csv")

## Richness and abundance ####

data_abundance <- read.csv("outputs/data/data_abundance.csv")
ab_fer <- data_abundance %>% filter(treatment == "Int")
ab_nat <- data_abundance %>% filter(treatment == "Ext") 



# II) Figure_envt ####

#specify path to save PDF to
destination = "draft/fig_envt.pdf"
# destination2 = "outputs/figures/1_abundance_richness_annuals.jpg"

#open PDF
pdf(file=destination,width = 6.5, height = 6)
# png(file=destination2)

#specify to save plots in 2x2 grid
par(mar = c(2, 6.5, 1.3, 1),mfrow = c(2,2))#,mpg = c(1,1,0)


## Boxplot of INN and INP ####
IN_toplot <- env_data  %>% 
  filter(variable %in% c("NNI", "PNI")) %>% 
  mutate(soil_index = paste(management,variable,sep="_")) 

IN_toplot$soil_index <- factor(IN_toplot$soil_index, levels = c("Int_NNI","Int_PNI","Ext_NNI","Ext_PNI"))

boxplot(IN_toplot$value ~ IN_toplot$soil_index,
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

biomass_toplot <- env_data  %>% 
  filter(variable == "Biomass")
biomass_toplot$management = factor(biomass_toplot$management, levels = c("Int", "Ext"))

boxplot(biomass_toplot$value ~ biomass_toplot$management,
        xlab = NA,
        ylab = "Maximum standing biomass (t/ha)",
        xaxt = "n",
        medlwd = 1
)
title("B - Biomass",adj = 0,line = 0.5)

## Disturbance ####

disturbance_toplot <- env_data  %>% 
  filter(variable == "Disturbance")
disturbance_toplot$management = factor(disturbance_toplot$management, levels = c("Int", "Ext"))


labels <- c("Intensive",
            "Extensive")

boxplot(disturbance_toplot$value ~ disturbance_toplot$management ,  
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


soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))

richness_per_guild_nat <- ab_nat %>% 
  count(depth,paddock,LifeHistory,line) 

richness_per_guild_fer <- ab_fer %>%
  count(depth,paddock,LifeHistory,line) 

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


# dev.off()




## Abundance of annual species ####
png(file="draft/boxplot abundance.png")

ann_fer <- ab_fer %>% 
  group_by(LifeHistory,paddock,line,depth,year) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  filter(LifeHistory =="annual") %>% 
  relocate(LifeHistory,depth,paddock,line,tot_relat_ab) %>% 
  mutate(year)

ann_nat <- ab_nat %>% 
  group_by(LifeHistory,depth,paddock,line,year) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  filter(LifeHistory =="annual") %>%
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

ppt <- F

if (ppt == T){
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
  
  
}
