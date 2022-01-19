source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
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


# 3) CSR and abundance ####

# 3.1) Maud####
# relative abundance of annuals (maud data: pin-point) ####
ab_maud_traits <- read.csv2("outputs/data/ab_Maud_traits.csv")
rel_ab_annuals <- ab_maud_traits %>% 
  group_by(depth,paddock,PC1score) %>% 
  select(PC1score,depth,paddock,species,abundance,LifeHistory) %>% 
  mutate(rel_ab = abundance/sum(abundance)) %>% 
  filter(LifeHistory == "annual") %>% 
  summarize(sum_ab_ann = sum(rel_ab)) %>% 
  ungroup() %>% 
  select(-c("depth"))

# Compute CWM of CSR: 1) CWM of traits ; 2) Pierce's algo for CSR
CWM_Maud_goodunit <- read.csv2("outputs/data/CWM_Maud.CSV") %>% 
  select(PC1score,depth,CWM_L_Area,CWM_LDMC,CWM_SLA) %>% 
  rename(SLA=CWM_SLA,LDMC=CWM_LDMC,L_Area=CWM_L_Area) %>% 
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(CWM_Maud_goodunit,"outputs/data/Pierce CSR/Traits_CWM_Maud.csv" ,row.names=F)
CSR_Maud <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_Maud_completed.csv") 

CSR_PC1_cover <- left_join(CSR_Maud,rel_ab_annuals,by = "PC1score") %>% 
  replace(is.na(.), 0)
# Rq: cover of annuals is 14% in parc 1, which is the most windy of La Fage

# en boxplot
ggplot(CSR_PC1_cover,aes(x=depth,y=sum_ab_ann))+
  geom_boxplot()

ggplot(CSR_PC1_cover,aes(x=depth,y=S))+
  geom_boxplot()
ggplot(CSR_PC1_cover,aes(x=depth,y=C))+
  geom_boxplot()
ggplot(CSR_PC1_cover,aes(x=depth,y=R))+
  geom_boxplot()

# 3.2) Adeline ####
adeline_ab_CSR_toplot_ann_fer <- read.csv2("outputs/data/Adeline_abundance_CSR.csv") %>% 
  filter(LifeHistory == "annual" & treatment == "Fer") %>% 
  rename(sum_ab_ann = sum_ab_relat)


# 3.3) Merge the two ####
CSR_ab_ann <- rbind(
  adeline_ab_CSR_toplot_ann_fer %>% 
  select(C,S,R,sum_ab_ann) %>% 
  mutate(depth = "Fer"),
  
  CSR_PC1_cover %>% 
  select(C,S,R,sum_ab_ann,depth)
)
CSR_ab_ann$depth <- factor(CSR_ab_ann$depth, levels = c("Fer","D","I","S"))







# II) Figure ####
# NB: Version pré-finale, à modifier sous inkscape
#specify path to save PDF to
destination = "outputs/figures/Fig_envir_cover/test3.pdf"

#open PDF
pdf(file=destination,width = 3, height = 7)

#specify to save plots in 2x2 grid
par(mar = c(2, 4.5, 0.1, 1),mfrow = c(4,1),mpg = c(1,1,0))

# 1) Boxplot of annual cover ####
boxplot(CSR_ab_ann$sum_ab_ann * 100 ~ CSR_ab_ann$depth,
        ylim=c(0,100),
        xlab = NA,
        ylab = "Relative abundance of annuals (%)",
        xaxt = "n",
        medlwd = 1)

# 2) Boxplots of CWM CSR scores ####

# boxplot(CSR_ab_ann$R ~ CSR_ab_ann$depth,        
#         xlab = NA,
#         ylab = "Community-weighted R score",
#         xaxt = "n")
# boxplot(CSR_ab_ann$S ~ CSR_ab_ann$depth,
#         xlab = NA,
#         ylab = "Community-weighted S score",
#         xaxt = "n")
# boxplot(CSR_ab_ann$C ~ CSR_ab_ann$depth,
#         xlab = NA,
#         ylab = "Community-weighted C score",
#         xaxt = "n")

CSR_ab_ann_gathered <- CSR_ab_ann %>% 
  gather(key = "score", value = "value", -c(depth,sum_ab_ann) )
CSR_ab_ann_gathered$score <- factor(CSR_ab_ann_gathered$score , levels = c("C","S","R"))

boxplot(
  value ~ score * depth , data = CSR_ab_ann_gathered, xaxt = "n",
  xlab = "", ylab = "Score (%)",
  col = c("black", "grey","white"),
  medlwd = 1
)

legend(
  "topleft", title = "Score",
  legend = c("C", "S", "R"), fill = c("black", "grey","white"),
  cex = 0.7
)

# 
# SandR <- CSR_ab_ann_gathered %>% 
#   filter(score %in% c("S","R"))
# SandR$score <- factor(SandR$score , levels = c("S","R"))
# 
# boxplot(
#   value ~ score * depth , data = SandR, 
#   xaxt = "n",
#   xlab = "",
#   ylab = "Score (%)",
#   col = c("white", "grey"),
#   medlwd = 1
# )
# 
# legend(
#   "topleft", title = "Score",
#   legend = c("S", "R"), fill = c("white", "grey"),
#   cex = 0.7
# )


# 3) Boxplot of biomass consumption ####
boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  
        width=c(1,4), 
        # col=c("orange" , "seagreen"),
        xlab = NA,
        ylab = "Proportion of biomass eaten",
        xaxt = "n",
        # log = "x",
        at = c(1,2),
        medlwd = 1
)

# 4) Bar plot of INN and INP ####
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

# 1) Annual cover ####
CSR_ab_ann_nat <- CSR_ab_ann %>% 
  filter(!(depth == "Fer"))

lm_cover <- lm(sum_ab_ann ~ depth, data = CSR_ab_ann_nat )
plot(lm_cover)

shapiro.test(residuals(lm_cover))
lmtest::bptest(lm_cover)
lmtest::dwtest(lm_cover)

anova(lm_cover)
summary(lm_cover)

multcomp::cld(emmeans::emmeans(lm_cover, specs = "depth", type = "response",
                               adjust = "tuckey"),
              Letters = "abcdefghi", details = T)
