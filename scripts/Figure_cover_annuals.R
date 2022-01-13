source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
# https://r-charts.com/base-r/margins/


pdf("outputs/figures/Fig_envir_cover/test2.pdf",         # File name
    width = 8, height = 20, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")   # Color model (cmyk is required for most publications)
# paper = "A4")

# par(mar = c(1, 1, 0.1, 1),mfrow = c(2,1))

# par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0,mfrow=c(1,1)) # default

par(mar = c(5.1, 4.1, 4.1, 2.1) ,mfrow = c(2,1))


# Vegetation consumption ####
disturbance <- read.table("data/environment/Disturbance_DivHerbe.txt",header=T,sep="\t",dec=",") %>% 
  filter(!(Trtmt == "Tem"))

boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  
        width=c(1,4), 
        col=c("orange" , "seagreen"),
        xlab = NA,
        ylab = "Proportion of biomass eaten",
        xaxt = "n",
        # log = "x",
        at = c(1,2)
)



# INN ####
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

BarPlot <- barplot(to_barplot,
        beside = TRUE, names.arg = c("Fer", "D","I","S"),#names.arg = IN2$soil, # yaxt = "n",
        ylim=c(0, max(c(IN2$meanINN,IN2$meanINP)+3)  ),  col = 3:4 )
legend("topright",
       legend = c("INN (%)","INP (%)"),
       pch = 15,
       col = 3:4)
error.bar(BarPlot,to_barplot,to_barplot_se)

# axis(4, at = myRightAxisAt, labels = myRightAxisLabs)

dev.off()






# Example from tuto
####
#specify path to save PDF to
destination = "outputs/figures/Fig_envir_cover/test3.pdf"

#open PDF
pdf(file=destination)

#specify to save plots in 2x2 grid
par(mfrow = c(2,2))

#save plots to PDF
for (i in 1:4) {   
  x=rnorm(i)  
  y=rnorm(i)  
  plot(x, y)   
}

#turn off PDF plotting
dev.off() 



# Version pré-finale, à modifier sous inkscape ####
#specify path to save PDF to
destination = "outputs/figures/Fig_envir_cover/test3.pdf"

#open PDF
pdf(file=destination)

#specify to save plots in 2x2 grid
par(mar = c(2, 10, 0.1, 1),mfrow = c(2,1),mpg = c(1,1,0))

#save plots to PDF

boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  
        width=c(1,4), 
        col=c("orange" , "seagreen"),
        xlab = NA,
        ylab = "Proportion of biomass eaten",
        xaxt = "n",
        # log = "x",
        at = c(1,2)
)

BarPlot <- barplot(to_barplot,
                   beside = TRUE, names.arg = c("Fer", "D","I","S"),#names.arg = IN2$soil, # yaxt = "n",
                   ylim=c(0, max(c(IN2$meanINN,IN2$meanINP)+3)  ),  col = 3:4,
                   ylab = "Nutrition index (%)")
legend("topright",
       legend = c("INN (%)","INP (%)"),
       pch = 15,
       col = 3:4)
error.bar(BarPlot,to_barplot,to_barplot_se)

#turn off PDF plotting
dev.off() 
















