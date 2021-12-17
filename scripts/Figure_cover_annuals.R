source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
# https://r-charts.com/base-r/margins/

par(mar = c(0.1, 1, 0.1, 1),
    mfrow = c(2,1))
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0,mfrow=c(1,1)) # default

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
  summarise(meanINN = mean(INN,na.rm = T),meanINP = mean(INP,na.rm=T))

myLeftAxisLabs <- pretty(seq(0, max(IN2$meanINN), length.out = 10))
myRightAxisLabs <- pretty(seq(0, max(IN2$meanINP), length.out = 10))
myLeftAxisAt <- myLeftAxisLabs/max(IN2$meanINN)
myRightAxisAt <- myRightAxisLabs/max(IN2$meanINP)

barplot(t(as.matrix(IN2[, c("meanINN", "meanINP")])),
        beside = TRUE, yaxt = "n", names.arg = IN2$soil,
        ylim=c(0, max(c(myLeftAxisLabs, myRightAxisLabs))))
axis(2, at = myLeftAxisLabs, labels = myLeftAxisLabs)
axis(4, at = myRightAxisLabs, labels = myRightAxisLabs)

legend("topright",
       legend = rownames(data),
       pch = 15,
       col = 1:nrow(data))



test <- data.frame(group = 1:10, var.a = rnorm(n = 10, mean = 500, sd = 20),
                   var.b = runif(10))
funProp <- function(testCol) {
  test[, testCol]/max(test[, testCol])
}
test$var.a.prop <- funProp("var.a")
test$var.b.prop <- funProp("var.b")

myLeftAxisLabs <- pretty(seq(0, max(test$var.a), length.out = 10))
myRightAxisLabs <- pretty(seq(0, max(test$var.b), length.out = 10))
myLeftAxisAt <- myLeftAxisLabs/max(test$var.a)
myRightAxisAt <- myRightAxisLabs/max(test$var.b)

barplot(t(as.matrix(test[, c("var.a.prop", "var.b.prop")])),
        beside = TRUE, yaxt = "n", names.arg = test$group,
        ylim=c(0, max(c(myLeftAxisAt, myRightAxisAt))))
axis(2, at = myLeftAxisAt, labels = myLeftAxisLabs)
axis(4, at = myRightAxisAt, labels = myRightAxisLabs)

# Vegetation consumption ####
disturbance <- read.table("data/environment/Disturbance_DivHerbe.txt",header=T,sep="\t",dec=",") %>% 
  filter(!(Trtmt == "Tem"))


boxplot(disturbance$Tx_CalcPic ~ disturbance$Trtmt ,  width=c(10,40), col=c("orange" , "seagreen"),
        xlab = NA,
        ylab = NA,
        xaxt = "n"
)





















