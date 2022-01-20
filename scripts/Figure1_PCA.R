source("scripts/1. Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C","LPC",
            "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
            "Flo","Disp","Mat_Per", #"Mat",
            "SeedMass"
)


data_traits_for_PCA <- MEAN %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) 



trtmt <- "Nat"
# Fertile ####
if(trtmt == "Fer"){
  data_traits_for_PCA2 <- data_traits_for_PCA %>% 
    filter(treatment == trtmt) %>% 
    select(-treatment) %>% 
    column_to_rownames("code_sp")
  
  ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
  factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30)) # percentage of variance explained
  
  coord_var <- data.frame(ACP1$var$coord) %>% 
    rownames_to_column()
  var.explain.dim1 <- round(ACP1$eig[1,2])
  var.explain.dim2 <- round(ACP1$eig[2,2])
  coord_ind <- data.frame(ACP1$ind$coord) %>% 
    rownames_to_column("Code_Sp") %>% 
    mutate(Lifelength = if_else(Code_Sp %in% annuals,"annual","perennial"))
  
  Lifelength <- coord_ind$Lifelength
  
  PCA1 <-
    ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=Lifelength))+
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    coord_equal() +
    geom_point() +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
    ggtitle(trtmt) +
    theme(legend.position = "none") +
    xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
    ylab(paste("Dim2",var.explain.dim2,"% variance explained"))
  
  # With species names
  # ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
  #   geom_text()
  
  plot_d1 <- ggplot(coord_ind,aes(x = Lifelength,y=Dim.1,color = Lifelength))+
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                          map_signif_level = TRUE,vjust = 0.5,col="black") +
    ggtitle("Dimension 1")
  
  
  plot_d2 <- ggplot(coord_ind,aes(x = Lifelength,y=Dim.2,color = Lifelength))+
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                          map_signif_level = TRUE,vjust = 0.5,col="black") +
    ggtitle("Dimension 2")
  
  
  PCA_fer <- ggarrange(
    PCA1,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(plot_d1,plot_d2, ncol = 2, labels = c("B", "C")), 
    nrow = 2, 
    labels = "A"       # Label of the line plot
  ) 
  
  ggsave("outputs/plots/PCA_fertile.png",PCA_fer,height = 20, width =20)
  
} else if (trtmt == "Nat"){
  # Natif ####
  trtmt <- "Nat"
  
  data_traits_for_PCA2 <- data_traits_for_PCA %>% 
    filter(treatment == trtmt) %>% 
    select(-treatment) %>% 
    column_to_rownames("code_sp")
  
  ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
  factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30)) # percentage of variance explained
  
  
  coord_var <- data.frame(ACP1$var$coord) %>% 
    rownames_to_column()
  var.explain.dim1 <- round(ACP1$eig[1,2])
  var.explain.dim2 <- round(ACP1$eig[2,2])
  coord_ind <- data.frame(ACP1$ind$coord) %>% 
    rownames_to_column("Code_Sp") %>% 
    mutate(Lifelength = if_else(Code_Sp %in% annuals,"annual","perennial"))
  
  Lifelength <- coord_ind$Lifelength
  # DISTINCTIVENESS measured as in ForCEEPS paper
  # forDi <- coord_ind %>% 
  #   select(-Lifelength) %>% 
  #   column_to_rownames("Code_Sp")
  # 
  # Dist <- traits_dist(forDi)  %>% 
  #   rownames_to_column() %>%
  #   rename(Code_Sp=rowname) %>% 
  #   mutate(Lifelength = Lifelength)
  
  coord_ind <- coord_ind %>% filter(!(Code_Sp %in% c("QUERPUBE","BUXUSEMP","CRATMONO","PRUNMAHA")))
  PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=Lifelength))+
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    coord_equal() +
    geom_point() +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
    ggtitle(trtmt) +
    theme(legend.position = "none")+
    xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
    ylab(paste("Dim2",var.explain.dim2,"% variance explained"))
  
  # With species names
  # ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
  #   geom_text()
  
  plot_d1 <- ggplot(coord_ind,aes(x = Lifelength,y=Dim.1,color = Lifelength))+
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                          map_signif_level = TRUE,vjust = 0.5,col="black") +
    ggtitle("Dimension 1")
  
  
  plot_d2 <- ggplot(coord_ind,aes(x = Lifelength,y=Dim.2,color = Lifelength))+
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                          map_signif_level = TRUE,vjust = 0.5,col="black") +
    ggtitle("Dimension 2")
  
  
  PCA_nat <- ggpubr::ggarrange(
    PCA2,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(plot_d1,plot_d2, ncol = 2, labels = c("B", "C")), 
    nrow = 2, 
    labels = "A"       # Label of the line plot
  ) 
  
  ggsave("outputs/plots/PCA_natif.png",PCA_nat,height = 20, width =20)
}





# Both PCA graphs ####
# grid.arrange(PCA1,PCA2,PCA1Di,PCA2Di, ncol=2)
PCA <- grid.arrange(PCA1,PCA2, ncol=2)
ggsave("outputs/plots/PCA.png",PCA,height = 20, width =20)


# ANOVA Position on axes ####
# NB: choose natif or fertile by running the corresponding code

# Axis 1
comp.dim.1 <- lm(Dim.1~Lifelength,data=coord_ind)
par(mfrow=c(2,2)) ; plot(comp.dim.1)

shapiro.test(residuals(comp.dim.1))
lmtest::bptest(comp.dim.1)
lmtest::dwtest(comp.dim.1)

par(mfrow= c(1,1)) ; plot(density(residuals(comp.dim.1)))

anova(comp.dim.1)
summary(comp.dim.1)

# Axis 2 not gignificant (and not normal...) in the natif
comp.dim.2 <- lm(Dim.2~Lifelength,data=coord_ind)
par(mfrow=c(2,2)) ; plot(comp.dim.2)

shapiro.test(residuals(comp.dim.2))
lmtest::bptest(comp.dim.2)
lmtest::dwtest(comp.dim.2)

par(mfrow= c(1,1)) ; plot(density(residuals(comp.dim.2)))

anova(comp.dim.2)
summary(comp.dim.2)