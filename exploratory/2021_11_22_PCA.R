source("scripts/1. Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
  # filter(!(Species == "Geranium dissectum - limbe")) %>% 
  # filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  # filter(!(Species == "Carex humilis?")) %>% 
  # filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer

# all the traits: 

# traits <- c("LDMC","SLA","L_Area",
#             "LCC","LNC","Ldelta13C","LPC",
#             "Hveg"  ,    "Hrepro"   , "Dmax"  ,    "Dmin" ,
#             "Mat","Flo","Disp","Mat_Per", 
#             "SeedMass"
# )

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C","LPC",
            "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
            "Flo","Disp","Mat_Per", #"Mat",
            "SeedMass"
)

#_______________________________________________________________________________
data_traits <- MEAN

# data_traits_for_PCA <- data_traits %>% 
#   filter(Trtmt=="Fer") %>%  # To compare with Garnier et al. 2018, where change in occurrence proba was computed on the G+F treatment
#   select(!!c("Code_Sp",traits)) %>% # subset of the traits that I want of analyse
#   column_to_rownames("Code_Sp")

# Version traits moyennés apr espèce sur les deux traitements
data_traits_for_PCA <- data_traits %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) 

#_____________________________________________________
# PCA ####
annuals <- MEAN %>% 
  filter(LifeHistory=="annual") %>% 
  pull(code_sp) %>% 
  unique()


# PCA on annuals f(pheno) ####
traits_selected_ann <- c("LDMC","SLA","L_Area",
                "LCC","LNC","LPC", #"Ldelta13C",
                "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
                # "Flo","Disp","Mat_Per", #"Mat",
                "SeedMass")

traits_ann <- MEAN %>% 
  filter(LifeHistory=="annual") %>% 
  select(!!c("code_sp","treatment",traits_selected_ann,"Disp","Mat")) %>% 
  group_by(code_sp) %>% 
  summarise(across(all_of(c(traits_selected_ann,"Disp","Mat")), mean, na.rm= TRUE)) 

ggplot(traits_ann,aes(x= Disp))+
  geom_density()

traits_ann %>% 
  filter(!(is.na(Disp))) %>%
  arrange(Disp) %>% 
  pull(code_sp)

data_PCA_ann <- traits_ann %>% 
  select(-c(Disp,Mat)) %>% 
  column_to_rownames("code_sp")

ACP1<-PCA(data_PCA_ann,graph = T)

coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
var.explain.dim3 <- round(ACP1$eig[3,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("Code_Sp") %>% 
  mutate(Disp = traits_ann$Disp) %>% 
  mutate(pheno = case_when(Disp<=160 ~"early",
                           Disp > 180 ~ "late",
                           TRUE ~ "middle"))
  

factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained

ggplot(coord_ind %>% filter(!is.na(pheno)),aes(x=Dim.1,y=Dim.2,colour=pheno))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained"))

  
ggplot(coord_ind %>% filter(!is.na(pheno)),aes(x=Dim.1,y=Dim.3,colour=pheno))+
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    coord_equal() +
    geom_point() +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.3*7-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.3*7, label=rowname), size = 4, vjust=1, color="black") +
    xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
    ylab(paste("Dim3",var.explain.dim3,"% variance explained"))

ggplot(coord_ind %>% filter(!is.na(pheno)),aes(x=Dim.2,y=Dim.3,colour=pheno))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.2*7-0.2, yend=Dim.3*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.2*7, Dim.3*7, label=rowname), size = 4, vjust=1, color="black") +
  xlab(paste("Dim2",var.explain.dim2,"% variance explained"))+
  ylab(paste("Dim3",var.explain.dim3,"% variance explained"))

#_______________________________________________________________________________
trtmt <- "Fer"
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
  
  ggsave("outputs/figures/PCA_fertile.png",PCA_fer,height = 20, width =20)
  
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

  ggsave("outputs/figures/PCA_natif.png",PCA_nat,height = 20, width =20)
}





# Both PCA graphs ####
# grid.arrange(PCA1,PCA2,PCA1Di,PCA2Di, ncol=2)
PCA <- grid.arrange(PCA1,PCA2, ncol=2)
ggsave("outputs/figures/PCA.png",PCA,height = 20, width =20)


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

#________________________________________________________________________________
# II) PCA on CWM ####
# Adeline ####
Adeline_ab_tr <- read.csv2("outputs/data/pooled_abundance_and_traits.csv") %>% 
  filter(dataset == "Adeline") %>% 
  group_by(paddock,id_transect_quadrat) # NB choose the level at which to compute moments. group, or  plot...

CWM_adeline0 <- Adeline_ab_tr %>% 
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()

CWM_adeline <- CWM_adeline0 %>% 
  select(contains("CWM")) %>% 
  unique() %>% 
  rename_at( vars( contains( "CWM_") ), list( ~paste(gsub("CWM_", "", .), sep = "_") ) ) %>% 
  mutate(Trtmt = case_when(
    str_detect(paddock, "C") ~ "Fer",
    str_detect(paddock, "N") ~ "Nat",
    str_detect(paddock, "T") ~ "Tem")) 


treatment <- CWM_adeline$Trtmt

trtmt="Fer"
data_traits_for_PCA2 <- CWM_adeline %>% 
  ungroup() %>% 
  select(-c(Trtmt,id_transect_quadrat,paddock)) %>% 
  select(all_of(traits))

ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column("trait")


var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  mutate(Trtmt = treatment)



CWM_Adeline <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2, color = Trtmt))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=trait), size = 4, vjust=1, color="black") +
  # theme(legend.position = "none") +
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) +
  ggtitle("Adeline CWM") +
  scale_color_manual(values=c("green", "orange", "black"))

ggsave("outputs/figures/PCA_CWM_Adeline.png",CWM_Adeline)
