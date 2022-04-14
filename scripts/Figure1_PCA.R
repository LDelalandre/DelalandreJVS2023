source("scripts/1. Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)

# MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")%>%
#   filter(!is.na(SLA))
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv")%>%
  filter(!is.na(SLA)) # keep only traits measured in the Nat_Sab

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",# "LPC",
            # "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
            "Flo","Disp"#,"Mat_Per", #"Mat",
            # "SeedMass"
)


data_traits_for_PCA <- MEAN %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE))
# NB : des fois, certaines espèces ont le même code, mais pas le même nom d'espèce.
# Il faut que je règle ça en faisant la moyenne dès  le départ (dans le script sur les traits)

annuals <- MEAN %>% 
  filter(LifeHistory=="annual") %>% 
  pull(code_sp) %>% 
  unique()

# Biplots ####
draw_curve <- function(x,a=113000,b=-1.58){
  a*x^b
}


ggplot(MEAN,aes(x=LDMC,y=SLA,color = LifeHistory))+
  geom_point() +
  # facet_wrap(~treatment) +
  xlim(c(0,800)) +
  ylim(c(0,70)) 
  geom_function(fun = draw_curve,color="black")

#________________________________________________________________________
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
    geom_point(size=4) +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 6, vjust=1, color="black") +
    ggtitle(trtmt) +
    theme(legend.position = "none") +
    xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
    ylab(paste("Dim2",var.explain.dim2,"% variance explained")) + 
    theme(text=element_text(size=20))
  
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
    geom_point( size = 4) +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 6, vjust=1, color="black") +
    ggtitle(trtmt) +
    theme(legend.position = "none")+
    xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
    ylab(paste("Dim2",var.explain.dim2,"% variance explained")) + 
    theme(text=element_text(size=20))
  
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

#_______________________________________________________________
# Annuals only ####
traits_pca_annuals <- MEAN %>% 
  filter(LifeHistory == "annual") %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) %>%  # JE POURRAIS NE PAS FAIRE DE SUMMARY !
  mutate(ddd = paste(code_sp,treatment,sep = "_")) %>% 
  column_to_rownames("ddd") %>% 
  select(-c(code_sp,treatment))
  
ACP1<-PCA(traits_pca_annuals,graph = FALSE)
factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained


coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("sp_tr") %>% 
  mutate(code_sp = str_sub(sp_tr,1L,-5L)) %>% 
  mutate(treatment = str_sub(sp_tr,-3L,-1L))

PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment,label = sp_tr))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) +
  scale_colour_manual(values=c("#009E73","#E69F00")) 

PCA_label <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment,label = sp_tr))+
  ggrepel::geom_label_repel()+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  # geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) +
  scale_colour_manual(values=c("#009E73","#E69F00")) 

# With species names
# ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
#   geom_text()

#___________________________________________________________________
# Annuals in both treatments ####
# traits <- c("LDMC","SLA","L_Area",
#             "LCC","LNC","Ldelta13C","LPC",
#             "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
#             "Flo","Disp","Mat_Per", #"Mat",
#             "SeedMass"
# )

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
# Prenons l'abondance dans le diachro sur plusieurs années pour voir si je capte plus d'sp
# Fertilized treatment
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))

ab_diachro <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year %in% c(2000,2001,2002,2003,2004,2005)) %>%  # /!\ No measure available after 2005! Should I keep 2006, or several years?
  filter(treatment == "Fer") %>% 
  # add relative abundance
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance))

# abondance dans le natif
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

ann_fer <- ab_fer %>% 
  filter(LifeForm1=="The") %>% 
  pull(code_sp) %>% 
  unique()
ann_nat <- ab_nat %>% 
  filter(LifeHistory=="annual") %>% 
  # filter(depth == "S") %>%
  pull(code_sp) %>% 
  unique()

# Species replacement
common <- intersect(ann_fer,ann_nat)
justfer <- setdiff(ann_fer,ann_nat) # species in fer not in nat_sup
justnat <- setdiff(ann_nat,ann_fer) # sp in nat_sup not in fer

# PCA
data_traits_for_PCA <- MEAN %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) %>% 
  filter(code_sp %in% c(justfer,justnat,common))
  
# Plot SLA LDMC
ggplot(data_traits_for_PCA,aes(x=LDMC,y=SLA,color = treatment))+
  geom_point() +
  # facet_wrap(~treatment) +
  xlim(c(0,800)) +
  ylim(c(0,70)) +
  scale_colour_manual(values=c("#009E73","#E69F00")) +
  geom_function(fun = draw_curve,color="black")

ggplot(data_traits_for_PCA,aes(x=LNC,y=SLA,color = treatment))+
  geom_point() +
  scale_colour_manual(values=c("#009E73","#E69F00"))

data_traits_for_PCA_fer <- data_traits_for_PCA %>% 
  ungroup() %>% 
  # filter(code_sp %in% annuals) %>%
  filter( code_sp %in% ann_fer & treatment == "Fer" ) %>%
  mutate(sp_tr = paste(code_sp,treatment,sep="_")) %>% 
  select(-c(treatment,code_sp)) %>% 
  column_to_rownames("sp_tr")

data_traits_for_PCA_nat <- data_traits_for_PCA %>% 
  ungroup() %>% 
  # filter(code_sp %in% annuals) %>%
  filter( code_sp %in% ann_nat & treatment == "Nat" ) %>%
  mutate(sp_tr = paste(code_sp,treatment,sep="_")) %>% 
  select(-c(treatment,code_sp)) %>% 
  column_to_rownames("sp_tr")

data_traits_for_PCA2 <- rbind(data_traits_for_PCA_fer,data_traits_for_PCA_nat)


ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained


coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("sp_tr") %>% 
  mutate(code_sp = str_sub(sp_tr,1L,-5L)) %>% 
  mutate(treatment = str_sub(sp_tr,-3L,-1L))

PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) +
  scale_colour_manual(values=c("#009E73","#E69F00")) 

PCA_label <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment,label = sp_tr))+
  ggrepel::geom_label_repel()+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  # geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) +
  scale_colour_manual(values=c("#009E73","#E69F00")) 

# With species names
# ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
#   geom_text()

plot_d1 <- ggplot(coord_ind,aes(x = treatment,y=Dim.1,color = treatment))+
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle("Dimension 1")


plot_d2 <- ggplot(coord_ind,aes(x = treatment,y=Dim.2,color = treatment))+
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle("Dimension 2")


PCA_annuals <- ggpubr::ggarrange(
  PCA2,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(plot_d1,plot_d2, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A"       # Label of the line plot
) 

PCA_annuals
PCA2






#______________________________________________________
# CWM of annuals in both treatments ####
CWM_annuals_fer <- read.csv2("outputs/data/CWM_annuals_fer.csv" ) %>% 
  column_to_rownames("id_transect_quadrat") %>% 
  select(-c(paddock))
CWM_annuals_nat <- read.csv2("outputs/data/CWM_annuals_nat.csv" ) %>% 
  mutate(id_transect_quadrat = paste(paddock,depth,line,sep="_")) %>% 
  column_to_rownames("id_transect_quadrat") %>% 
  select(-c(paddock,depth,line))

CWM_annuals <- rbind(CWM_annuals_fer,CWM_annuals_nat) %>% 
  select(all_of(traits))





ACP1<-PCA(CWM_annuals,graph = FALSE)
factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained


coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("transect") %>% 
  mutate(treatment = if_else(str_detect(transect,"F"),"Fer","Nat"))

treatment <- coord_ind$teatment


PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  ylab(paste("Dim2",var.explain.dim2,"% variance explained")) + 
  scale_colour_manual(values=c("#009E73","#E69F00"))

# With species names
# ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
#   geom_text()

plot_d1 <- ggplot(coord_ind,aes(x = treatment,y=Dim.1,color = treatment))+
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle("Dimension 1")


plot_d2 <- ggplot(coord_ind,aes(x = treatment,y=Dim.2,color = treatment))+
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle("Dimension 2")


PCA_annuals <- ggpubr::ggarrange(
  PCA2,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(plot_d1,plot_d2, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A"       # Label of the line plot
) 

PCA_annuals
PCA2
