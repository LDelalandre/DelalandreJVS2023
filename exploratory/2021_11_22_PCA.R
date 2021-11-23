source("scripts/1. Packages.R")
source("scripts/2. Import files.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer
  
traits <- c("LDMC","SLA","LCC","LNC","Ldelta13C","LPC","SeedMass","Mat","Flo","Disp","Mat_Per", "L_Area", 
            "Hveg"  ,    "Hrepro"   , "Dmax"  ,    "Dmin" )
c1 <-
  MEAN %>% 
  mutate(sp_tr = paste(Code_Sp,Trtmt, sep = "_")) %>% 
  column_to_rownames("sp_tr") %>% 
  ungroup() %>% 
  select(all_of(traits)) 

ACP1 <- FactoMineR::PCA(c1,graph = FALSE)
factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30))


c1$Dim.1 <- ACP1$ind$coord[, 1] 
c1$Dim.2 <- ACP1$ind$coord[, 2]
c1$Dim.3 <- ACP1$ind$coord[, 3] 
c1$SName_trtmt <- rownames(c1)

c1[c('SName', 'Trtmt')] <- str_split_fixed(c1$SName_trtmt, '_', 2)


axis <- ACP1$var$coord[,c(1,2)] %>% 
  data.frame() %>% 
  rownames_to_column(var="varnames")

axis2 <- transform(axis,
                   Dim.1 = 4.7 * Dim.1,
                   Dim.2 = 4.7* Dim.2)

var.explain.dim1 <- round(ACP1$eig[1,2],digits=1)
var.explain.dim2 <- round(ACP1$eig[2,2],digits=1)
var.explain.dim2 <- round(ACP1$eig[3,2],digits=1)


ggplot(data=c1,aes(x=Dim.1,y=Dim.2)) + 
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() + 
  geom_text(data=axis2, aes(x=Dim.1, Dim.2, label=varnames), size = 5, vjust=1, color="black")+
  geom_segment(data=axis2, aes(x=0, y=0, xend=Dim.1-0.2, yend=Dim.2-0.2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  # ggrepel::geom_label_repel(aes(label = SName),size = 5) + # or SName +
  geom_point() +
  scale_colour_gradient(low="#00BFC4",high="#F8766D") +
  # scale_colour_gradient(low = "#132B43",high = "#56B1F7" ) +
  labs(x=paste0("Dim 1 (",var.explain.dim1,"%)"),
       y=paste0("Dim 2 (",var.explain.dim2,"%)") ) +
  theme_classic() +
  theme(axis.title=element_text(size=15),axis.text=element_text(size=15)) +
  facet_wrap(~Trtmt)



#_______________________________________________________________________________
# Version comité ####
data_traits <- read.csv2("data/processed/2021_03/all traits_mean.csv")
traits <- c("LDMC","SLA","LCC","LNC","Ldelta13C","LPC","SeedMass","Disp","Mat_Per","Hrepro") # "Mat","Flo","Hveg",
# data_traits_for_PCA <- data_traits %>% 
#   filter(Trtmt=="Fer") %>%  # To compare with Garnier et al. 2018, where change in occurrence proba was computed on the G+F treatment
#   select(!!c("Code_Sp",traits)) %>% # subset of the traits that I want of analyse
#   column_to_rownames("Code_Sp")

# Version traits moyennés apr espèce sur les deux traitements
data_traits_for_PCA <- data_traits %>% 
  select(!!c("Code_Sp","Trtmt",traits)) %>%  # subset of the traits that I want of analyse
  group_by(Code_Sp,Trtmt) %>% 
  summarise(across(LDMC:Hrepro, mean, na.rm= TRUE)) 

#_____________________________________________________
# PCA ####
sp <- read.csv2("data/processed/2021_02_annuals/cle_taxo_withLifeForm.csv")
annuals <- sp %>% 
  filter(LifeForm1=="The") %>% 
  pull(CODE_ESP)

# Fertile
trtmt <- "Fer"

data_traits_for_PCA2 <- data_traits_for_PCA %>% 
  filter(Trtmt == trtmt) %>% 
  select(-Trtmt) %>% 
  column_to_rownames("Code_Sp")

ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("Code_Sp") %>% 
  mutate(Lifelength = if_else(Code_Sp %in% annuals,"annual","perennial"))

Lifelength <- coord_ind$Lifelength
forDi <- coord_ind %>% 
  select(-Lifelength) %>% 
  column_to_rownames("Code_Sp")

Dist <- traits_dist(forDi)  %>% 
  rownames_to_column() %>%
  rename(Code_Sp=rowname) %>% 
  mutate(Lifelength = Lifelength)

PCA1 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=Lifelength))+
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
  ylab(paste("Dim1",var.explain.dim2,"% variance explained"))

PCA1Di <- ggplot(Dist,aes(x=Dim.1,y=Dim.2,shape=Lifelength,colour = Di))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  # geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
  #              arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  # geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle(trtmt) +
  scale_color_gradientn(colours = rainbow(5))


# Natif
trtmt <- "Nat"

data_traits_for_PCA2 <- data_traits_for_PCA %>% 
  filter(Trtmt == trtmt) %>% 
  select(-Trtmt) %>% 
  column_to_rownames("Code_Sp")

ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("Code_Sp") %>% 
  mutate(Lifelength = if_else(Code_Sp %in% annuals,"annual","perennial"))

Lifelength <- coord_ind$Lifelength
forDi <- coord_ind %>% 
  select(-Lifelength) %>% 
  column_to_rownames("Code_Sp")

Dist <- traits_dist(forDi)  %>% 
  rownames_to_column() %>%
  rename(Code_Sp=rowname) %>% 
  mutate(Lifelength = Lifelength)

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
  ylab(paste("Dim1",var.explain.dim2,"% variance explained"))

PCA2Di <- ggplot(Dist,aes(x=Dim.1,y=Dim.2,shape=Lifelength,colour = Di))+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  # geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
  #              arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  # geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle(trtmt) +
  scale_color_gradientn(colours = rainbow(5))

# Both PCA graphs
# grid.arrange(PCA1,PCA2,PCA1Di,PCA2Di, ncol=2)
grid.arrange(PCA1,PCA2, ncol=2)
