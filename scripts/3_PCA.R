source("scripts/1. Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)

# MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")%>%
#   filter(!is.na(SLA))

# keep only traits measured in the Nat_Sab
# = compare trait values in Nat_Sab and in fertile
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv")%>%
  filter(!is.na(SLA)) 

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",# "LPC",
            # "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
            "Disp"#,"Mat_Per", #"Mat","Flo",
            # "SeedMass"
)



annuals <- MEAN %>% 
  filter(LifeHistory=="annual") %>% 
  pull(code_sp) %>% 
  unique()



#________________________________________________________________________
perform_pca <- function(data_traits_for_PCA){
  ACP1<-PCA(data_traits_for_PCA2,graph = FALSE)
  factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30)) # percentage of variance explained
  
  coord_var <- data.frame(ACP1$var$coord) %>% 
    rownames_to_column()
  var.explain.dim1 <- round(ACP1$eig[1,2])
  var.explain.dim2 <- round(ACP1$eig[2,2])
  coord_ind <- data.frame(ACP1$ind$coord) %>% 
    rownames_to_column("Code_Sp") %>% 
    mutate(Lifelength = if_else(Code_Sp %in% annuals,"annual","perennial"))
  
  list(coord_ind,coord_var,var.explain.dim1,var.explain.dim2)
}

plot_pca <- function(coord_ind,coord_var,var.explain.dim1,var.explain.dim2){
  ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=Lifelength))+
    theme_classic() +
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    coord_equal() +
    geom_point(size=4) +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*5.5-0.2, yend=Dim.2*5.5-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=Dim.1*5.5, Dim.2*5.5, label=rowname), size = 6, vjust=1, color="black") +
    ggtitle(trtmt) +
    theme(legend.position = "none") +
    xlab(paste0("Dim1 (",var.explain.dim1,"%)"))+
    ylab(paste0("Dim2 (",var.explain.dim2,"%)")) + 
    theme(text=element_text(size=20)) 
    
}

boxplot_dimension <- function(coord_ind,dim){
  ggplot(coord_ind,aes_string(x = "Lifelength",y=dim,color = "Lifelength"))+
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                          map_signif_level = TRUE,vjust = 0.5,col="black") +
    ggtitle(dim)
}

plot_pca_boxplot <- function(PCA1,plot_d1,plot_d2){
  ggarrange(
    PCA1,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(plot_d1,plot_d2, ncol = 2, labels = c("B", "C")), 
    nrow = 2, 
    labels = "A"       # Label of the line plot
  ) 
}


compute_anova_dim_x <- function(coord_ind,dimension){
  mod_dim <- lm(coord_ind[,dimension + 1] ~ coord_ind[, length(colnames(coord_ind))])
  # predict score in dimension wanted as a function of LifeLength
}

info_anova_dim_x <- function(mod_dim,dimension){
  normality_test <- shapiro.test(residuals(mod_dim))
  homoscedasticity_test <- lmtest::bptest(mod_dim)
  independence_test <- lmtest::dwtest(mod_dim)
  
  anov <- anova(mod_dim)
  sum <- summary(mod_dim)
  
  list(normality_test,homoscedasticity_test,independence_test,
       anov,sum)
}

draw_curve <- function(x,a=113000,b=-1.58){ # to link SLA to LDMC (papier Eric)
  a*x^b
}


# 1) PCA comparing annuals and perennial ####

data_traits_for_PCA <- MEAN %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE))
# NB : des fois, certaines espèces ont le même code, mais pas le même nom d'espèce.
# Il faut que je règle ça en faisant la moyenne dès  le départ (dans le script sur les traits)

# i) Fertile ####
trtmt <- "Fer"
data_traits_for_PCA2 <- data_traits_for_PCA %>% 
  filter(treatment == trtmt) %>% 
  select(-treatment) %>% 
  column_to_rownames("code_sp")

pca_output <- perform_pca(data_traits_for_PCA2)
coord_ind <- pca_output[[1]]
coord_var <- pca_output[[2]]
var.explain.dim1 <- pca_output[[3]]
var.explain.dim2 <- pca_output[[4]]

PCA_fer <- plot_pca(coord_ind,coord_var,var.explain.dim1,var.explain.dim2 ) +
  ggtitle((expression(paste("G"^'+',"F",sep=''))))

PCA_sp_names <- ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
  ggrepel::geom_label_repel() # With species names

plot_d1 <- boxplot_dimension(coord_ind,dim = "Dim.1")
plot_d2 <- boxplot_dimension(coord_ind,dim = "Dim.2")
PCA_fer_boxplot <- plot_pca_boxplot(PCA_fer,plot_d1,plot_d2)
ggsave("figures/PCA_fertile.png",PCA_fer_boxplot,height = 20, width =20)


# ii) Natif ####
trtmt <- "Nat"
data_traits_for_PCA2 <- data_traits_for_PCA %>% 
  filter(treatment == trtmt) %>% 
  select(-treatment) %>% 
  column_to_rownames("code_sp")

pca_output <- perform_pca(data_traits_for_PCA2)
coord_ind <- pca_output[[1]]
coord_var <- pca_output[[2]]
var.explain.dim1 <- pca_output[[3]]
var.explain.dim2 <- pca_output[[4]]

PCA_nat <- plot_pca(coord_ind,coord_var,var.explain.dim1,var.explain.dim2 ) +
  ggtitle((expression(paste("GU"[S],sep=''))))

PCA_sp_names <- ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
  ggrepel::geom_label_repel() # With species names

plot_d1 <- boxplot_dimension(coord_ind,dim = "Dim.1")
plot_d2 <- boxplot_dimension(coord_ind,dim = "Dim.2")
PCA_nat_boxplot <- plot_pca_boxplot(PCA_nat,plot_d1,plot_d2)

ggsave("outputs/figures/PCA_natif.png",PCA_nat_boxplot,height = 20, width =20)

# iii) Both PCA graphs ####
PCA <- grid.arrange(PCA_fer,PCA_nat, ncol=2,heights = 100)
ggsave("draft/PCA_annuals_perennials.png",PCA,height = 20, width =20)


# iv) ANOVA Position on axes ####
dimension <- 2
anov_dim <- compute_anova_dim_x(coord_ind,dimension)

par(mfrow=c(2,2)) ; plot(anov_dim) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(anov_dim))) # normality_graph

# Compare scores on axis "dimension" between annuals and perennials
info_anova <- info_anova_dim_x(anov_dim,dimension)
info_anova[[1]] # normality
info_anova[[2]] # homoscedasticity
info_anova[[3]] # independence of residuals
info_anova[[4]] # anova
info_anova[[5]] # summary



#_______________________________________________________________
# 2) Annuals whoses traits were measured ####
# Reflects partly abundance, since the most abundant were measured (even in my case).
# Rq here on Fer and Nat_Sab

# J'ai enlevé des outliers (deux). Est-ce justifié ? Je pense que oui.
traits_pca_annuals <- MEAN %>% 
  filter(!(code_sp=="FILAPYRA")) %>% # /!\  outlier (lien SLA-dispersion) !
  filter(!(code_sp == "CREPVESI-HAE")) %>% # /!\  outlier
  filter(LifeHistory == "annual") %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) %>%  # /!\ JE POURRAIS NE PAS FAIRE DE MOYENNE, et prendre les données brutes !
  mutate(ddd = paste(code_sp,treatment,sep = "_")) %>% 
  column_to_rownames("ddd") %>% 
  select(-c(code_sp,treatment))

# compute CSR scores on these plants
forCSR_ann_database <- traits_pca_annuals %>%
  rownames_to_column("sp_tr") %>% 
  select(sp_tr,L_Area,LDMC,SLA) %>%
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) # change from mg/g to %
write.csv2(forCSR_ann_releves,
           "outputs/data/Pierce CSR/Traits_annuals_database.csv" ,row.names=F)

CSR_ann_database <- read.csv2("outputs/data/Pierce CSR/Traits_annuals_database_completed.csv" ) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric()) %>% 
  mutate(treatment = if_else(str_detect(sp_tr,"F"),"Fer","Nat"))

CSR_ann_database_gathered <- CSR_ann_database %>% 
  select(C,S,R,treatment) %>% 
  gather(key = score,value = value,-treatment) 

CSR_ann_database_gathered %>% 
  ggplot(aes(x= score,y=value, color = treatment)) +
  theme_classic()+
  geom_boxplot()   +
  scale_colour_manual(values=c("#009E73","#E69F00"),
                      labels = c(expression(paste("G"^'+',"F",sep='')),
                                 expression(paste("GU"[S],sep=''))) ) 
####
  
  
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

PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment))+ #,label = sp_tr
  theme_classic()+
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
  scale_colour_manual(values=c("#009E73","#E69F00"),
                      labels = c(expression(paste("G"^'+',"F",sep='')),
                                 expression(paste("GU"[S],sep=''))) )+
  scale_fill_manual('Program Type', values=c('pink','blue')) +
  labs(color = "Origin") +
  ggtitle("Species traits")

ggsave("draft/PCA_annuals.jpg",PCA2)

PCA_label <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment,label = sp_tr))+
  ggrepel::geom_label_repel() +
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() 
  # geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
  #              arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  # # geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  # ggtitle("Annuals in fer and nat sup") +
  # # theme(legend.position = "none")+
  # xlab(paste("Dim1",var.explain.dim1,"% variance explained"))+
  # ylab(paste("Dim2",var.explain.dim2,"% variance explained")) +
  # scale_colour_manual(values=c("#009E73","#E69F00")) 

# With species names
# ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
#   geom_text()
dimension <- 1
anov_dim <- compute_anova_dim_x(coord_ind,dimension)

par(mfrow=c(2,2)) ; plot(anov_dim) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(anov_dim))) # normality_graph

# Compare scores on axis "dimension" between annuals and perennials
info_anova <- info_anova_dim_x(anov_dim,dimension)
info_anova[[1]] # normality
info_anova[[2]] # homoscedasticity
info_anova[[3]] # independence of residuals
info_anova[[4]] # anova
info_anova[[5]] # summary
#___________________________________________________________________
# 3) Annuals in abundance relevés in Fer and Nat_Sab ####

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

# compute CSR scores on these plants
forCSR_ann_releves <- data_traits_for_PCA2 %>%
  rownames_to_column("sp_tr") %>% 
  select(sp_tr,L_Area,LDMC,SLA) %>%
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) # change from mg/g to %
write.csv2(forCSR_ann_releves,
           "outputs/data/Pierce CSR/Traits_annuals_releves.csv" ,row.names=F)

CSR_ann_releves <- read.csv2("outputs/data/Pierce CSR/Traits_annuals_releves_completed.csv" ) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric()) %>% 
  mutate(treatment = if_else(str_detect(sp_tr,"F"),"Fer","Nat"))

CSR_ann_releves_gathered <- CSR_ann_releves %>% 
  select(C,S,R,treatment) %>% 
  gather(key = score,value = value,-treatment) 

CSR_ann_releves_gathered %>% 
  ggplot(aes(x= score,y=value, color = treatment)) +
  theme_classic()+
  geom_boxplot()   +
  scale_colour_manual(values=c("#009E73","#E69F00"),
                    labels = c(expression(paste("G"^'+',"F",sep='')),
                               expression(paste("GU"[S],sep=''))) ) 

t.test(CSR_ann_releves$C ~ CSR_ann_releves$treatment)
t.test(CSR_ann_releves$S ~ CSR_ann_releves$treatment)
t.test(CSR_ann_releves$R ~ CSR_ann_releves$treatment)


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
PCA2

PCA_label <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment,label = sp_tr))+
  ggrepel::geom_label_repel() +
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


# ANOVA Position on axes
dimension <- 2
anov_dim <- compute_anova_dim_x(coord_ind,dimension)

par(mfrow=c(2,2)) ; plot(anov_dim) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(anov_dim))) # normality_graph

# Compare scores on axis "dimension" between annuals and perennials
info_anova <- info_anova_dim_x(anov_dim,dimension)
info_anova[[1]] # normality
info_anova[[2]] # homoscedasticity
info_anova[[3]] # independence of residuals
info_anova[[4]] # anova
info_anova[[5]] # summary

# Différenciation légère sur axe 1 (Lié notamment à la dispersion), et sur l'axe 2.

#______________________________________________________
# 4) GWM of annuals in both treatments ####
# Guild Weighted Mean
CWM_annuals_fer <- read.csv2("outputs/data/CWM_annuals_fer.csv" ) %>% 
  mutate(id_transect_quadrat = paste0("annual",id_transect_quadrat)) %>% 
  column_to_rownames("id_transect_quadrat") %>% 
  select(-c(paddock)) 
CWM_annuals_nat <- read.csv2("outputs/data/CWM_annuals_nat.csv" ) %>% 
  mutate(id_transect_quadrat = paste("annual",paddock,depth,line,sep="_")) %>% 
  column_to_rownames("id_transect_quadrat") %>% 
  select(-c(paddock,depth,line))

CWM_annuals <- rbind(CWM_annuals_fer,CWM_annuals_nat) %>% 
  select(all_of(traits))

CWM_perennials_fer <- read.csv2("outputs/data/CWM_perennials_fer.csv" ) %>% 
  mutate(id_transect_quadrat = paste0("perennial",id_transect_quadrat)) %>%
  column_to_rownames("id_transect_quadrat") %>% 
  select(-c(paddock))
CWM_perennials_nat <- read.csv2("outputs/data/CWM_perennials_nat.csv" ) %>% 
  mutate(id_transect_quadrat = paste("perennial",paddock,depth,line,sep="_")) %>% 
  column_to_rownames("id_transect_quadrat") %>% 
  select(-c(paddock,depth,line))

CWM_perennials <- rbind(CWM_perennials_fer,CWM_perennials_nat) %>% 
  select(all_of(traits)) 

CWM_all <- rbind(CWM_perennials,CWM_annuals)


# Plot CSR ####

rbind(CWM_perennials_fer,CWM_perennials_nat,CWM_annuals_fer,CWM_annuals_nat) %>% 
  rownames_to_column("transect") %>% 
  mutate(guild = if_else(str_detect(transect,"annual"),"annual","perennial")) %>%
  mutate(treatment = if_else(str_detect(transect,"F"),"Fer","Nat")) %>% 
  select(C,S,R,guild,treatment) %>% 
  gather(key = score,value = value,-c(guild,treatment)) %>% 
  ggplot(aes(x= score,y=value,color=guild)) +
    geom_boxplot() +
    facet_wrap(~ treatment)

# ACP ####
ACP1<-PCA(CWM_all,graph = FALSE)
factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained


coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("transect") %>% 
  mutate(treatment = if_else(str_detect(transect,"F"),"Fer","Nat")) %>% 
  mutate(guild = if_else(str_detect(transect,"annual"),"annual","perennial")) %>% 
  mutate(treatment_guild = paste(treatment,guild,sep="_"))

treatment <- coord_ind$treatment
guild <- coord_ind$guild


PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour= guild))+
  # scale_shape_manual(values = c(1,2))+
  theme_classic()+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point(aes(shape = treatment)) +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*5-0.2, yend=Dim.2*5-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  
  geom_text_repel(data=coord_var, aes(x=Dim.1*5, Dim.2*5, label=rowname), 
                  size = 4, vjust=1, color="black") +
  # theme(legend.position = "none")+
  xlab(paste("Dim1 (",var.explain.dim1,"%)")) +
  ylab(paste("Dim2 (",var.explain.dim2,"%)")) + 
  # scale_colour_manual( # values=c("#009E73","#E69F00"),
  #                      values = c("#F8766D","#00BFC4" ),
  #                     labels = c(expression(paste("G"^'+',"F",sep='')),
  #                                expression(paste("GU"[S],sep=''))) ) +
  scale_shape_manual(values = c(1,2 ),
                     labels = c(expression(paste("G"^'+',"F",sep='')),
                                expression(paste("GU"[S],sep=''))) ) +
  labs(shape = "Origin", color = "Guild") +
  ggtitle("Traits aggregated at the guild level")
PCA2

ggsave("draft/PCA_annuals_CWM.jpg",PCA2)




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


dimension <- 1
anov_dim <- compute_anova_dim_x(coord_ind,dimension)

par(mfrow=c(2,2)) ; plot(anov_dim) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(anov_dim))) # normality_graph

# Compare scores on axis "dimension" between annuals and perennials
info_anova <- info_anova_dim_x(anov_dim,dimension)
info_anova[[1]] # normality
info_anova[[2]] # homoscedasticity
info_anova[[3]] # independence of residuals
info_anova[[4]] # anova
info_anova[[5]] # summary


# Relationships SLA - LDMC ####

relat_sla_ldmc <- ggplot(MEAN,aes(x=LDMC,y=SLA,color = LifeHistory))+
  geom_point() +
  # facet_wrap(~treatment) +
  xlim(c(0,800)) +
  ylim(c(0,70)) +
  geom_function(fun = draw_curve,color="black") 
  facet_grid(vars(treatment),vars(LifeHistory)) 

ggsave("outputs/figures/Appendix/Relationship_SLA_LDMC.jpg",relat_sla_ldmc)
