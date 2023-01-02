source("scripts/Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)

# keep only traits measured in the Nat_Sab
# = compare trait values in Nat_Sab and in fertile
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  filter(!species == "Geranium dissectum - pÃ©tiole")

code_sp_lifeform <- read.csv2("data/species_names_lifehistory.csv")


#_____________________________________________

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "H_FLORE",#"Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "FLO_FLORE", #Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)

MEAN %>% 
  # filter(treatment=="Nat") %>%
  select(all_of(traits)) %>% 
  psych::corPlot()


annuals <- MEAN %>% 
  filter(LifeHistory=="annual") %>% 
  pull(code_sp) %>% 
  unique()



#________________________________________________________________________
perform_pca <- function(data_traits_for_PCA){
  ACP1<-PCA(data_traits_for_PCA,graph = FALSE)
  # NB: scale.unit = TRUE by default. Data are thus scaled.
  plot_var_explained <- factoextra::fviz_eig(ACP1, addlabels = TRUE, ylim = c(0, 30)) # percentage of variance explained
  
  coord_var <- data.frame(ACP1$var$coord) %>% 
    rownames_to_column()
  var.explain.dim1 <- round(ACP1$eig[1,2])
  var.explain.dim2 <- round(ACP1$eig[2,2])
  var.explain.dim3 <- round(ACP1$eig[3,2])
  var.explain.dim4 <- round(ACP1$eig[4,2])
  coord_ind <- data.frame(ACP1$ind$coord) %>% 
    rownames_to_column("Code_Sp") %>% 
    mutate(Lifelength = if_else(Code_Sp %in% annuals,"annual","perennial"))
  
  list(coord_ind,coord_var,var.explain.dim1,var.explain.dim2,var.explain.dim3,var.explain.dim4,plot_var_explained)
}

plot_pca <- function(coord_ind,coord_var,DimA,DimB,var.explain.dimA,var.explain.dimB){
  # dimA and dimB can be the first and second dimension, or the first and third, etc.
  ggplot(coord_ind,aes_string(x=DimA,y=DimB), width = 10, height = 10)+
    theme_classic() +
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    coord_equal() +
    # geom_point(size=4,aes(shape = Lifelength,colour=Lifelength)) +
    geom_point(size=4,aes(shape = LifeForm1,colour=LifeForm1)) +
    geom_segment(data=coord_var, aes(x=0, y=0, xend=get(DimA)*5.5-0.2, yend=get(DimB)*5.5-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=get(DimA)*5.5, get(DimB)*5.5, label=rowname), size = 6, vjust=1, color="black") +
    ggtitle(trtmt) +
    # theme(legend.position = "none") +
    xlab(paste0(DimA," (",var.explain.dimA,"%)"))+
    ylab(paste0(DimB," (",var.explain.dimB,"%)")) + 
    theme(text=element_text(size=20)) 
    
}


# ggplot(coord_var,aes(x=get("Dim.1"),y=Dim.2))+
#   geom_point()
# ggplot(data=coord_var, aes_string(x=0, y=0, xend=Dim.1*5.5-0.2, yend=Dim.2*5.5-0.2), 
#              arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
#   geom_point()


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
  mod_dim <- lm(coord_ind[,dimension + 1] ~ coord_ind[, length(colnames(coord_ind))-1])
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

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(depth == "S")

data_traits_for_PCA <- MEAN %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) 
# NB : des fois, certaines espèces ont le même code, mais pas le même nom d'espèce.
# Il faut que je règle ça en faisant la moyenne dès  le départ (dans le script sur les traits)

## SI ANALYSE DE SENSIBILITE AUX ESPECES ####
# pour analyse avec sp dans relevés bota

data_fer <- data_traits_for_PCA %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- data_traits_for_PCA %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
data_traits_for_PCA <- rbind(data_fer,data_nat)



# Clustering ####
# https://www.statology.org/k-means-clustering-in-r/
library("cluster")
library(factoextra)

ftreatment <- "Nat"

data_clustering <- data_traits_for_PCA %>% 
  filter(treatment == ftreatment) %>% 
  # mutate(sp_trtmt = paste(code_sp,treatment,sep="_")) %>% 
  column_to_rownames("code_sp") %>% 
  select(-c(treatment)) %>% 
  select(SLA,LDMC,L_Area,LCC,LNC)
  # select(-c(Ldelta13C,Hrepro,Dmax))

# clustering sur dimensions d'ACP (et ensuite le faire sur traits bruts en annexe)
# data_clustering <- coord_ind %>% 
#   column_to_rownames("Code_Sp") %>% 
#   select(-c(Lifelength,LifeForm1,cluster))

data_clustering2 <- data_clustering %>% 
  na.omit() %>%
  scale()

# data_clustering2 <- data_traits_for_PCA2 # celui où j'ai ajouté les moyennes par colonne

# find the optiomal number of clusters (2 examples)
dev.off()
fviz_nbclust(data_clustering2, kmeans, method = "wss")

gap_stat <- clusGap(data_clustering2,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
fviz_gap_stat(gap_stat)



#make this example reproducible
set.seed(1)

#perform k-means clustering with k = 2 clusters
km <- kmeans(data_clustering2, centers = 2, nstart = 25)

#view results
data_clusters <- km$cluster%>% 
  as.data.frame() %>% 
  rownames_to_column("code_sp")
colnames(data_clusters) <- c("code_sp", "cluster")
data_clusters <- data_clusters%>% 
  mutate(cluster = as.factor(cluster))

fviz_cluster(km, data = data_clustering2) 

lifeform_per_cluster <- coord_ind %>% 
  select(Code_Sp,Lifelength) %>% 
  rename(code_sp = Code_Sp) %>% 
  merge(data_clusters)
# nb of life hisotries in each cluster
cont_table <- table( lifeform_per_cluster$Lifelength,lifeform_per_cluster$cluster)
chisq.test(cont_table)

# Fisher’s exact test is an alternative to chi-squared test used mainly when a chi-square 
# approximation is not satisfactory (i.e., small sample size, and you get the warning message) 
fisher.test(cont_table)


## SI ANALYSE SENSIBILITE AUX TRAITS ####

## i) Fertile ####
trtmt <- "Fer"
data_traits_for_PCA2 <-data_traits_for_PCA %>% 
  filter(treatment == trtmt) %>%
  select(-treatment) %>% 
  column_to_rownames("code_sp")

for(i in 1:ncol(data_traits_for_PCA2)){
  data_traits_for_PCA2[is.na(data_traits_for_PCA2[,i]), i] <- 
    mean(data_traits_for_PCA2[,i], na.rm = TRUE)
}


# PCA avec fer et nat
# data_traits_for_PCA2 <- data_traits_for_PCA %>% 
#   ungroup() %>% 
#   mutate(pop = paste(code_sp,treatment,sep="_")) %>% 
#   select(-c(code_sp,treatment)) %>% 
#   column_to_rownames("pop")

# coord_ind <- pca_output[[1]] %>% 
#   separate(col = Code_Sp, into = c("Code_sp","treatment"), sep = "_") %>% 
#   merge(code_sp_lifeform) %>% 
#   mutate(LifeLenght = LifeForm1)


pca_output <- perform_pca(data_traits_for_PCA2)
# pca_output[[7]]
coord_ind <- pca_output[[1]] %>% 
  merge(code_sp_lifeform,by="Code_Sp") %>% 
  merge(data_clusters %>% rename(Code_Sp = code_sp))
  
coord_var <- pca_output[[2]]
var.explain.dim1 <- pca_output[[3]]
var.explain.dim2 <- pca_output[[4]]
var.explain.dim3 <- pca_output[[5]]
var.explain.dim4 <- pca_output[[6]]

Dim.A <- "Dim.1"
Dim.B <- "Dim.2"
Var.A <- var.explain.dim1
Var.B <- var.explain.dim2

# Faire des graphes de même taille
etendue_dim <- coord_ind %>% 
  gather(key=dim,value=coordinate,-c(Code_Sp,Lifelength)) %>% 
  group_by(dim) %>% 
  summarize(min = min(coordinate),max=max(coordinate)) %>% 
  filter(!(dim == "LifeForm1")) %>%
  mutate(min = as.numeric(min),max = as.numeric(max)) %>% 
  mutate(etendue= max-min)

ratio <- etendue_dim %>% 
  select(dim,etendue) %>% 
  spread(key = dim, value= etendue) %>% 
  mutate(ratio12=Dim.1/Dim.2,
         ratio13 = Dim.1/Dim.3)
dev.off()
PCA_fer12 <- plot_pca(coord_ind,coord_var,DimA=Dim.A,DimB=Dim.B ,Var.A,Var.B ) +
  ggtitle((expression(paste("G"^'+',"F",sep='')))) +
  # coord_fixed(ratio = ratio$ratio12) +
  # ggrepel::geom_label_repel(aes(label=Code_Sp)) + theme(legend.position = "none")
  # ellipse
  ggforce::geom_mark_ellipse(data = coord_ind,aes(fill = cluster), # ,label = cluster
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm')) +
  scale_fill_brewer() 
  # theme(legend.position = "none")
  # scale_fill_manual( values = c("#9933FF","#33FFFF" ) )

PCA_fer13 <- plot_pca(coord_ind,coord_var,DimA=Dim.A,DimB= "Dim.3" ,Var.A,var.explain.dim3 ) +
  ggtitle((expression(paste("G"^'+',"F",sep=''))))+
  # coord_fixed(ratio = ratio$ratio13) +
  ggforce::geom_mark_ellipse(data = coord_ind,aes(fill = cluster), # ,label = cluster
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm')) +
  scale_fill_brewer() 
  

PCA_fer12
PCA_fer13

PCA_sp_names <- ggplot (coord_ind,aes(x=Dim.1,y=Dim.2,label = Code_Sp,colour=Lifelength))+
  ggrepel::geom_label_repel() # With species names

plot_d1 <- boxplot_dimension(coord_ind,dim = "Dim.1")
plot_d2 <- boxplot_dimension(coord_ind,dim = "Dim.2")
PCA_fer_boxplot <- plot_pca_boxplot(PCA_fer12,plot_d1,plot_d2)
ggsave("output/figures/PCA_fertile.png",PCA_fer_boxplot,height = 20, width =20)

coord_ind_fer <- coord_ind



# nb of life hisotries in each cluster
cont_table <- table(coord_ind$cluster,coord_ind$Lifelength)
chisq.test(cont_table)

# Fisher’s exact test is an alternative to chi-squared test used mainly when a chi-square 
# approximation is not satisfactory (i.e., small sample size, and you get the warning message) 
fisher.test(cont_table)



## ii) Natif ####
trtmt <- "Nat"
data_traits_for_PCA2 <- data_traits_for_PCA %>% 
  filter(treatment == trtmt) %>% 
  select(-treatment) %>% 
  # filter(!(code_sp %in%  c("BROMEREC","POTENEUM"))) %>%  # Ces deux espèces tirent fortement l'ACP ! Pourquoi ?
  column_to_rownames("code_sp") 

# data_traits_for_PCA2 <- data_traits_for_PCA2%>% select(SLA,LNC, L_Area)  # /!\

for(i in 1:ncol(data_traits_for_PCA2)){
  data_traits_for_PCA2[is.na(data_traits_for_PCA2[,i]), i] <- 
    mean(data_traits_for_PCA2[,i], na.rm = TRUE)
}

pca_output <- perform_pca(data_traits_for_PCA2 )
pca_output[[7]]
coord_ind <- pca_output[[1]] %>% 
  merge(code_sp_lifeform,by="Code_Sp") 
  # full_join(data_clusters %>% rename(Code_Sp = code_sp))
coord_var <- pca_output[[2]]
var.explain.dim1 <- pca_output[[3]]
var.explain.dim2 <- pca_output[[4]]
var.explain.dim3 <- pca_output[[5]]
var.explain.dim4 <- pca_output[[6]]


Dim.A <- "Dim.1"
Dim.B <- "Dim.2"
Var.A <- var.explain.dim1
Var.B <- var.explain.dim2

# Faire des graphes de même taille
etendue_dim <- coord_ind %>% 
  gather(key=dim,value=coordinate,-c(Code_Sp,Lifelength)) %>% 
  group_by(dim) %>% 
  summarize(min = min(coordinate),max=max(coordinate)) %>% 
  filter(!(dim == "LifeForm1")) %>%
  mutate(min = as.numeric(min),max = as.numeric(max)) %>%
  mutate(etendue= max-min)

ratio <- etendue_dim %>% 
  select(dim,etendue) %>% 
  spread(key = dim, value= etendue) %>% 
  mutate(ratio12=Dim.1/Dim.2,
         ratio13 = Dim.1/Dim.3)

PCA_nat12 <-coord_ind %>% 
  # filter(!is.na(cluster)) %>% # ATTENTION juste pour simplifier, virer les NA
  plot_pca(coord_var,DimA= Dim.A ,DimB= Dim.B , Var.A,Var.B ) +
  ggtitle((expression(paste("GU"[S],sep='')))) +
  coord_fixed(ratio = ratio$ratio12) +
  # ellipse
  # ggforce::geom_mark_ellipse(data = coord_ind,aes(fill = cluster), # ,label = cluster
  #                            expand = unit(0.5,"mm"),
  #                            label.buffer = unit(-5, 'mm')) +
  scale_fill_brewer() 

PCA_nat13 <-coord_ind %>% 
  plot_pca(coord_var,DimA= Dim.A ,DimB= "Dim.3" , Var.A,var.explain.dim3 ) +
  ggtitle((expression(paste("GU"[S],sep='')))) +
  coord_fixed(ratio = ratio$ratio13)+
  # ellipse
  ggforce::geom_mark_ellipse(data = coord_ind,aes(fill = cluster), # ,label = cluster
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm')) +
  scale_fill_brewer() 

PCA_nat12
PCA_nat13

PCA_sp_names <-coord_ind %>% 
  ggplot (aes_string(x=Dim.A,y=Dim.B,label = "Code_Sp",colour="Lifelength"))+
  ggrepel::geom_label_repel() # With species names

plot_d1 <- boxplot_dimension(coord_ind,dim = "Dim.1")
plot_d2 <- boxplot_dimension(coord_ind,dim = "Dim.2")
PCA_nat_boxplot <- plot_pca_boxplot(PCA_nat12,plot_d1,plot_d2)

ggsave("outputs/figures/PCA_natif.png",PCA_nat_boxplot,height = 20, width =20)

coord_ind_nat <- coord_ind

# nb of life hisotries in each cluster
cont_table <- table(coord_ind$cluster,coord_ind$Lifelength)
chisq.test(cont_table)

## iii) Group PCA graphs ####
# PCA <- grid.arrange(PCA_fer12,PCA_nat12,PCA_fer13,PCA_nat13,
#                       layout_matrix=rbind(c(1,2),c(3,4)) )

PCA <- ggarrange(PCA_fer12,PCA_nat12,PCA_fer13,PCA_nat13,
          labels = c("A","B","C","D"))



# Extract the legend alone, from the data frame of species removal expe
leg <- ggpubr::get_legend(PCA_fer12)
legend <- ggpubr::as_ggplot(leg)

rapport <- 1
PCA12 <- ggarrange( PCA_fer12 +
                     theme(legend.position = c(0.9,0.9)) +
                     xlim(c(-4,7)) + 
                     ylim(c(-4,6))  +
                     coord_fixed(ratio = rapport) ,
                   PCA_nat12 +
                     theme(legend.position = "none") +
                     xlim(c(-4,7))+ 
                     ylim(c(-4,6)) +
                     coord_fixed(ratio = rapport),
                   # legend,
                 labels = c("A","B"),
                 height = 6,
                 row = 1)
PCA12

ggsave("draft/PCA_annuals_perennials.png",PCA12,height = 20, width =20)
ggsave("draft/PCA_annuals_perennials_legend.png",legend)


rbind(coord_ind_fer %>% mutate(treatment ="Fer"),
      coord_ind_nat %>% mutate(treatment = "Nat"))



plot <- MEAN_CSR_shallow_in_abundance %>% 
  filter(treatment%in% c("Nat","Fer")) %>% 
  mutate(zone = if_else(treatment == "Fer", "G+F","GU-S")) %>% 
  
  ggplot(aes_string(x="zone", y=trait, label = "code_sp",fill = "zone",shape = "LifeHistory")) +
  theme_classic()+
  theme(axis.title.x=element_blank())+
  # geom_boxplot(aes(color = LifeHistory)) +
  geom_boxplot() +
  geom_point(aes(color = LifeHistory,shape = LifeHistory),
             size = 2,
             position = position_dodge(width = .75)) +
  # scale_fill_manual(values = c("grey30", "grey80")) +
  scale_fill_manual(values = c("grey", "white")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ,
        axis.title.x = element_blank() )

## iv) ANOVA Position on axes ####
dimension <- 1
anov_dim <- compute_anova_dim_x(coord_ind,dimension)

par(mfrow=c(2,2)) ; plot(anov_dim) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(anov_dim))) # normality_graph

# Compare scores on axis "dimension" between annuals and perennials
info_anova <- info_anova_dim_x(anov_dim,dimension)
info_anova[[1]] # normality
info_anova[[2]] # homoscedasticity
info_anova[[3]] # independence of residuals
anov <- info_anova[[4]] # anova
sum <- info_anova[[5]] # summary
sum$coefficients
anov$`Pr(>F)`[1]

sum$adj.r.squared 


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
var_explained_pca_annual<-factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained
ggsave("draft/var_explained_pca_annuals.png",var_explained_pca_annual)

coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
var.explain.dim3 <- round(ACP1$eig[3,2])
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
  xlab(paste0("Dim1 (",var.explain.dim1,"%)"))+
  ylab(paste0("Dim2 (",var.explain.dim2,"%)")) +
  scale_colour_manual(values=c("#009E73","#E69F00"),
                      labels = c(expression(paste("G"^'+',"F",sep='')),
                                 expression(paste("GU"[S],sep=''))) )+
  scale_fill_manual('Program Type', values=c('pink','blue')) +
  labs(color = "Origin") +
  theme(plot.title = element_blank())


PCA2_dim3 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.3,colour=treatment))+ #,label = sp_tr
  theme_classic()+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.3*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.3*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste0("Dim1 (",var.explain.dim1,"%)"))+
  ylab(paste0("Dim3 (",var.explain.dim3,"%)")) +
  scale_colour_manual(values=c("#009E73","#E69F00"),
                      labels = c(expression(paste("G"^'+',"F",sep='')),
                                 expression(paste("GU"[S],sep=''))) )+
  scale_fill_manual('Program Type', values=c('pink','blue')) +
  labs(color = "Origin") +
  theme(plot.title = element_blank())
  
PCA_annuals <- gridExtra::grid.arrange(PCA2,PCA2_dim3, ncol=2,heights = 100)
ggsave("draft/PCA_annuals.jpg",PCA_annuals,width = 15, height = 11)

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
ggplot(data_traits_for_PCA,aes(x=LDMC,y=SLA,color = treatment)) +
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
# 4) CWM of annuals in both treatments ####
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

ggplot(MEAN,aes(x=log(SeedMass),y=LDMC,color = LifeHistory))+
  geom_point()

#_______________________________________________________
# Comparison of perennials ####
# PCA

traits_pca_perennials <- MEAN %>% 
  filter(LifeHistory == "perennial") %>% 
  select(!!c("code_sp","treatment",traits)) %>%  # subset of the traits that I want of analyse
  group_by(code_sp,treatment) %>% 
  summarise(across(all_of(traits), mean, na.rm= TRUE)) %>%  # /!\ JE POURRAIS NE PAS FAIRE DE MOYENNE, et prendre les données brutes !
  mutate(ddd = paste(code_sp,treatment,sep = "_")) %>% 
  column_to_rownames("ddd") %>% 
  select(-c(code_sp,treatment))


ACP1<-PCA(traits_pca_perennials,graph = FALSE)
var_explained_pca_perennials<-factoextra::fviz_eig(ACP1, addlabels = TRUE) # percentage of variance explained
var_explained_pca_perennials

coord_var <- data.frame(ACP1$var$coord) %>% 
  rownames_to_column()
var.explain.dim1 <- round(ACP1$eig[1,2])
var.explain.dim2 <- round(ACP1$eig[2,2])
var.explain.dim3 <- round(ACP1$eig[3,2])
coord_ind <- data.frame(ACP1$ind$coord) %>% 
  rownames_to_column("sp_tr") %>% 
  mutate(code_sp = str_sub(sp_tr,1L,-5L)) %>% 
  mutate(treatment = str_sub(sp_tr,-3L,-1L))

PCA2 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.2,colour=treatment)) + #,label = sp_tr
  theme_classic() +
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.2*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.2*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste0("Dim1 (",var.explain.dim1,"%)"))+
  ylab(paste0("Dim2 (",var.explain.dim2,"%)")) +
  scale_colour_manual(values=c("#009E73","#E69F00"),
                      labels = c(expression(paste("G"^'+',"F",sep='')),
                                 expression(paste("GU"[S],sep=''))) )+
  scale_fill_manual('Program Type', values=c('pink','blue')) +
  labs(color = "Origin") +
  theme(plot.title = element_blank())


PCA2_dim3 <- ggplot(coord_ind,aes(x=Dim.1,y=Dim.3,colour=treatment))+ #,label = sp_tr
  theme_classic()+
  geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
  geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
  coord_equal() +
  geom_point() +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1*7-0.2, yend=Dim.3*7-0.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  geom_text_repel(data=coord_var, aes(x=Dim.1*7, Dim.3*7, label=rowname), size = 4, vjust=1, color="black") +
  ggtitle("Annuals in fer and nat sup") +
  # theme(legend.position = "none")+
  xlab(paste0("Dim1 (",var.explain.dim1,"%)"))+
  ylab(paste0("Dim3 (",var.explain.dim3,"%)")) +
  scale_colour_manual(values=c("#009E73","#E69F00"),
                      labels = c(expression(paste("G"^'+',"F",sep='')),
                                 expression(paste("GU"[S],sep=''))) )+
  scale_fill_manual('Program Type', values=c('pink','blue')) +
  labs(color = "Origin") +
  theme(plot.title = element_blank())

PCA_annuals <- gridExtra::grid.arrange(PCA2,PCA2_dim3, ncol=2,heights = 100)
ggsave("draft/PCA_annuals.jpg",PCA_annuals,width = 15, height = 11)

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





# Clustering tests ####
library(cluster)
ftreatment <- "Nat"


  
data_clustering <- data_traits_for_PCA %>% 
  filter(treatment == ftreatment) %>% 
  # mutate(sp_trtmt = paste(code_sp,treatment,sep="_")) %>% 
  column_to_rownames("code_sp") %>% 
  select(-c(treatment)) 
  # select(!!ftrait)
  # select(SLA,LDMC,L_Area,LCC,LNC)
# select(c(Ldelta13C,Hrepro,Dmax))

# clustering sur dimensions d'ACP (et ensuite le faire sur traits bruts en annexe)
# data_clustering <- coord_ind %>% 
#   column_to_rownames("Code_Sp") %>% 
#   select(-c(Lifelength,LifeForm1,cluster))

# impute values with mice
# https://datascienceplus.com/handling-missing-data-with-mice-package-a-simple-approach/
library(mice)
init = mice(data_clustering, maxit=0) 
meth = init$method
predM = init$predictorMatrix

# predM[, c("BMI")]=0 # I select the BMI variable to not be included as predictor during imputation
# meth[c("Age")]=""# kip a variable from imputation (e.g. Ager). This variable will be used for prediction.

# meth[(c("SeedMass"))] = "norm"
# I set different methods for each variable. You can add more than one variable in each method.
# "norm" ; "logreg" ; "polyreg"
# set.seed(103)
imputed = mice(data_clustering, method=meth, predictorMatrix=predM, m=5)
imputed <- complete(imputed)

imputed_PCA <- data_clustering %>% 
  mutate_at(vars(colnames(data_clustering)), ~replace_na(.,mean(., na.rm = TRUE)))


# clustering
data_clustering2 <- imputed %>% 
  na.omit() %>%
  scale()

#WITH PCA 
data_clustering2 <- imputed_PCA %>% 
  scale()

# find the optiomal number of clusters (2 examples)
fviz_nbclust(data_clustering2, kmeans, method = "wss")

gap_stat <- clusGap(data_clustering2,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
fviz_gap_stat(gap_stat)


#make this example reproducible
#perform k-means clustering with k = 2 clusters
km <- kmeans(data_clustering2, centers = 2, nstart = 25)

#view results
data_clusters <- km$cluster%>% 
  as.data.frame() %>% 
  rownames_to_column("code_sp")
colnames(data_clusters) <- c("code_sp", "cluster")
data_clusters <- data_clusters%>% 
  mutate(cluster = as.factor(cluster))

fviz_cluster(km, data = data_clustering2)




# nb of life hisotries in each cluster
clust_res <- imputed %>% 
  na.omit() %>% 
  rownames_to_column("Code_Sp") %>% 
  merge(data_clusters %>% rename(Code_Sp = code_sp) ) %>% 
  merge(code_sp_lifeform) %>% 
  mutate(Lifehistory = case_when(LifeForm1 == "The" ~ "annual",
                                 TRUE ~"perennial")) 


cont_table <- table(clust_res$cluster,clust_res$Lifehistory)
# chisq.test(cont_table)
cont_table

# Fisher’s exact test is an alternative to chi-squared test used mainly when a chi-square 
# approximation is not satisfactory (i.e., small sample size, and you get the warning message) 
ftest <- fisher.test(cont_table)
pval <- ftest$p.value
ftest


cluster_table_nat  <- cont_table %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster") %>% 
  mutate(p.val = c(round(pval,digits =3),""))

cluster_table <- merge(cluster_table_fer,cluster_table_nat,by="cluster")

TABLE <- cluster_table %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Cluster", "Nb of annuals", "Nb of perennials","p.val",
                                   "Nb of annuals", "Nb of perennials","p.val")) %>%
  kableExtra::kable_styling("hover", full_width = F) %>% 
  kableExtra::add_header_above(c(" "=1,"G+F" = 3,"GUS" = 3))

cat(TABLE, file = "draft/clustering.doc")


