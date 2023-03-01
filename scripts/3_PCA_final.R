library(tidyverse)
library(FactoMineR)
library(ggpubr)
library(rstatix)
library(cowplot)
library(funrar)
library(mFD)


MEAN_multivar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv") %>% 
  select(-Disp) %>% 
  dplyr::rename(LCCm = LCC) %>% 
  dplyr::rename(LNCm = LNC)
# MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")

code_sp_lifeform <- read.csv2("data/species_names_lifehistory.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCCm","LNCm","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

fMEAN <- MEAN_multivar %>%
  select(species,code_sp,LifeHistory,treatment,any_of(traits)) %>% 
  na.omit() %>% # WE DON'T FILL THE MATRIX
  mutate(sp_trt = paste(code_sp,treatment,sep="_"))
rownames(fMEAN) <- NULL

# fMEAN %>% 
#   group_by(LifeHistory,treatment) %>% 
#   summarize(n = n())


data_hypervolume <- fMEAN %>% 
  # filter(treatment == "Nat") %>%
  # filter(LifeHistory == "annual") %>%
  column_to_rownames("sp_trt") %>% 
  select(-c(species,code_sp,LifeHistory,treatment))

# indices <- data_hypervolume %>%
#   dbFD()
# indices$FRic

data_hypervolume_Nmass <- data_hypervolume %>% 
  mutate(Nmass = LNC * 10) %>% # change unit from mg/g to %
  mutate(LMA = 1/SLA * 10^7) %>% # convert kg/cm² to g/m² (1kg sur 1cm² = 1^4 kg sur 1m², = 10^7 g sur 1m²)
  mutate(log10_Amass = 0.74 * log10(Nmass) - 0.54 * log10(LMA) + 2.96) %>% #equation from supp. data Wright et al. 2003
  mutate(Amass = 10^log10_Amass) %>% 
  select(-c(Nmass,LMA,log10_Amass))

# PCA_hypervolume <- PCA(data_hypervolume_Nmass,scale.unit=TRUE,graph=T,quanti.sup = 9) # for A maxx
PCA_hypervolume <- PCA(data_hypervolume,scale.unit=TRUE,graph=F) 
percent_var <- factoextra::fviz_eig(PCA_hypervolume, addlabels = TRUE, ylim = c(0, 30))
# point d'inflexion sur le troisième axe. Présenter les trois axes
var.explain.dim1 <- round(PCA_hypervolume$eig[1,2])
var.explain.dim2 <- round(PCA_hypervolume$eig[2,2])
var.explain.dim3 <- round(PCA_hypervolume$eig[3,2])


# Distinctiveness ####
traits_dist<-function(traits){
  # Computes functional distinctiveness on all the species and adds a new columns with it on the input data frame (i.e. traits, here).
  dist_tot<-compute_dist_matrix(traits,metric='euclidean')
  tidy<-as.data.frame(rownames(traits)) # tidy format for computing distinctiveness in the fonction below
  colnames(tidy)<-"SName"
  distinct_tot<-distinctiveness_com(com_df=tidy,
                                    sp_col=colnames(tidy),abund=NULL,
                                    dist_matrix=dist_tot,relative=F)
  traits$Di<-distinct_tot$Di
  traits
}

# avec les deux premières dim de l'ACP
dist <- PCA_hypervolume$ind$coord %>% 
  as.data.frame() %>%
  select(Dim.1    ,   Dim.2) %>%
  traits_dist() %>% 
  rownames_to_column("sp_tr") %>% 
  separate(sp_tr, into = c("code_sp","treatment"),sep = "_") %>% 
  merge(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","Annuals","Perennials")) %>% 
  mutate(treatment = if_else(treatment == "Fer", "Intensive","Extensive"))

plot_di <- dist  %>% 
  ggplot(aes(x = LifeHistory,y=Di)) +
  geom_boxplot() +
  facet_wrap(~treatment) +
  geom_point() +
  ggrepel::geom_text_repel(data = dist %>% filter(code_sp %in% c("CARDNUTA","RUMEACET","STIPPENN","TEUCMONT","FESTCHRI","BELLPERE")),
                            aes(x = LifeHistory,y=Di,label = species)) +
  theme_classic() +
  ylab("Functional distinctiveness") +
  xlab("")

ggsave("draft/sp_distinctes.png",plot_di)


# Functional diversity ####

# abundance
# ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
#   rename(Ligne = id_transect_quadrat) %>% 
#   select(species,Ligne,relat_ab)
# ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
#   filter(depth == "S") %>% 
#   select(species,Ligne,relat_ab)
# abundance_studied <- rbind(ab_fer,ab_nat) %>% 
#   unique()
# 
# ab_sp <- abundance_studied %>%
#   spread(species,relat_ab) %>%
#   column_to_rownames("Ligne") %>%
#   replace(is.na(.),0)
# # Prendre l'intersection (espèces en commun) des data.frames des traits des abondances
# sp_common <- intersect( colnames(ab_sp) , fMEAN$species )
# 
# data_abundance <- ab_sp[,sp_common] %>%
#   as.matrix()

# 
# fspace <- tr.cont.fspace(
#   sp_tr = data_hypervolume,
#   pca = TRUE,
#   nb_dim = 3, # Nombre d'axes que l'on souhaite conserver.
#   scaling = "scale_center",
#   compute_corr = "pearson")
# 
# data_abundance <- matrix(nrow = 2,ncol = dim(data_hypervolume[2]),
#                          dimnames = list(NULL,rownames(data_hypervolume)))
# data_abundance[1,] <- rep(1, times = dim(data_hypervolume)[1] )
# data_abundance[2,] <- rep(1, times = dim(data_hypervolume)[1] )
# 
# alpha_fd <- alpha.fd.multidim(fspace$sp_faxes_coord ,
#                               asb_sp_w = data_abundance ,
#                               ind_vect = c("fric"))
# indices <- alpha_fd$functional_diversity_indices
# head(indices)


# Number of dimensions ####
# source("scripts/functions/dimensionality_script_mouillot_V1.R")
# AUC <- dimension_funct_taina(trait_df = data.frame(scale(data_hypervolume, center=T, scale=T)), dim_pcoa = 10, rep = 99, cores = 1,
#                       metric_scaled = TRUE)
# # three dimensions selected!
# ggplot(AUC, aes(x = dim, y = AUC))+
#   geom_point() + 
#   geom_smooth()


coord_ind <- PCA_hypervolume$ind$coord %>% 
  as.data.frame() %>% 
  rownames_to_column("sp_trt") %>% 
  separate("sp_trt",into = c("code_sp","treatment"),sep="_") %>% 
  left_join(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

coord_var <- PCA_hypervolume$var$coord %>% 
  as.data.frame() %>% 
  rownames_to_column("trait")

# identify species on both treatment, or in one treatment only
zones <- coord_ind %>% 
  select(code_sp,treatment,LifeHistory) %>% 
  mutate(presence = 1) %>% 
  spread(key = treatment,value = presence) %>% 
  mutate(zone = case_when(is.na(Fer) & Nat ==1 ~ "Nat",
                          Fer == 1 & Nat == 1 ~ "Both",
                          TRUE ~ "Fer")) %>% 
  select(code_sp,zone,LifeHistory) %>% 
  mutate(zone2 = if_else(zone == "Both","Both","One"))
coord_ind <- full_join(coord_ind,zones)

fLH <- "annual"

DimA <- "Dim.1"
DimB <- "Dim.2"
var.explain.dimA <- var.explain.dim1
var.explain.dimB <- var.explain.dim2



PLOT <- NULL
i <- 0
for (fLH in c("annual","perennial")){
  i <- i+1
  plot_ACP <- coord_ind %>% 
    filter(LifeHistory == fLH) %>%
    ggplot(aes_string(x=DimA,y=DimB), width = 10, height = 10)+
    theme_classic() +
    geom_hline(aes(yintercept=0), size=.2,linetype="longdash") + 
    geom_vline(aes(xintercept = 0),linetype = "longdash", size=.2)+
    # coord_equal() +
    geom_point(size=2,aes(color = treatment,shape = zone2)) + # ,colour=treatment
    # geom_segment(data=coord_var, aes(x=0, y=0, xend=get(DimA)*5.5-0.2, yend=get(DimB)*5.5-0.2), 
    #              arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    # ggrepel::geom_text_repel(data=coord_var, aes(x=get(DimA)*5.5, get(DimB)*5.5, label=trait), size = 6, vjust=1, color="black") +
    theme(legend.position = "none") +
    xlab("Dim.1") + ylab("Dim.2")+
    theme(text=element_text(size=10)) +
    ggforce::geom_mark_ellipse(data = coord_ind %>% 
                                 filter(LifeHistory == fLH),aes(fill = treatment), # ,label = cluster
                               expand = unit(0.5,"mm"),
                               label.buffer = unit(-5, 'mm')) +
    # {if(fLH == "annual") ggtitle("Annuals") } +
    # {if(fLH == "perennial") ggtitle("Perennials") } +
    xlim(c(-4,7.5)) +
    ylim(c(-4,6)) +
    {if(fLH == "annual") scale_shape_manual(values = c(1,19)) } + # pour les annuelles 
    {if(fLH == "perennial") scale_shape_manual(values = c(2,17)) } # pour les pérennes METTRE UN TRIANGLE A L'ENVERS POUR PERENNES MESUREES DANS LES DEUX
  
  if (fLH == "annual"){
    data_sp_names <- coord_ind %>% 
      filter(LifeHistory == fLH) %>%
      filter(code_sp %in% c("BUPLBALD","MINUHYBR","HORNPETR",
                            "VICISATI-SATI","CREPVESI-HAE","HORDMURI","TRIFSCAB","BROMDIAN"))
  } else{
    data_sp_names <- coord_ind %>% 
      filter(LifeHistory == fLH) %>%
      filter(code_sp %in% c("STIPPENN","TEUCMONT","FESTCHRI",
                            "CARDNUTA","RUMEACET","TRIFREPE","POATRIV","BELLPERE"))
  }

  
  plot_ACP2 <- plot_ACP +
    ggrepel::geom_text_repel(data = data_sp_names,
                             aes(x = Dim.1, y = Dim.2, label = code_sp),
                             max.overlaps = 10000)
                             # min.segment.length = 0) # ,box.padding = 2
  plot_ACP2
  
  PLOT[[i]] <- plot_ACP
  

  # https://www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets/
  
  # https://ggrepel.slowkow.com/articles/examples.html
  
  to_boxplot <- coord_ind %>% 
    filter(LifeHistory == fLH) %>% 
    mutate(Management = if_else(treatment == "Fer","Int.","Ext.")) %>% 
    mutate(Management = factor(Management, levels = c("Int.","Ext."))) %>% 
    # select(-treatment) %>% 
    gather(key = dim, value = coordinate, - c(code_sp,Management,treatment,LifeHistory,species,LifeForm1,zone2)) %>% 
    filter(dim %in% c("Dim.1","Dim.2")) %>% 
    mutate(Coordinate = as.numeric(coordinate))
  

  
  stat.test <- to_boxplot %>%
    group_by(dim) %>%
    t_test(Coordinate ~ Management) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance() %>% 
    add_xy_position(x = "Management")
  
  bxp <- ggboxplot(
    to_boxplot, x = "Management", y = "Coordinate", fill = "#00AFBB", 
    facet.by = "dim"
  )
  
  bxp +stat_pvalue_manual(stat.test)
  
  boxplot <- to_boxplot %>% 
    ggplot(aes(x = Management, y = Coordinate,color = Management)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~dim) +
    theme_classic() +
    # ggtitle(fLH) +
    ylab("Coordinate") +
    theme(legend.position  = "none") +
    # ggsignif::geom_signif(
    #   data = annotation_df,
    #   aes(xmin = start, xmax = end, annotations = label, y_position = y),
    #   textsize = 3, vjust = -0.2,
    #   manual = TRUE
    # )
    
    # MANUAL ANNOTATION ggsignif!!
    # ggsignif::geom_signif(data = data.frame(Group = c("Dim.1","Dim.2")),
    #                         aes(y_position=c(5.3, 6.3), xmin=c(0.8, 0.8), xmax=c(1.2, 1.2),
    #                             annotations=c("**", "NS")), tip_length=0, manual = T)
    # ggsignif::geom_signif(comparisons = list(c("annual", "perennial","Extensive","Intensive")), 
    #                       map_signif_level=TRUE) +
    geom_line(aes(group = code_sp),
              alpha = 0.4, color = "black") +
    {if(fLH == "annual") geom_point(aes(shape = zone2),size = 2) } +
    {if(fLH == "perennial") geom_point(aes(shape = zone2),size = 2) } + # position = position_dodge(width = 0.75)) SSi pas shape = zone2
    {if(fLH == "annual") scale_shape_manual(values = c(1,19)) } + # pour les annuelles 
    {if(fLH == "perennial") scale_shape_manual(values = c(2,17)) } + # pour les pérennes METTRE UN TRIANGLE A L'ENVERS POUR PERENNES MESUREES DANS LES DEUX
    ylim(c(-3,8)) +
    stat_pvalue_manual(stat.test)
    

  
  i <- i+1
  PLOT[[i]] <- boxplot
  
  
}

## Plot variables ####
coord_axes <- PCA_hypervolume$var$coord %>%
  as.data.frame() %>% 
  rownames_to_column("trait") %>% 
  mutate(trait = case_when(trait == "L_Area" ~ "LA",
                           trait == "Hrepro" ~ "H",
                           TRUE ~ trait))

plot_axis_pca <- ggplot(coord_axes) +
  geom_segment( aes(x=0, y=0, xend=Dim.1, yend=Dim.2), 
                arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  ggrepel::geom_label_repel( aes(x=Dim.1, Dim.2, label=trait), size = 4, vjust=1, color="black")  +
  theme_classic() +
  xlab(paste0(DimA," (",var.explain.dimA,"%)"))+
  ylab(paste0(DimB," (",var.explain.dimB,"%)")) 


## plot legend ####
plot <- coord_ind %>% 
  rename(Management = treatment) %>% 
  mutate(Management = if_else(Management == "Fer","Intensive","Extensive")) %>% 
  filter(LifeHistory == "annual") %>%
  ggplot(aes_string(x=DimA,y=DimB), width = 10, height = 10)+
  theme_classic() +
  geom_point(size=4,aes(color = Management)) + # ,shape = zone2# ,colour=treatment
  theme(text=element_text(size=15)) +
  ggforce::geom_mark_ellipse(data = coord_ind %>%
                               rename(Management = treatment) %>% 
                               mutate(Management = if_else(Management == "Fer","Intensive","Extensive")) %>%
                               filter(LifeHistory == fLH),aes(fill = Management),
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm')) +
  labs(color = 'Management',fill='Management') +
  scale_color_manual(labels = c("Extensive", "Intensive"), values = c("#00BFC4","#F8766D"))+
  scale_fill_manual(labels = c("Extensive", "Intensive"), values = c("#00BFC4","#F8766D"))
  

# CHANGE LEGEND NAMES

leg <- ggpubr::get_legend(plot)
legend <- ggpubr::as_ggplot(leg)

# 
# PCA <- ggpubr::ggarrange(PLOT[[1]],PLOT[[3]],plot_var,
#                          PLOT[[2]],PLOT[[4]],legend,ncol = 3,nrow = 2)

# PCA <- gridExtra::grid.arrange(plot_axis_pca, legend,PLOT[[1]],PLOT[[3]],
#                          PLOT[[2]],PLOT[[4]],ncol = 2,nrow = 3)
# 

PCA <- plot_grid(plot_axis_pca, legend,
          PLOT[[1]],PLOT[[3]],
          PLOT[[2]],PLOT[[4]],
          ncol = 2, align = "hv",
          labels = c("A","","B","C","D","E"))

ggsave("draft/PCA.png",PCA,height = 10, width =6)


# trash
plot1 <- ggplot(coord_var,aes(x=Dim.1,y=Dim.2,label=trait)) +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1, yend=Dim.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  ggrepel::geom_text_repel(data=coord_var, aes(x=Dim.1, Dim.2, label=trait), size = 6, vjust=1, color="black") +
  # geom_point(data = coord_ind)+
  theme_classic()

plot2 <- ggplot(coord_var,aes(x=Dim.1,y=Dim.3,label=trait)) +
  geom_segment(data=coord_var, aes(x=0, y=0, xend=Dim.1, yend=Dim.3), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  ggrepel::geom_text_repel(data=coord_var, aes(x=Dim.1, Dim.3, label=trait), size = 6, vjust=1, color="black") +
  theme_classic()

pcas <- gridExtra::grid.arrange(plot1, plot2, ncol=2)





# ggsave("outputs/plots/PCA_fer.png",pcas,height = 10, width = 14)

# correlation

panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = rgb(0, 1, 1, alpha = 0.5), ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

plot_cor <- fMEAN %>% 
  # filter(treatment == "Fer") %>%
  column_to_rownames("sp_trt") %>% 
  select(-c(code_sp,LifeHistory,treatment)) %>% 
  mutate(L_Area = log(L_Area)) %>% 
  pairs(diag.panel = panel.hist,
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines)
plot_cor



SLA_LDMC <- function(x) 11.3 * 10^4 * x ^(-1.58)
# (Equation de Garnier, 2001, functional ecology)


fMEAN %>% 
  # filter(treatment == "Fer") %>% 
  # filter(LifeHistory == "annual") %>% 
  ggplot(aes(x = LDMC, y = SLA,color = LifeHistory)) +
  geom_point() +
  geom_function(fun = SLA_LDMC,color = "black")


fMEAN %>% 
  filter(treatment == "Nat") %>%
  # filter(LifeHistory == "annual") %>%
  ggplot(aes(x = log(L_Area), y = LNC)) +
  geom_point(aes(color = LifeHistory)) 
  # geom_smooth(method = "lm")

fMEAN %>% 
  # filter(treatment == "Fer") %>%
  # filter(LifeHistory == "annual") %>%
  ggplot(aes(x = SLA, y = LNC)) +
  geom_point(aes(color = LifeHistory))

