library(tidyverse)
library(FactoMineR)


MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")

code_sp_lifeform <- read.csv2("data/species_names_lifehistory.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

fMEAN <- MEAN %>% 
  select(code_sp,LifeHistory,treatment,all_of(traits)) %>% 
  select(-Disp) %>% 
  na.omit() %>% # WE DON'T FILL THE MATRIX
  mutate(sp_trt = paste(code_sp,treatment,sep="_"))
rownames(fMEAN) <- NULL


trt <- "Nat"
data_hypervolume <- fMEAN %>% 
  # filter(treatment == trt) %>%
  # filter(LifeHistory == "annual") %>% 
  column_to_rownames("sp_trt") %>% 
  select(-c(code_sp,LifeHistory,treatment)) 


PCA_hypervolume <- PCA(data_hypervolume,scale.unit=TRUE,graph=F)
percent_var <- factoextra::fviz_eig(PCA_hypervolume, addlabels = TRUE, ylim = c(0, 30))
# point d'inflexion sur le troisième axe. Présenter les trois axes
var.explain.dim1 <- round(PCA_hypervolume$eig[1,2])
var.explain.dim2 <- round(PCA_hypervolume$eig[2,2])
var.explain.dim3 <- round(PCA_hypervolume$eig[3,2])

# Number of dimensions ####
source("scripts/functions/dimensionality_script_mouillot_V1.R")
AUC <- dimension_funct_taina(trait_df = data.frame(scale(data_hypervolume, center=T, scale=T)), dim_pcoa = 10, rep = 99, cores = 1,
                      metric_scaled = TRUE)
# three dimensions selected!
ggplot(AUC, aes(x = dim, y = AUC))+
  geom_point() + 
  geom_smooth()



coord_ind <- PCA_hypervolume$ind$coord %>% 
  as.data.frame() %>% 
  rownames_to_column("sp_trt") %>% 
  separate("sp_trt",into = c("code_sp","treatment"),sep="_") %>% 
  left_join(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

coord_var <- PCA_hypervolume$var$coord %>% 
  as.data.frame() %>% 
  rownames_to_column("trait")

fLH <- "perennial"

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
    coord_equal() +
    geom_point(size=4,aes(shape = treatment)) + # ,colour=treatment
    geom_segment(data=coord_var, aes(x=0, y=0, xend=get(DimA)*5.5-0.2, yend=get(DimB)*5.5-0.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
    geom_text_repel(data=coord_var, aes(x=get(DimA)*5.5, get(DimB)*5.5, label=trait), size = 6, vjust=1, color="black") +
    # theme(legend.position = "none") +
    xlab(paste0(DimA," (",var.explain.dimA,"%)"))+
    ylab(paste0(DimB," (",var.explain.dimB,"%)")) + 
    theme(text=element_text(size=20)) +
    theme(legend.position = "none")+
    {if(fLH == "annual") scale_shape_manual(values = c(1,19)) } + # pour les annuelles
    {if(fLH == "perennial") scale_shape_manual(values = c(2,17)) } # pour les pérennes
  
  PLOT[[i]] <- plot_ACP
  
  boxplot <- coord_ind %>% 
    filter(LifeHistory == fLH) %>% 
    gather(key = dim, value = coordinate, - c(code_sp,treatment,LifeHistory,species,LifeForm1)) %>% 
    filter(dim %in% c("Dim.1","Dim.2")) %>% 
    ggplot(aes(x = treatment, y = coordinate,shape = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~dim) +
    theme_classic() +
    geom_point(position = position_dodge(width = 0.75)) +
    ggtitle(fLH) +
    theme(legend.position = "none")+
    {if(fLH == "annual") scale_shape_manual(values = c(1,19)) } + # pour les annuelles
    {if(fLH == "perennial") scale_shape_manual(values = c(2,17)) } # pour les pérennes
  
  i <- i+1
  PLOT[[i]] <- boxplot
  
  
}



PCA <- ggarrange(PLOT[[1]],PLOT[[3]],PLOT[[2]],PLOT[[4]],
                 labels = c("A","B","C","D"))

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
