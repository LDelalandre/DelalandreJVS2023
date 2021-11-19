source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
source("scripts/functions/traitRange.test.R")
source("scripts/functions/traitMoments.test.R")

library(funrar)

ABUNDANCE_traits <- read.csv2("outputs/data/pooled_abundance_and_traits.csv")
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
# Maud <- ABUNDANCE_traits %>% filter(dataset == "Maud")

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8")) %>% 
  mutate(plot = paste(depth,paddock,sep="_"))

# Maud's data ####
# Preparing data ####

# stacked abundance
ab_Maud_unclean <- ABUNDANCE %>% 
  filter(dataset=="Maud") %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  filter(!(is.na(LifeHistory))) %>% 
  # NB : /!\ I sould keep family and lifeform info, but only when these are variables are clean in the dataset.
  # (sinon, ça découple une espèce en deux artificiellement dans le calcul de la moyenne).
  mutate(plot=paste(depth,paddock,sep="_"))



ab_Maud <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", sheet = "abondances par parcelle", 
                   startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(Species,abundance,-plot) %>%
  mutate(Species=str_replace(Species,"_"," ")) %>% 
  full_join(MEAN,by="Species") %>% 
  dplyr::rename(species=Species,code_sp=Code_Sp) %>% 
  filter(!is.na(abundance)) %>% 
  mutate(depth = str_sub(plot,start = 1L,end=1L)) %>% 
  mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
  select(-plot) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth"))



# traits
traits <- MEAN %>% 
  filter(Trtmt == "Nat" & Species %in% ab_Maud$species) %>% 
  # select(code_sp,Nb_Lf:Mat) %>%
  select(Code_Sp,Nb_Lf:Mat ) %>%
  unique() %>% 
  arrange(Code_Sp) %>% 
  column_to_rownames("Code_Sp") %>% 
  filter(!is.na(SLA  ))


# Distinctiveness ####
distance <- compute_dist_matrix(traits_table = traits, metric = "euclidean") 

# Species in line and plots in column
# and summarize at the scale desired
df_comm <- ab_Maud %>% 
  group_by(code_sp,plot) %>% 
  select(code_sp,plot,abundance) %>%
  summarize(abundance = sum(abundance)) %>% 
  ungroup() %>% 
  group_by(plot) %>% 
  mutate(rel_ab = abundance/sum(abundance)) %>% 
  select(code_sp,plot,rel_ab) %>% 
  filter(code_sp %in% colnames(distance)) %>% # select the species for which we have the traits
  spread(code_sp,rel_ab) %>%
  column_to_rownames("plot")

comm <- df_comm %>% 
  replace(is.na(.),0) %>%
  as.matrix() 

# comm[comm>0] <-1 # Presence-absence only, instead of abundance

distin <- distinctiveness(pres_matrix = comm,
                          relative=F,
                          dist_matrix = distance)


distin_stack <- distin %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("code_sp") %>% 
  gather(key = plot,value = Di,-code_sp) %>% 
  filter(!is.na(Di)) 

Maud_di <- ab_Maud %>% 
  select(PC1score,depth,paddock,species,code_sp,plot,LifeHistory) %>% 
  filter(!is.na(LifeHistory)) %>% 
  full_join(distin_stack,.,by=c("code_sp","plot"))


ggplot(Maud_di,aes(x=PC1score,y=Di,color =LifeHistory,label=code_sp))+
  geom_point()

ggplot(Maud_di,aes(x=LifeHistory,y=Di))+
  geom_boxplot() +
  facet_wrap(~plot)

# Distinctiveness without abundance
tidy<-as.data.frame(rownames(traits)) # tidy format for computing distinctiveness in the fonction below
colnames(tidy)<-"code_sp"
distinct_tot<-distinctiveness_com(com_df=tidy,
                                  sp_col=colnames(tidy),abund=NULL,
                                  dist_matrix=distance,relative=F) %>% 
  merge(Maud_di %>% select(code_sp,LifeHistory),by=c("code_sp")) %>% 
  unique()

ggplot(distinct_tot,aes(x=LifeHistory,y=Di))+
  geom_point()

# Reagrder une communauté en particulier
Maud <- Maud %>% 
  filter(community == "S_P10")

com <- Maud %>% 
  select(code_sp,community,rel_ab)

traits <- Maud %>% 
  select(code_sp,Nb_Lf:Mat) %>%
  # select(code_sp,SLA) %>%
  unique() %>% 
  column_to_rownames("code_sp")

distance <- compute_dist_matrix(traits_table = traits, metric = "euclidean")
distin <- distinctiveness_stack(com_df = com,
                                sp_col = c("code_sp"), 
                                com = c("community"), 
                                abund = c("rel_ab"),
                                dist_matrix = distance)

Maud_di <- full_join(Maud,distin,by=c("code_sp","community","rel_ab"))

ggplot(Maud_di,aes(x=LifeHistory,y=Di))+
  geom_boxplot() +
  facet_wrap(~community)


# Trait range ####
comm <- as.data.frame(comm)

null_range <- null.range(comm = comm,trait = "SLA", df = traits, nreps = 9999) %>% 
  mutate(ES = 2*(P.lower.range-0.5)) %>% 
  mutate(plot=rownames(comm)) %>% 
  full_join(soil_Maud,by="plot")

ggplot(null_range,aes(x=PC1score,y=ES))+
  geom_point() +
  geom_hline(yintercept=0) +
  geom_segment(aes(x = PC1score, y = 0, xend = PC1score, yend = ES))+
  ggtitle(trait)
