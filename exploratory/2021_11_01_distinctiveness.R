source("scripts/1. Packages.R")
source("scripts/2. Import files.R")

ABUNDANCE_traits <- read.csv2("outputs/data/pooled_abundance_and_traits.csv")
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

# On SLA only
library(funrar)

maud <- ABUNDANCE_traits %>%
  filter(dataset=="Maud") %>% 
  group_by(depth) %>% 
  mutate(relat_abundance = abundance/sum(abundance))

deep_maud <- maud %>% 
  filter(depth=="D")

superficial_maud <- maud %>% 
  filter(depth=="S")

traits_table <- MEAN %>% 
  filter(Trtmt == "Nat") %>% 
  unique() %>% 
  column_to_rownames("Code_Sp") 
  # select(-c(Trtmt,Species,Samily,LifeForm1,LifeForm2))
  
# /!\ refaire le tableau MEAN : j'ai des espèces dupliquées, e.g. parce qu'à un moment, ce n'est pas la bonne famille qui est rentrée.
# En attendant

traits_table <- LeafMorpho %>% 
  filter(Treatment %in% c("Nat_Dlm","Nat_Dol","Nat_Sab")) %>% 
  select(Code_Sp,SLA) %>% 
  filter(Code_Sp %in% maud$code_sp) %>% 
  group_by(Code_Sp) %>% 
  mutate(SLA=mean(SLA)) %>% 
  unique() %>% 
  column_to_rownames("Code_Sp")
  
distance <- compute_dist_matrix(traits_table = traits_table, metric = "euclidean") %>% 
  as.matrix()
  

maud2 <- maud %>% 
  filter(code_sp %in% rownames(distance))

superficial_maud2 <- superficial_maud %>% 
  filter(code_sp %in% rownames(distance))

distin <- distinctiveness_stack(com_df = superficial_maud2,
                      sp_col = "code_sp", com = "depth", abund = "relat_abundance",
                      dist_matrix = distance,
                      relative = F) # a logical indicating if distinctiveness should be scaled relatively to the community (scaled by max functional distance among the species of the targeted community)


ggplot(distin,aes(x=LifeForm1,y=Di))+
  geom_boxplot()
