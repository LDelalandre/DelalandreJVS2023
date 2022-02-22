source("scripts/1. Packages.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
MEAN_pheno <- MEAN %>% 
  select(species,treatment,code_sp,LifeHistory,LifeForm1,Flo,Disp,Mat,Mat_Per)

#_________________________________________________________________________
# Annuals in the natif
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")
annuals_nat <- ab_nat %>% 
  filter(LifeHistory=="annual") %>% 
  pull(code_sp) %>% 
  unique()

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
annuals_fer <- ab_fer %>% 
  filter(LifeForm1=="The") %>% 
  pull(code_sp) %>% 
  unique()

# Mat and Disp are similar. Use Disp (?) = day of seed dispersal.
# NB: Mat est un trait issu de la feuille Seed, et pas Pheno
MEAN_pheno %>% 
  ggplot(aes(x=Mat,y=Disp))+
  geom_point()


# Difference of phenology of annuals between Fer and Nat
#  /!\ C'est peut-être un effet du moment de la mesure d'Adeline : elle a peut-être noté 
# les espèces en fleur ou qui dispersaient au moment où elle était là !!
MEAN_pheno %>% 
  filter(LifeHistory == "annual") %>%
  ggplot(aes(x=treatment,y= Mat_Per))+
  geom_boxplot()

annuals_nat_pheno <- MEAN_pheno %>% 
  filter(LifeHistory=="annual") %>% 
  filter(treatment == "Nat") %>% 
  filter(!is.na(Disp)) %>% 
  pull(code_sp) %>% 
  unique()

annuals_fer_pheno <- MEAN_pheno %>% 
  filter(LifeHistory=="annual") %>% 
  filter(treatment == "Fer") %>% 
  filter(!is.na(Disp)) %>% 
  pull(code_sp) %>% 
  unique()


intersect(annuals_nat_pheno,annuals_nat)
intersect(annuals_fer_pheno,annuals_fer)

# annuals in relevés, but pheno not measured
setdiff(annuals_nat,annuals_nat_pheno)
pheno_lacking_fer <- setdiff(annuals_fer,annuals_fer_pheno)

sp_leo <- pheno_leo %>% 
  pull(code_sp) %>% 
  unique() %>% 
  as.vector()

setdiff(annuals_nat,sp_leo) # J'ai quasiment tout relevé.
setdiff(annuals_fer,sp_leo) # J'ai quasiment tout relevé.
setdiff(pheno_lacking_fer,sp_leo)

#______________________________________________________
# Données Léo graines ####
day_month <- read.table("data/phenology/day_month.txt",header=T)
get_day_of_year <- function(Month,Day2){
  day_month %>% 
    filter(month == Month) %>% 
    filter(day_of_month == Day2) %>% 
    pull(day_of_year)
}

pheno_leo <- read.xlsx("data/phenology/Pheno_leo.xlsx",sheet="rawdata") %>% 
  mutate(code_sp = toupper(code_sp)) %>% 
  mutate(Day = as.Date(date- 25569, origin = "1970-01-01")) %>% 
  separate(col = Day,into = c("Year","Month","Day2")) %>% 
  mutate(Year = as.numeric(Year), Month = as.numeric(Month), Day2 = as.numeric(Day2)) %>% 
  mutate(Disp = map2_dbl(Month,Day2,get_day_of_year)) %>% 
  mutate(treatment = if_else(plot %in% c("C1","C2"),true = "Fer", false = "Nat")) %>% 
  mutate(Treatment = if_else(plot %in% c("C1","C2"),true = "Fer_Clc", false = "Nat_Sab")) 

ggplot(pheno_leo,aes(x=treatment,y=Disp))+
  geom_boxplot()

sp_info <- MEAN_pheno %>% 
  mutate(from = "database") %>% 
  select(species,code_sp,LifeHistory,LifeForm1) %>% 
  unique()

pheno_leo2 <- pheno_leo %>% 
  select(code_sp,treatment,Disp) %>% 
  group_by(code_sp, treatment) %>% 
  summarize(Disp_leo = mean(Disp))

pheno_merged <- MEAN_pheno %>% 
  merge(pheno_leo2,by=c("code_sp","treatment"))

ggplot(pheno_merged,aes(x=Disp,y=Disp_leo,label = code_sp))+
  geom_point() +
  geom_label()+
  facet_wrap(~treatment) +
  geom_abline(slope = 1, intercept = 0) 

cor()

ggsave("outputs/plots/Phenology_of_dispersal.png",height = 12, width = 12)


# Mesure par mesure #####
data_file <- "LaFage_PlantTraitsDP_vp.xlsx"
Pheno <- read.xlsx(paste0("data/traits/",data_file), sheet = "Pheno", startRow = 1, colNames = TRUE) %>% 
  select(Treatment,Year,Code_Sp,LifeForm1,Disp,measurementDeterminedBy)

PHENO <- pheno_leo %>% 
  rename(Code_Sp = code_sp) %>% 
  mutate(measurementDeterminedBy = "Leo Delalandre", LifeForm1 = "The") %>% 
  select(Treatment,Year,Code_Sp,LifeForm1,Disp,measurementDeterminedBy) %>% 
  rbind(Pheno) %>% 
  mutate(treatment = str_sub(Treatment,1,3)) %>% 
  filter(!(treatment == "Tem"))

PHENO_sp_leo <- PHENO %>% 
  filter(Code_Sp %in% sp_leo)

ggplot(PHENO_sp_leo, aes(x= measurementDeterminedBy, y = Disp, color = treatment ) ) +
  geom_point() +
  facet_wrap(~Code_Sp) +
  theme(axis.text.x = element_text(angle = 90))


# Environment ####
day_mean_month <- day_month %>% 
  group_by(month) %>% 
  summarize(mean_day_of_year = mean(day_of_year))

envir <- read.table("data/phenology/monthly_data_envir.txt",header=T, sep= "\t") %>% 
  mutate(month = c(1:12)) %>% 
  merge(day_mean_month,by="month") 


# Regarder le lien entre le moment de dispersion et le moment de l'année





 
      