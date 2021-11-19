source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
source("scripts/functions/traitRange.test.R")
source("scripts/functions/traitMoments.test.R")

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8")) %>% 
  mutate(plot = paste(depth,paddock,sep="_"))

df <- read.xlsx("data/traits/data traits pour article JEcol 2012.xlsx", sheet = 1, startRow = 1, colNames = TRUE) %>% 
  filter(!lifeform=="Th") %>% 
  remove_rownames() %>% 
  column_to_rownames("sp") %>% 
  select(-c(code_sp,family,lifeform,lifestyle))

comm <-  read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", sheet = "abondances par parcelle", 
                   startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  column_to_rownames("plot") %>% 
  select(!!rownames(df))

dim(comm) # 77 comme dans le papier


  
nreps = 9999

trait <- "SLA"

null_range <- null.range(comm,trait,df,nreps)
  
null_range2 <- null_range %>% 
  mutate(ES = 2*(P.lower.range-0.5)) %>% 
  mutate(plot=rownames(comm)) %>% 
  full_join(soil_Maud,by="plot")

ggplot(null_range2,aes(x=PC1score,y=ES)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_segment(aes(x = PC1score, y = 0, xend = PC1score, yend = ES))+
  ggtitle(trait)
