# 0) Number of annual species ####
# NB : must be done at the level of plot, I guess.
richness_per_guild_nat <- ab_nat %>% 
  count(depth,paddock,LifeHistory) %>% 
  merge(soil_Maud,by=c("depth","paddock"))

richness_per_guild_fer <- ab_fer %>%
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  mutate(depth = "Fer") %>% 
  count(depth,paddock,LifeHistory) %>% 
  mutate(PC1score = NA)

richness_per_guild <- rbind(richness_per_guild_nat,richness_per_guild_fer)
richness_per_guild$depth <- factor(richness_per_guild$depth , levels = c("Fer","D","I","S"))

richness_per_guild_toplot <- richness_per_guild %>% 
  spread(LifeHistory,n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(relative_richness_annual = annual / (annual + perennial))

boxplot(relative_richness_annual *100 ~ depth,
        data = richness_per_guild_toplot,
        # ylim=c(0,100),
        xlab = NA,
        ylab = "Relative abundance of annuals (%)",
        xaxt = "n",
        medlwd = 1)

# 1) Boxplot of annual cover ####
ann_fer <- ab_fer %>% 
  mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
  group_by(LifeHistory,year,paddock,id_transect_quadrat) %>% 
  dplyr::rename(line = id_transect_quadrat) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  filter(LifeHistory =="annual") %>% 
  mutate(depth = "Fer") %>% 
  relocate(LifeHistory,year,depth,paddock,line,tot_relat_ab)

ann_nat <- ab_nat %>% 
  group_by(LifeHistory,depth,paddock,line) %>% 
  summarise(tot_relat_ab = sum(relat_ab)) %>% 
  filter(LifeHistory =="annual") %>% 
  mutate(year = 2009) %>% 
  relocate(LifeHistory,year,depth,paddock,line,tot_relat_ab) %>% 
  mutate(line = as.character(line))

cover_annuals <- rbind(ann_fer,ann_nat)
cover_annuals$depth <- factor(cover_annuals$depth , levels = c("Fer","D","I","S"))


boxplot(tot_relat_ab * 100 ~ depth,
        data = cover_annuals %>% filter(!(depth == "Fer")),
        # ylim=c(0,100),
        xlab = NA,
        ylab = "Relative abundance of annuals (%)",
        xaxt = "n",
        medlwd = 1)
# Est-ce qu'il faut que je bosse au niveau de la ligne, ou du parc? Au niveau du parc, j'ai deux points pour le fertile,
# sauf si je prends plusieurs années.

# 2) Boxplots of CWM CSR scores ####
CSR_toplot_nat <- CWM2_nat %>% 
  select(depth,CWM_C,CWM_S,CWM_R) %>% 
  gather(key = "score", value = "value", -c(depth) )

CSR_toplot_fer <- CWM2_fer %>% 
  mutate(depth = "Fer") %>% 
  select(depth,CWM_C,CWM_S,CWM_R) %>% 
  gather(key = "score", value = "value", -c(depth) )

CSR_toplot <- rbind(CSR_toplot_nat,CSR_toplot_fer)
CSR_toplot$score <- factor(CSR_toplot$score , levels = c("CWM_C","CWM_S","CWM_R"))
CSR_toplot$depth <- factor(CSR_toplot$depth , levels = c("Fer","D","I","S"))


boxplot(
  value ~ score * depth , data = CSR_toplot, xaxt = "n",
  xlab = "", ylab = "Score (%)",
  col = c("black", "grey","white"),
  medlwd = 1
)

legend(
  "topleft", title = "Score",
  legend = c("C", "S", "R"), fill = c("black", "grey","white"),
  cex = 0.7
)


