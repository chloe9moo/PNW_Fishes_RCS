## 10. MaxEnt Model Creation
## Project: PNW Fishes SDMs
## Author: CE Moore
## Date: 8 AUG 2022

## Purpose: Run species distribution models with MAXENT

#set up ----

.libPaths(.libPaths()[3:1]) #for ARC only
library(tidyverse); library(sf); library(ENMeval)

#home
# PATH <- getwd() #local
PATH <- "/data" #ARC

##load data ----
species.list <- read.csv(paste0(PATH, "/occurrence_data/RCS_FocalSpecies.csv")) %>% filter(N_final >= 20)#species w/ >10 pts & >25% of their range in the study extent
occ <- read.csv(paste0(PATH, "/occurrence_data/occurrence_data_finalfilter.csv")) %>%
  filter(species %in% species.list$Species)

env.var <- read.csv(paste0(PATH, "/SDM/Input_Data/EnvPredictor_Table.csv")) %>% select(Variable.Name) %>% filter(!Variable.Name %in% c("smmr_fl", "wntr_fl"))
env <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")
env.swd <- env %>% 
  st_drop_geometry() %>%
  mutate(x = as.integer(COMID),
         y = as.integer(COMID)) %>%
  select(x, y, all_of(env.var$Variable.Name)) %>%
  relocate(x, y) %>%
  replace(is.na(.), 0)

occ <- occ %>% filter(COMID %in% env$COMID)

##select spp & bg type to run ----
spp <- species.list[1, "Species"]
b.t <- "SpCont"

cat("Running sdm for", spp, "using", b.t, "background points... \n")

##prepare inputs ----

sp.occ <- occ %>% 
  filter(species == spp) %>%
  dplyr::select(COMID) %>%
  left_join(., env.swd, by=c("COMID"="x")) %>%
  mutate(x = COMID,
         VA_MA = ifelse(VA_MA <= -9998, 0, VA_MA)) %>%
  select(-COMID) %>%
  relocate(x, y) %>%
  replace(is.na(.), 0)
  
bg <- st_read(dsn=paste0(PATH, "/SDM/Input_Data/Bg_Points"), layer=paste0(gsub(" ", "_", spp), "_", b.t, "_BG")) %>%
  st_drop_geometry() %>%
  mutate(x = as.integer(COMID),
         y = as.integer(COMID),
         VA_MA = ifelse(VA_MA <= -9998, 0, VA_MA)) %>%
  select(x, y, all_of(env.var$Variable.Name)) %>%
  relocate(x, y) %>%
  replace(is.na(.), 0)

#set parameters
p.arg <- list(fc = c("L", "LQ", "Q", "P", "T", "H"), rm = 1:5)

#run model ----
mx <- ENMevaluate(occs = sp.occ,
                  bg = bg,
                  tune.args = p.arg,
                  partitions = "block",
                  algorithm = "maxent.jar",
                  taxon.name = spp)
saveRDS(mx, file = paste0(PATH, "/SDM/Output/Model_RDS/", gsub(" ", "_", spp), "_ENMeval_", b.t, ".rds"))

cat("ENM evaluation complete + saved... \n")

res <- mx@results
write.csv(res, file = paste0(PATH, "/SDM/Output/Results_Tables/", gsub(" ", "_", spp), "_res_", b.t, ".csv"), row.names = F)
cat("Results saved... \n")

##select top model ----

if (any(res$cbi.val.avg > 0.6, na.rm = T)) {
  
  #if cbi values are > 0.7, select model based on highest CBI, with tie breaker of highest AUC
  cat("selecting based on CBI \n")
  opt.mod <- res %>%
    filter(cbi.val.avg == max(cbi.val.avg, na.rm = T)) %>%
    filter(auc.val.avg == max(auc.val.avg, na.rm = T))

   } else {
     
     #if cbi values are low, the model probably needs to be reworked, but pick the highest AUC
     cat("selecting based on AUC only \n")
     opt.mod <- res %>%
       filter(auc.val.avg == max(auc.val.avg))

   }

opt.mod <- opt.mod[1,] #if there's still a tie, just take the first one

###pull out the model ----

top.mod <- mx@models[[opt.mod$tune.args]]
var.imp <- mx@variable.importance[[opt.mod$tune.args]]

write.csv(var.imp, file = paste0(PATH, "/SDM/Output/Results_Tables/", gsub(" ", "_", spp), "_", b.t, "_VarImp.csv"), row.names = F)
cat("Variable importance table saved... \n")

#visualize ----

# #example data for testing
# p <- read.csv(paste0(PATH, "/ExData/gadbis.csv")) %>% select(-sp)
# b <- read.csv(paste0(PATH, "/ExData/background.csv")) %>% select(-sp)
# e <- read.csv(paste0(PATH, "/ExData/rivers.csv"))

# evalplot.stats(e = mx, stats = "or.mtp", color = "fc", x.var = "rm")

