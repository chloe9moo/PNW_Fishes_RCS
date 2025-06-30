## 02. Background Pt Creation
## Project: PNW Fishes SDMs
## Author: CE Moore
## Date: 4 AUG 2022

## Purpose: Create background 'points' matching the SWD format of the occurrence records

#SET UP ----

.libPaths(.libPaths()[3:1]) #for ARC only
library(tidyverse); library(sf); library(terra)

#home
# PATH <- "~/Documents/projects/PNW_fishes" #local
PATH <- "/data" #ARC

#load data
species.list <- read.csv(paste0(PATH, "/occurrence_data/RCS_FocalSpecies.csv")) #species w/ >10 pts & >25% of their range in the study extent
occ <- read.csv(paste0(PATH, "/occurrence_data/occurrence_data_finalfilter.csv")) %>%
  filter(species %in% species.list$Species)

# env <- st_read(dsn = paste0(PATH, "/PROSPER/combined_streamdata"), layer = "prosper_nhd_flowline")
env <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")
huc <- st_read(dsn = paste0(PATH, "/HUC"), layer = "PNW_huc12") %>% 
  st_transform(st_crs(env)) %>%
  st_filter(., env)


#make bg points ----
pt.track <- data.frame(Species = species.list$Species) %>% mutate(N_occ = NA, N_bg_spc = NA, N_bg_r = NA)

for (i in 1:nrow(species.list)) {
  spp <- species.list[i, "Species"] #select species
  spp.occ <- occ %>% filter(species == spp) %>% st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(env)) #make spatial points

  pt.track$N_occ[pt.track$Species == spp] <- nrow(spp.occ)
  
  spp.huc <- st_filter(huc, spp.occ) #select hucs that species points fall in
  sp.huc.bf <- spp.huc %>% st_union() %>% st_buffer(dist = 5000) #add 5 km buffer to HUCs (sort of arbitrary, looked up dispersal of river fishes, found Comte & Olden paper which saw range of 1 km to 12 km :|)
  nhd.huc <- st_filter(env, sp.huc.bf) #select nhd flowlines that fall within the HUC buffer

##watershed spatially controlled bg pts ----

  if (nrow(nhd.huc) > 10000) { 
  
    nhd.bg <- nhd.huc %>%
      filter(!COMID %in% spp.occ$COMID) %>% #don't include 'presence' flowlines 
      slice_sample(n = 10000)
  
  } else { 
  
    nhd.bg <- nhd.huc %>%
      filter(!COMID %in% spp.occ$COMID)
    
  }
  pt.track$N_bg_spc[pt.track$Species == spp] <- nrow(nhd.bg)

  st_write(nhd.bg, dsn = paste0(PATH, "/SDM/Input_Data/Bg_Points/", gsub(" ", "_", spp), "_SpCont_BG.shp"), delete_layer = TRUE, quiet = T)

##randomly selected bg pts ----

  spp.hull <- spp.occ %>% st_union() %>% st_convex_hull() %>% st_buffer(dist = 50000)
  nhd.hull <- st_filter(env, spp.hull)
  
  if (nrow(nhd.hull) > 10000) { 
  
    r.nhd.bg <- nhd.hull %>%
      filter(!COMID %in% spp.occ$COMID) %>%
      slice_sample(n = 10000)
    
  } else { 
    
    r.nhd.bg <- nhd.hull %>%
      filter(!COMID %in% spp.occ$COMID)
    
  }
  pt.track$N_bg_r[pt.track$Species == spp] <- nrow(r.nhd.bg)

  st_write(r.nhd.bg, dsn = paste0(PATH, "/SDM/Input_Data/Bg_Points/", gsub(" ", "_", spp), "_Rand_BG.shp"), delete_layer = TRUE, quiet = T)
  
  cat( ifelse(i == nrow(species.list), "=| complete \n", ifelse(i == 1, "|=", "=")) )

}

write.csv(pt.track, file=paste0(PATH, "/SDM/Input_Data/Occ_Bg_Ct_Tracksheet.csv"), row.names = F)
