## 01. Environmental Data Set Up
## Project: PNW Fishes SDMs
## Author: CE Moore
## Date: 25 AUG 2022

## Purpose: Summarize environmental data to HYBAS_ID for predictor set to run MaxEnt SDMs

#set up ----

# .libPaths(.libPaths()[3:1]) #for ARC only
library(tidyverse); library(sf); library(terra); library(exactextractr)

#HOME
PATH <- getwd() #local
# PATH <- "/data" #ARC

occ <- read.csv(paste0(PATH,"/occurrence_data/occ_500mdist_no_outliers.csv"))
# occ <- occ %>% select(dataset_ID, species, speciesKey, decimalLatitude, decimalLongitude, year, coordinateUncertaintyInMeters,
#                       individualCount, orig.species, family, coordinatePrecision, institutionCode, gbifID, COMID, dist2strm_km)
# # write.csv(occ, file = paste0(PATH, "/occurrence_data/occurrence_data_finalfilter.csv"), row.names = F)

#read in PROSPER CSVs
prosp <- lapply(list.files(paste0(PATH, "/PROSPER/FinalRevisedPublishedTables"), ".csv", full.names = TRUE), read.csv)
#subset for testing
# prosp <- prosp[c(1,2)]

#combine
prosp.only <- bind_rows(prosp)

##read in NHD data
nhd <- st_read(paste0(PATH, "/NHD"), layer = "PNW_flowline")

#sum. prosper ----

prosp.only <- prosp.only %>%
  select(-c(ROWID, COUNT, Year, OBSPRED_ID, region, AREA, MAJORITYspc, MEDIANspc, STDspc, STDspp)) %>%
  na_if(-9999) %>%
  mutate(COMID = as.character(COMID),
          summer_flow = rowMeans(select(., JunStreamflow, JulStreamflow, AugStreamflow, SepStreamflow), na.rm=T),
          winter_flow = rowMeans(select(., DecStreamflow, JanStreamflow, FebStreamflow), na.rm=T),
          annual_flow = rowSums(select(., matches('*Streamflow'))),
          lowQ = ifelse(annual_flow == 0, 0, 
                        ((JunStreamflow + JulStreamflow + AugStreamflow + SepStreamflow) / annual_flow))) %>%
  group_by(COMID) %>%
  mutate(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  ungroup() %>%
  unique()
  
##save with flow values
write.csv(prosp.only, file = paste0(PATH, "/SDM/Input_Data/prosper_all_SDMs.csv"), row.names = F)

#combine NHD+PROSPER ----

#NHD variables to include: Stream Order, catchment area, dist. to network end, dist upstream, slope (gradient), velocity
prosp.nhd <- nhd %>%
 st_zm() %>%
 select(COMID, StreamOrde, AreaSqKM, Pathlength, ArbolateSu, SLOPE, VA_MA) %>%
 mutate(COMID = as.character(COMID),
        SLOPE = na_if(SLOPE, -9998),
        VA_MA = na_if(VA_MA, -9998)) %>%
 right_join(., prosp.only, by = "COMID")

rm(prosp, nhd, prosp.only)

# prosp.nhd[is.na(prosp.nhd)] <- 0

st_write(prosp.nhd, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", append = FALSE)
cat("Prosper, nhd saved ... \n")

#FLASH prep ----
w95 <- st_read(dsn = paste0(PATH, "/NHD/S_USA.Hydro_FlowMet_1990s.gdb"))

nhd.sdm <- w95 %>%
 st_drop_geometry() %>%
 select(COMID, W95_HIST) %>%
 right_join(., prosp.nhd, by = "COMID") %>%
 st_as_sf()
# nhd.sdm <- st_as_sf(nhd.sdm)

rm(w95, prosp.nhd)

st_write(nhd.sdm, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", delete_layer = T)
cat("Prosper, nhd, W95 saved ... \n")

#PRISM prep ----

nhd.sdm <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")

##+ watershed avg. ----
huc <- st_read(dsn = paste0(PATH, "/HUC"), layer = "PNW_huc12") %>% st_transform(st_crs(nhd.sdm))
huc <- st_filter(huc, nhd.sdm)
# nhd.sdm <- st_filter(nhd.sdm, huc)

pri <- rast(list.files(paste0(PATH, "/PRISM"), pattern = "asc.asc$", full.names = T))
names(pri) <- c("ppt", "tmax", "tmin")
pnw_prism <- crop(pri, huc)

huc.val <- exact_extract(pnw_prism, huc, fun = 'mean', weights = "area", append_cols = "huc12", force_df = TRUE, stack_apply = TRUE)
huc <- left_join(huc, huc.val, by = "huc12")

n.h <- data.frame(COMID=character(), huc12=character())

for (i in 1:nrow(huc)) {

 h <- huc[i,] %>% select(huc12)
 n <- st_filter((nhd.sdm %>% select(COMID)), h) %>%
   st_drop_geometry() %>%
   mutate(huc12 = h$huc12)
 n.h <- rbind(n.h, n)

}

write.csv(n.h, file = paste0(PATH, "/SDM/Input_Data/HUCxNHD_PRISM_overlap.csv"), row.names = F)
n.h <- read.csv(paste0(PATH, "/SDM/Input_Data/HUCxNHD_PRISM_overlap.csv"))

n.h <- n.h %>%
 left_join(., huc, by="huc12") %>%
 select(-geometry, -huc12) %>%
 group_by(COMID) %>%
 mutate(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
 unique()

nhd.sdm <- nhd.sdm %>%
 left_join(., n.h, by="COMID") %>%
 select(-matches("*Strmf"), -AprStrm, -OctStrm)
names(nhd.sdm) <- c("COMID", "W95_HIS", "StrmOrd", "AreSqKM", "Pthlngt", "ArboltS", "SLOPE", "VA_MA", "CUMDRAI",
                   "MEANspc", "MEANspp", "AugTemp", "smmr_fl", "wntr_fl", "annl_fl", "lowQ", "WSAREA",
                   "mn.ppt.ws", "mn.tmax.ws", "mn.tmin.ws", "geometry")

st_write(nhd.sdm, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", delete_layer = T)
cat("Prosper, nhd, W95, prism/wtrshd, saved ... \n")

##+ buffer avg. ----
# nhd.sdm <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")

nhd.buff <- nhd.sdm %>%
 st_buffer(dist = 500)
pnw_prism <- crop(pri, nhd.buff)

buff.val <- exact_extract(pnw_prism, nhd.buff, fun = 'mean', weights = "area",
                         append_cols = "COMID", force_df = TRUE, stack_apply = TRUE)
names(buff.val) <- c("COMID", "mn.ppt.bf", "mn.tmax.bf", "mn.tmin.bf")

nhd.sdm <- left_join(nhd.sdm, buff.val, by="COMID")

st_write(nhd.sdm, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", delete_layer = T)
cat("Prosper, nhd, W95, prism/wtrshd, prism/buff saved ... \n")

#TCC prep ----
# nhd.sdm <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")

#full US
tcc <- rast(paste0(PATH,"/NLCD/usfs_carto_CONUS_2011/usfs_2011_treecanopy_cartographic_12-14-2018.img")) %>% project(crs(vect(nhd.sdm)))

#crop
library(maps)
pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon"))
pnw <- pnw %>% st_union() %>% st_buffer(dist=1000)
p.tcc <- crop(tcc, pnw)
rm(tcc)
#save
writeRaster(p.tcc, paste0(PATH,"/NLCD/PNW_tcc_cat.asc"), overwrite=T)
activeCat(p.tcc) <- 0
p.tcc <- as.numeric(p.tcc)
writeRaster(p.tcc, paste0(PATH,"/NLCD/PNW_tcc.asc"), overwrite=T)

p.tcc <- rast(paste0(PATH,"/NLCD/PNW_tcc.asc"))

buff.val <- exact_extract(p.tcc, nhd.buff, fun = 'mean', weights = "area",
                         append_cols = "COMID", force_df = TRUE, stack_apply = TRUE)
names(buff.val) <- c("COMID", "mn.tcc")

nhd.sdm <- left_join(nhd.sdm, buff.val, by="COMID")

st_write(nhd.sdm, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", delete_layer = T)
cat("Prosper, nhd, W95, prism/wtrshd, prism/buff, canopy cover saved ... \n")

#NLCD prep ----
nhd.sdm <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")

p.nlcd <- rast(paste0(PATH, "/NLCD/PNW_NLCD.asc"))
activeCat(p.nlcd) <- 0
# rcl.mat <- cbind(class=unique(p.nlcd$value), new.class=c(1, 2, 3, 4, 4, 4, 4, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9))
rcl.mat <- cbind(class=c(0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95), 
                 new.class=c(1, 2, 3, 4, 4, 4, 4, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9))

c.nlcd <- classify(p.nlcd, as.matrix(rcl.mat)) 
c.nlcd <- as.factor(c.nlcd) #reclassify and make single category raster

# huc1 <- huc[c(1,100),] #subset for testing

freqs <- exact_extract(c.nlcd, huc, function(value, coverage_fraction) {
  data.frame(value = value,
             frac = (coverage_fraction / sum(coverage_fraction, na.rm = T))*100) %>%
    group_by(value) %>%
    summarize(freq = sum(frac, na.rm = T), .groups = 'drop') %>%
    pivot_wider(names_from = 'value',
                names_prefix = 'freq_',
                values_from = 'freq')
}, append_cols = "huc12") %>% 
  mutate(across(starts_with('freq'), replace_na, 0))

huc <- left_join(huc, freqs, by = "huc12") %>%
  rename_with(., recode, 
              freq_1 = "p_uncl",
              freq_2 = "p_water",
              freq_3 = "p_snow",
              freq_4 = "p_dev",
              freq_5 = "p_barr",
              freq_6 = "p_forst",
              freq_7 = "p_grass",
              freq_8 = "p_agr",
              freq_9 = "p_wetl")

#n.h <- data.frame(COMID=character(), huc12=character())
#
#for (i in 1:nrow(huc)) {
#  
#  h <- huc[i,] %>% select(huc12)
#  n <- st_filter((nhd.sdm %>% select(COMID)), h) %>%
#    st_drop_geometry() %>%
#    mutate(huc12 = h$huc12)
#  n.h <- rbind(n.h, n)
#  
#}
#
#write.csv(n.h, file = paste0(PATH, "/SDM/Input_Data/HUCxNHD_overlap.csv"), row.names = F)
n.h <- read.csv(paste0(PATH, "/SDM/Input_Data/HUCxNHD_overlap.csv"))

n.h <- n.h %>%
  left_join(., huc, by="huc12") %>%
  select(-geometry, -huc12, -WSAREA, -mean.ppt, -mean.tmax, -mean.tmin) %>%
  group_by(COMID) %>%
  mutate(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  unique()

nhd.sdm <- nhd.sdm %>%
  left_join(., n.h, by="COMID")

st_write(nhd.sdm, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", delete_layer = T)
cat("Prosper, nhd, W95, prism/wtrshd, tcc, and nlcd saved ... \n")


#update PROSPER ----
nhd.sdm <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline") #saved data w/ extracted variables

#updated prosper data
prosp <- lapply(list.files(paste0(PATH, "/PROSPER/FinalRevisedPublishedTables"), ".csv", full.names = TRUE), read.csv)
#combine
prosp.only <- bind_rows(prosp)
# prosp.only <- prosp.only %>% filter(COMID %in% nhd.sdm$COMID)

#summarize the same way
prosp.only <- prosp.only %>%
  select(-c(X, ROWID, COUNT, Year, OBSPRED_ID, region, AREA, MAJORITYspc, MEDIANspc, STDspc, STDspp)) %>%
  na_if(-9999) %>%
  mutate(COMID = as.character(COMID),
         summer_flow = rowMeans(select(., JunStreamflow, JulStreamflow, AugStreamflow, SepStreamflow), na.rm=T),
         winter_flow = rowMeans(select(., DecStreamflow, JanStreamflow, FebStreamflow), na.rm=T),
         annual_flow = rowSums(select(., matches('*Streamflow')), na.rm = T),
         lowQ = ifelse(annual_flow == 0, 0, 
                       ((JunStreamflow + JulStreamflow + AugStreamflow + SepStreamflow) / annual_flow))) %>%
  group_by(COMID) %>%
  mutate(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  ungroup() %>%
  distinct()

##save with flow values
write.csv(prosp.only, file = paste0(PATH, "/SDM/Input_Data/prosper_all_SDMs.csv"), row.names = F)
rm(prosp, nhd, prosp.only)

#combine new prosp with all other env var
prosp.only <- read.csv(paste0(PATH, "/SDM/Input_Data/prosper_all_SDMs.csv"))
prosp.only <- prosp.only %>% filter(COMID %in% nhd.sdm$COMID) %>% 
  select(-contains("Streamflow"), -AugTemp_SE) %>%
  mutate(COMID = as.character(COMID)) %>%
  rename(CUMDRAI=CUMDRAINAG, smmr_fl=summer_flow, wntr_fl=winter_flow, annl_fl=annual_flow)

#the new prosper data has fewer flowlines, however this shouldn't be a problem as only 1 occurrence point is nearest to them
nhd.sdm <- nhd.sdm %>% filter(COMID %in% prosp.only$COMID) %>%
  select(-colnames(prosp.only %>% select(-COMID))) %>%
  left_join(., prosp.only, by="COMID")

st_write(nhd.sdm, paste0(PATH, "/SDM/Input_Data"), "prosper_nhd_flowline", driver = "ESRI Shapefile", delete_layer = T)


#variable Corr ----
# nhd.sdm <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")
c <- cor((nhd.sdm %>% st_drop_geometry() %>% select(-COMID)), use = "complete.obs", method = "pearson") %>%
  as.data.frame() %>%
  round(., digits = 4)
write.csv(c, file = paste0(PATH, "/SDM/Input_Data/EnvVar_Correlations.csv"), row.names = T)

cat("Correlation csv saved... \n")

#check CRS by plot ----
# library(maps)
# sf.pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon"))
# png(paste0(PATH,"/SDM/BUFvsNHD2.png"), width = 1000, height = 1000)
# # plot(huc1$geometry, add=T)
# plot(sf.pnw$geom, col="blue")
# plot(nhd.sdm$geometry, col="red", lwd=12, add=T)
# plot(huc$geometry, col="black", add=T)
# dev.off()
