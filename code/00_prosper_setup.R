##TITLE: PROSPER Data Set-Up
##AUTHOR: C. E. Moore
##Updated on AUG 2022

#PURPOSE: Format and prepare PROSPER dataset for use in RCS calculation

# SET UP ----

# .libPaths(.libPaths()[3:1]) #ARC
library(tidyverse); library(sf) 

#paths
PATH <- "~/Documents/Projects/PNW_fishes" #local
# PATH <- "/data" #ARC

PATH_prosp <- paste0(PATH, "/PROSPER")
PATH_nhd <- paste0(PATH, "/NHD")

#read in PROSPER CSVs
prosp <- lapply(list.files(paste0(PATH, "/PROSPER/FinalRevisedPublishedTables"), ".csv", full.names = TRUE), read.csv)
# subset for testing
prosp <- prosp[c(1,2)]

#combine
prosp.only <- bind_rows(prosp)
cat(length(unique(prosp.only$COMID)), "unique flowlines with prosper data \n")

#read in NHD data
nhd <- st_read(PATH_nhd, layer = "PNW_flowline")

#SUM + FILT. PROSPER ----

prosp.only <- prosp.only %>%
  mutate(summer_flow = JulStreamflow + AugStreamflow + SepStreamflow,
         annual_flow = rowSums(select(., matches('*Streamflow'))),
         lowQ = summer_flow / annual_flow)
#save with flow values
write.csv(prosp.only, file = paste0(PATH_prosp, "/prosper_baseflow_", Sys.Date(), ".csv"), row.names = F)

#summarize to across years at each flowline
p <- prosp.only %>%
  select(COMID, AugTemp, MEANspc, STDspc, MEANspp, STDspp, lowQ) %>%
  na_if(-9999) %>%
  group_by(COMID) %>%
  mutate(yr_MEANspc = mean(MEANspc, na.rm = T),
         yr_STDspc = sd(MEANspc, na.rm = T),
         yr_MEANspp = mean(MEANspp, na.rm = T),
         yr_STDspp = sd(MEANspp, na.rm = T),
         yr_mn_AugTemp = mean(AugTemp, na.rm = T),
         yr_sd_AugTemp = sd(AugTemp, na.rm = T),
         yr_mn_lowQ = mean(lowQ, na.rm = T),
         yr_sd_lowQ = sd(lowQ), na.rm = T) %>%
  select(COMID, matches('yr*')) %>%
  distinct()

#COMB. NHD + SUMPROSP ----

prosp.nhd <- nhd %>% select(COMID, LENGTHKM, AreaSqKM, TotDASqKM, geometry) %>%
  right_join(., p, by = "COMID") %>% st_zm()
# prosp.nhd[is.na(prosp.nhd)] <- 0
st_write(prosp.nhd, paste0(PATH_prosp, "/combined_streamdata"), "prosper_nhd_flowline", driver = "ESRI Shapefile", append = FALSE, delete_layer = TRUE)
cat("Averaged prosper data w/ nhd saved ...  /n")
