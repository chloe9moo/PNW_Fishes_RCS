##TITLE: Area of Occupancy Code
##AUTHOR: C. E. Moore; modified from code by S. Silknetter + T. DuBose
##Updated on 22 AUG 2022

#PURPOSE: Using species occurrence records and HUC12s, calculate the area of occurrence for each species. 
##To be used to calculate RCS index with climate sensitivity

#set up ----

#packages
library(sf); library(tidyverse)

#paths
PATH <- "~/Documents/Projects/PNW_fishes"
##in
PATH_occdat <- paste0(PATH, "/occurrence_data")
PATH_HUC12 <- paste0(PATH, "/HUC/PNW_huc12.rds")
##out
PATH_AOO_out <- paste0(PATH,"/RCS_results/AOO")
PATH_AOO_huc12 <- paste0(PATH_AOO_out, "/huc12_sp")
PATH_AOO_buffer <- paste0(PATH_AOO_out, "/buffer_sp")

#spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

#read in data
##HUC
huc12 <- readRDS(PATH_HUC12)

##occurrence data -- using data that includes expert review, pts <500 m from prosper streams, outliers removed
species.list <- read.csv(paste0(PATH, "/occurrence_data/RCS_FocalSpecies.csv")) #target species
occ.dat <- read.csv(paste0(PATH_occdat, "/occ_500mdist_no_outliers.csv")) %>%
  filter(species %in% species.list$Species) #only spp with > 10 occurrences

#set AOO function ----

find_area <- function(x){
  
  taxa <- unique(x$species)
  
  #read in sp occurrence data
  geodata <- occ.dat %>%
    filter(species == taxa) %>%
    dplyr::select(species, decimalLongitude, decimalLatitude, year, dataset_ID)
  tot.pt <- nrow(geodata)
  
  #convert from 'geodata' CSV to a sf object. 
  dat_sp <- st_as_sf(geodata, coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo) 
  
  #calculate area occupied per species.
  ##HUC12 Area
  dat_sp <- st_join(dat_sp, huc12)  %>% # Store the watershed name as an attribute of the data.
    rename(watershed.name=all_of('huc12')) # renaming watershed name column
  N_watersheds <- unique(dat_sp$watershed.name) #identify which watersheds are occupied
  na.pts <- length(is.na(dat_sp$watershed.name)[is.na(dat_sp$watershed.name) == TRUE])
    
  occHUC_sp <- huc12 %>%
    rename(watershed.name=all_of("huc12")) %>%
    filter(watershed.name %in% N_watersheds) # Extract unique HUC12 data for the species.
    
  huc.area <- sum(occHUC_sp$WSAREA, na.rm=T)  # Sum the total area (square kilometer) for each unique HUC12.
  occHUC_sp %>%
    select("WSAREA","watershed.name","shape_Length","shape_Area", "shape") %>%
    st_write(., dsn = paste0(PATH_AOO_huc12, sep = "/", gsub(" ", "_", taxa), '.shp'), append=F)

  na.wsh <- sum(is.na(occHUC_sp$WSAREA))
  # print(paste(length(is.na(N_watersheds)[is.na(N_watersheds) == FALSE]), class(st_geometry(occHUC_sp))))
  
  ##Buffer area
  #buffer expects projected CRS
  dat_sp <- dat_sp %>% filter(!is.na(dat_sp$watershed.name)) #remove any points that don't fall within a HUC
  dat_sp <- st_transform(dat_sp, crs.albers)
  
  #1 km buffer, union buffer, calc total area
  buff_sp <- st_union(st_buffer(dat_sp, dist = 1))
  st_write(buff_sp, dsn = paste0(PATH_AOO_buffer, sep = "/", gsub(" ", "_", taxa), '.shp'), append=F) #save shapefile
  
  buff.area <- as.numeric(st_area(buff_sp))
    
  AOOs <- data.frame(scientific_name = taxa,
                     huc_area_sqkm = huc.area,
                     buff_area_sqkm = buff.area,
                     total_pts = tot.pt,
                     na_pts = na.pts,
                     na_ws = na.wsh,
                     stringsAsFactors = F)
  return(AOOs)
}

#calc AOO ----

occ.list <- split(occ.dat, occ.dat$species)
aoo.out <- NULL
# #subset for testing
# occ.list <- occ.list[c(1:5)]

#system.time(
aoo.out <- lapply(occ.list, find_area)
#)
#full species run locally = ~15 min, 7.6G mem

## parallel option // more mem, less time
# library(parallel)
# 
# n.cores <- detectCores() - 2
# 
# system.time(
# aoo.out <- mclapply(occ.list, find_area, mc.cores = n.cores)
# )
##5 species 15 sec, up to 13G mem

##unlist and save output
aoo.out <- bind_rows(aoo.out)
write.csv(aoo.out, file = paste0(PATH_AOO_out, "/AOO_", Sys.Date(), ".csv"), row.names = FALSE)
# aoo.out <- read.csv(paste0(PATH_AOO_out,"/AOO_2022-06-08.csv"))

#summarize aoo ----

library(scales)

AOOs_summarized <- aoo.out %>%
  #add ranks for quick comparison
  mutate(rank_WS=rank(huc_area_sqkm),
         rank_buff=rank(buff_area_sqkm),
         rescale(aoo.out$mean_Rank, to = c(0,1))) %>%
  rowwise() %>%
  mutate(mean_rank=mean(c(rank_WS, rank_buff)),
         sd_rank=sd(c(rank_WS, rank_buff)),
         log10_buff_sqkm=log10(buff_area_sqkm),
         log10_WS_sqkm=log10(huc_area_sqkm)) %>%
  select(-c(huc_area_sqkm, buff_area_sqkm, total_pts, na_pts, na_ws))
write.csv(AOOs_summarized, file = paste0(PATH_AOO_out, "/AOO_sumranks_", Sys.Date(), ".csv"), row.names = FALSE)

