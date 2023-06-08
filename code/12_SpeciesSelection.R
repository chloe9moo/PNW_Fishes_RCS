##TITLE: Calc info helpful for species selection
##AUTHOR: C. E. Moore
##Updated on 22 AUG 2022

#PURPOSE: 1) Find pct of each species range in study area

# set up ----

library(tidyverse); library(sf) #local

#paths
PATH <- "~/Documents/Projects/PNW_fishes" #local
# PATH <- "/home/chloe9mo/PNW_fishes" #ARC
PATH_occdat <- paste0(PATH, "/occurrence_data")


#calc pct range in study area ----

##load in species and make dataframe
occ.dat <- read.csv(paste0(PATH_occdat, "/occ_500mdist_no_outliers.csv"))
species.list <- unique(occ.dat$species) %>% as.data.frame() %>% add_column(tot_range_km = NA, pct_range = NA) %>% rename(species = '.')

##load in study area
library(maps)
pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon")) %>% 
  st_union() %>% st_sf()
# pnw.area <- as.numeric(st_area(pnw)) / 1e6 #calculate area of study area and convert to km2, but now I'm realizing this wasn't necessary

##+ calc for species with iucn range maps ----
for(i in 27:44) {
  
  shp <- st_read(dsn=paste0(PATH, "/range_maps/", i), layer = "data_0") 
  sp <- unique(shp$BINOMIAL) #pull out species name
  shp <- shp %>% st_union() %>% st_sf()
  range.area <- as.numeric(st_area(shp)) / 1e6 #calc total range area in km2
  species.list[species.list$species == sp, 'tot_range_km'] <- range.area
  pnw.shp <- st_intersection(pnw, shp) %>% 
    mutate(intersect_area = as.numeric(st_area(.))/1e6) %>% #get total overlapping area
    st_drop_geometry() %>% select(intersect_area) 
  species.list[species.list$species == sp, 'pct_range'] <- pnw.shp$intersect_area / range.area * 100
  
}
write.csv(species.list, file=paste0(PATH_occdat, "/pct_range_studyarea.csv"), row.names = FALSE)

##+ calc for species without iucn range maps ----
library(rgbif)

species.list <- read.csv(paste0(PATH_occdat, "/pct_range_studyarea.csv")) #reload just in case

#get species keys
spec.keys <- occ.dat %>% select(species, speciesKey) %>% 
  unique() %>% 
  left_join(., species.list, by = 'species') %>% 
  filter(is.na(pct_range)) %>%
  select(species, speciesKey)

#download gbif data
gbif_download <- occ_download(
  #step 1 - get occ for a species
  pred_in("taxonKey",  spec.keys$speciesKey), 
  pred("hasCoordinate", TRUE), #step 2 (coordinates avail.)
  pred("hasGeospatialIssue", FALSE), #partial step 3 (remove issues)
  pred_not(pred("issue", "IDENTIFIED_DATE_UNLIKELY")), #step 3 (remove issues)
  pred_not(pred("issue", "RECORDED_DATE_MISMATCH")), #step 3 (remove issues)
  pred_not(pred("occurrenceStatus", "ABSENT")), #step 4 (keep only 'presence' occurrences, also 'US'?)
  pred_notnull("country"), #step 5
  #pred_within(na_land_wkt), #step 6 - to pull occurrences within the NA shapefile, leftover from anuran RCS, can add back in
  #pred_lte("eventDate", date), #pull up to current date
  #pred_lt("coordinateUncertaintyInMeters", 2000), #removing uncertain pts, not doing this for now
  format = "SIMPLE_CSV"
)
#user=user, pwd=pwd, email=email) set these in the Renvironment

occ_download_wait(gbif_download) 
occ_download_get(gbif_download[1], path = PATH_occdat)
unzip(paste0(PATH_occdat, sep = "/", gbif_download[1],".zip"), exdir=PATH_occdat)
fish_occ<-read_tsv(paste0(PATH_occdat, sep="/", gbif_download[1],".csv"),
                   col_types="cccccccccccccccccccicnnnnnnnnTnnnccccccccTcccccTcc")
# fish_occ<-read_tsv(paste0(PATH_occdat, "/0069282-210914110416597.csv"),
#                    col_types="cccccccccccccccccccicnnnnnnnnTnnnccccccccTcccccTcc")

#calc area w/ convex hull for each species
for (i in 1:nrow(spec.keys)) {
  
  sp <- spec.keys[i, 'species'] #get the species
  sp.occ <- fish_occ %>% filter(species == sp) %>% select(species, decimalLongitude, decimalLatitude) %>%
    st_as_sf(coords=c('decimalLongitude', 'decimalLatitude')) %>% st_set_crs(., st_crs(pnw)) %>% #convert to sf pts
    st_buffer(50000) %>% #do 50 km buffer on the pts
    st_union() %>% st_sf() %>% st_make_valid() %>% filter(!st_is_empty(.))
  range.area <- as.numeric(st_area(sp.occ)) / 1e6 #calc total range area in km2
  species.list[species.list$species == sp, 'tot_range_km'] <- range.area
  pnw.shp <- st_intersection(sp.occ, pnw) %>% 
    mutate(intersect_area = as.numeric(st_area(.))/1e6) %>% #get total overlapping area
    st_drop_geometry() %>% select(intersect_area) 
  species.list[species.list$species == sp, 'pct_range'] <- pnw.shp$intersect_area / range.area * 100
  write.csv(species.list, file=paste0(PATH_occdat, "/pct_range_studyarea.csv"), row.names = FALSE)
  cat(sp, " complete...")
  cat("\n")
  
}


# make spp selection table ----

## order of filtering
#1. Remove points outside PNW extent
#2. Remove points outside year range (1995-2015)
#3. Removal of points based on expert review
#4. Removal of points > 500 m from a PROSPER stream
#5. Removal of points flagged as > 3 SD from mean of all env. variables used

ran <- read.csv(paste0(PATH_occdat, "/pct_range_studyarea.csv"))

#1. Points in PNW extent
pnw.filt <- read.csv(paste0(PATH_occdat,"/filter_tracking_datasets/PNWextent_ALLocc_2021-08-19.csv")) %>%
  select(species) %>%
  group_by(species) %>%
  mutate(one_pnw_ct = n()) %>%
  distinct()
#2. Points in year range
year.filt <- read.csv(paste0(PATH_occdat,"/filter_tracking_datasets/1995_2015_ALLocc_2021-08-19.csv")) %>%
  select(species) %>%
  group_by(species) %>%
  mutate(two_year_ct = n()) %>%
  distinct()
#3. Points confirmed by expert review
exp.filt <- read.csv(paste0(PATH_occdat,"/filter_tracking_datasets/Occ_PostExpQAQC.csv")) %>%
  select(species) %>%
  group_by(species) %>%
  mutate(three_expert_ct = n()) %>%
  distinct()
#4. Points within 500 m of stream
dist.filt <- read.csv(paste0(PATH_occdat,"/filter_tracking_datasets/Occ_PostExpQAQC.csv")) %>%
  filter(dist2strm_km < 0.5) %>%
  select(species) %>%
  group_by(species) %>%
  mutate(four_dist_ct = n()) %>%
  distinct()
#5. Points not > 3 SD of environmental variable means (also the final count!)
spp.cts <- read.csv(paste0(PATH_occdat, "/occ_500mdist_no_outliers.csv")) %>%
  select(species) %>%
  group_by(species) %>%
  mutate(final_no_out_ct = n()) %>%
  distinct() %>% 
  right_join(., ran, by="species") %>%
  right_join(., dist.filt) %>%
  right_join(., exp.filt) %>%
  right_join(., year.filt) %>%
  right_join(., pnw.filt) %>%
  select(-c(tot_range_km)) %>%
  arrange()
# spp.cts <- spp.cts[,c(1:3,7,6,5,4)]

write.csv(spp.cts, file=paste0(PATH_occdat, "/RecordCts_SpeciesSelection.csv"), row.names = F)


#visualize ----

##+ make histograms of species x year ----
# occ.dat <- read.csv(paste0(PATH_occdat, "/NHD_near_distance_occurrences_2021-11-02.csv"))
# occ.dat <- occ.dat %>% filter(dist2strm_km < 0.5)
# species.list <- unique(occ.dat$species) %>% as.data.frame()
# 
# pdf(file = paste0(PATH_occdat, "/SpeciesxYear_Histograms.pdf"))
# for(i in 1:nrow(species.list)) {
# 
#   sp <- species.list[i,]
#   p <- occ.dat %>% select(species, year) %>% filter(species == sp) %>%
#   ggplot(., aes(x=year)) + geom_histogram() + facet_wrap(species ~ .)
#   plot(p)
# 
# }
# dev.off()
##+ count of pts lost at diff. stream distance cutoffs ----

#read in occurrence data
occ.dat <- read.csv(paste0(PATH_occdat, "/filter_tracking_datasets/Occ_PostExpQAQC.csv"))

#make table
occ.lost <- as.data.frame(matrix(nrow = 5, ncol = 4, dimnames=list(NULL, c("filter_distance", "n_species", "pts_lost", "tot_pts"))))
occ.lost$filter_distance <- c("100 m", "500 m", "1 km", "5km", "none")

#no filter
occ.lost$n_species[occ.lost$filter_distance == "none"] <- length(unique(occ.dat$species))
occ.lost$pts_lost[occ.lost$filter_distance == "none"] <- 0
occ.lost$tot_pts[occ.lost$filter_distance == "none"] <- nrow(occ.dat)

#5 km
occ.5km <- occ.dat %>% filter(dist2strm_km < 5)
occ.lost$n_species[occ.lost$filter_distance == "5km"] <- length(unique(occ.5km$species))
occ.lost$pts_lost[occ.lost$filter_distance == "5km"] <- nrow(occ.dat) - nrow(occ.5km)
occ.lost$tot_pts[occ.lost$filter_distance == "5km"] <- nrow(occ.5km)

#1 km
occ.1km <- occ.dat %>% filter(dist2strm_km < 1)
occ.lost$n_species[occ.lost$filter_distance == "1 km"] <- length(unique(occ.1km$species))
occ.lost$pts_lost[occ.lost$filter_distance == "1 km"] <- nrow(occ.dat) - nrow(occ.1km)
occ.lost$tot_pts[occ.lost$filter_distance == "1 km"] <- nrow(occ.1km)

#500 m
occ.500m <- occ.dat %>% filter(dist2strm_km < 0.5)
occ.lost$n_species[occ.lost$filter_distance == "500 m"] <- length(unique(occ.500m$species))
occ.lost$pts_lost[occ.lost$filter_distance == "500 m"] <- nrow(occ.dat) - nrow(occ.500m)
occ.lost$tot_pts[occ.lost$filter_distance == "500 m"] <- nrow(occ.500m)

#100 m
occ.100m <- occ.dat %>% filter(dist2strm_km < 0.1)
occ.lost$n_species[occ.lost$filter_distance == "100 m"] <- length(unique(occ.100m$species))
occ.lost$pts_lost[occ.lost$filter_distance == "100 m"] <- nrow(occ.dat) - nrow(occ.100m)
occ.lost$tot_pts[occ.lost$filter_distance == "100 m"] <- nrow(occ.100m)

write.csv(occ.lost, file = "dist_nearest_strm_filtertracking.csv")  

#plot
ol <- occ.lost %>% 
  mutate(distance = factor(filter_distance, levels = filter_distance)) %>%
  select(-filter_distance) %>% 
  gather(type, count, tot_pts:pts_lost) %>% 
  ggplot(., aes(x=distance, y=count, fill=type)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_discrete(name=" ", labels=c("Points Lost", "Total Points")) +
  xlab("Filtered Distance to Nearest Stream") + ylab("N Occurrence Pts")
ggsave("filter_distance_barplot.png", plot=ol)


##+ where are outliers occurring? ----
library(maps); library(mapview)
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
occ.dat <- read.csv(paste0(PATH_occdat, "/occ_500mdist_outliers_flagged.csv"))
dat.sp <- st_as_sf(occ.dat, coords=c("decimalLongitude", "decimalLatitude"), crs=crs.geo)

pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon"))

spp <- "Cottus gulosus"

ggplot() +
  geom_sf(data = pnw, fill = "gray") +
  geom_sf(data = dat.sp %>% filter(species==spp), aes(fill = out_sd), shape=21, size=3, alpha=0.7) +
  scale_fill_manual(values = c("TRUE"= "#FF0000", "FASLE"="#808080")) +
  theme_bw() + xlab("Longitude") + ylab("Latitude")
mapview((dat.sp %>% filter(species==spp))["out_sd"])
