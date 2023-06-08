##TITLE: GBIF + Other Occurrence Download & Filtering
##AUTHOR: C. E. Moore
##Updated on 22 AUG 2022

#PURPOSE: Pull GBIF occurrence records for target PNW fish species, filter, and record counts
#Also combine all datasets, and filter to PNW and within year range, keep record of counts
##modified from code by T. DuBose, S. Silknetter

#set up ----

library(rgbif); library(tidyverse); library(sf); library(readxl)

PATH <- "~/Documents/Projects/PNW_fishes"
PATH_occ <- paste0(PATH, "/occurrence_data")
date <- Sys.Date() #for pulling data up to current date and saving dated data

species.list <- read_excel(paste0(PATH_occ, "/Other_Datasets/Mimsetal2018_SelectedVertNetData.xlsx"), sheet = "6")
names(species.list) <- c("common_name", "species", "notes", "stream_species")

#include some non-stream species (ie listed as 'no') based on notes to keep
species.list <- species.list %>% mutate(stream_species = ifelse(species == "Prosopium coulterii", "Yes", stream_species))

#remove the 'no's and 'na's
species.list <- species.list[species.list$stream_species %in% "Yes",]


#gbif ----
##+ get species keys ----

pnw.gbif.taxa <- NULL
for(u in 1:nrow(species.list)){
  pnw.gbif.taxa1 <- name_backbone(name=species.list[u,2], 
                                   rank='species', kingdom='animals')
  pnw.gbif.taxa <- bind_rows(pnw.gbif.taxa1, pnw.gbif.taxa)
}

#remove synonyms with repeat species keys
pnw.gbif.taxa <- pnw.gbif.taxa %>% filter(!duplicated(speciesKey))


#+ pull data, run initial filter ----

ind_mat<-matrix(c(1,34,35,68), ncol=2, byrow = T) #using this to break up download into 2 chunks
fish_occ <-NULL #for occurrence data
gbif_codes <- NULL #for download codes

for(j in 1:2){
  gbif_download = occ_download(
    #step 1 - get occ for a species
    pred_in("taxonKey",  pnw.gbif.taxa$speciesKey[ind_mat[j,1]:ind_mat[j,2]]), 
    pred("hasCoordinate", TRUE), #step 2 (coordinates avail.)
    pred("hasGeospatialIssue", FALSE), #partial step 3 (remove issues)
    pred_not(pred("issue", "IDENTIFIED_DATE_UNLIKELY")), #step 3 (remove issues)
    pred_not(pred("issue", "RECORDED_DATE_MISMATCH")), #step 3 (remove issues)
    pred_not(pred("occurrenceStatus", "ABSENT")), #step 4 (keep only 'presence' occurrences, also 'US'?)
    pred_notnull("country"), #step 5
    #pred_within(na_land_wkt), #step 6 - to pull occurrences within the NA shapefile, leftover from anuran RCS, can add back in
    pred_lte("eventDate", date), #pull up to current date
    #pred_lt("coordinateUncertaintyInMeters", 2000), #removing uncertain pts, not doing this for now
    format = "SIMPLE_CSV"
  )
    #user=user, pwd=pwd, email=email) set these in the Renvironment
  
  occ_download_wait(gbif_download) 
  occ_download_get(gbif_download[1], path = PATH_occ)
  unzip(paste0(PATH_occ, sep = "/", gbif_download[1],".zip"), exdir=PATH_occ)
  occ1<-read_tsv(paste0(PATH_occ, sep="/", gbif_download[1],".csv"),
                 col_types="cccccccccccccccccccicnnnnnnnnTnnnccccccccTcccccTcc")
  
  fish_occ<-bind_rows(fish_occ, occ1)
  gbif_codes<-append(gbif_codes, gbif_download[1])
}

write.table(gbif_codes,
            paste0(PATH, sep = "/", "GBIF_codes_", date, ".txt"))

View(head(fish_occ)) #check

# #to load in saved csv's from download separately:
# X1<-read_tsv(paste0(PATH_occ, sep="/", "0325726-200613084148143.csv"),
#                col_types="cccccccccccccccc-ccicnnnnnnnnTnnnccccccccTcccccTcc")
# X2<-read_tsv(paste0(PATH_occ, sep="/", "0325730-200613084148143.csv"),
#                col_types="cccccccccccccccc-ccicnnnnnnnnTnnnccccccccTcccccTcc")
# fish_occ <- rbind(X1, X2)

# remove occ with NA coordinates
fish_occ <- fish_occ %>% filter(!is.na(decimalLatitude)) #na seems to be data incorrectly uploaded/shifted for future reference

# count occurrences pulled down
filter_count <- fish_occ %>% 
  group_by(species) %>%
  mutate(basic_filter = n()) %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) %>%
  mutate(lat_long_norep = n()) %>%
  slice(1) %>% select(species, basic_filter, lat_long_norep)
#lost 1 species in the pull!

# remove the duplicates of lat + long + year
fish_occ <- fish_occ %>% 
  group_by(species) %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE)

# add the lost species into the count
pnw.gbif.taxa$species[!(pnw.gbif.taxa$species %in% filter_count$species)] #Cottus bendirei lost
filter_count[nrow(filter_count) + 1, 1] <- "Cottus bendirei"
filter_count[is.na(filter_count)] <- 0

filter_count <- filter_count %>% arrange(species)

#save full set w/ basic filter applied
write.csv(fish_occ, paste0(PATH_occ, "/fish_GBIFocc_full_", date, ".csv"), row.names = FALSE)

fish_occ <- read.csv(paste0(PATH_occ, "/GBIF/fish_GBIFocc_full_2021-07-16.csv"))


##+ filter PNW extent only ####

library(maps)

#spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

#prep boundary + background
pnw <- maps::map('state', regions = c("Oregon", "Washington", "Idaho"), fill = TRUE, col="transparent") #load states of interest
sf.pnw <- st_as_sf(pnw, crs.geo)

#make data spatial
dat_sp <- st_as_sf(fish_occ, coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo)
dat_alb <- dat_sp %>% st_transform(crs.albers) #N. American projected crs

#10 km buffer around extent
pnw.buff <- sf.pnw %>% st_union() %>% st_transform(crs.albers) %>% st_buffer(dist=10)

#clip data
fish_occ.pnw <- dat_alb[pnw.buff,]  

#save counts
pnw_count <- fish_occ.pnw %>% st_drop_geometry() %>% group_by(species) %>% count(.) 
filter_count <- merge(filter_count, pnw_count, by="species", all.x=TRUE)
colnames(filter_count)[4] <- "pnw_only"

#save filtered occurrences only in PNW extent
fish_occ.pnw <- st_transform(fish_occ.pnw, crs.geo)
st_write(fish_occ.pnw, dsn = paste0(PATH_occ, "/GBIF/PNWextent_GBIFocc_", date, ".csv"), delete_layer = T, layer_options = "GEOMETRY=AS_XY")


##+ filter 1995 - 2015 ----

occ <- read.csv(paste0(PATH_occ, "/GBIF/PNWextent_GBIFocc_2021-07-16.csv"))
yr.occ <- occ %>% filter(year >= 1995 & year <= 2015)

#save counts
yr_count <- yr.occ %>% group_by(species) %>% count(.) 
filter_count <- merge(filter_count, yr_count, by="species", all.x=TRUE)
colnames(filter_count)[5] <- "yr_1995_2015"

#save filtered occurrences only in PNW extent and 1995 - 2015
write.csv(yr.occ, file = paste0(PATH_occ, "/GBIF/1995_2015_GBIFocc_", date, ".csv"), row.names=FALSE)


##+ filter 2004 - 2015 ####

yr.occ2 <- yr.occ %>% filter(year >= 2004 & year <= 2015)

#save counts
yr_count2 <- yr.occ2 %>% group_by(species) %>% count(.)
filter_count <- merge(filter_count, yr_count2, by="species", all.x=TRUE)
colnames(filter_count)[6] <- "yr_2004_2015"

#save filtered occurrences only in PNW extent and 1995 - 2015
write.csv(yr.occ2, file = paste0(PATH_occ, "/GBIF/2004_2015_GBIFocc_", date, ".csv"), row.names=FALSE)

#save filter counts
filter_count[is.na(filter_count)] <- 0
write.csv(filter_count, file = paste0(PATH, "/GBIF_filter_tracksheet_", date, ".csv"))


#combined ----

#read in datasets
pnw.gbif <- read.csv(paste0(PATH_occ, "/GBIF/PNWextent_GBIFocc_2021-07-16.csv"))
columbia.basin <- read_excel(paste0(PATH_occ, "/Other_Datasets/OJSMNH Columbia Basin specimens_15Mar2021.xlsx"))
vertnet18 <- read_excel(paste0(PATH_occ, "/Other_Datasets/Mimsetal2018_SelectedVertNetData.xlsx"), sheet = "1", 
                        col_types = c("text", "text", "text", "numeric", "text", "numeric", "numeric", "numeric", "text", "text",
                                      "text", "text", "text", "text", "text","text", "text", "text", "text", "text","text", "text"))
WA.EMAP <- read.csv(paste0(PATH_occ, "/Other_Datasets/fishDataWA_EMAPandNRSAandEcology.csv"))
conner.mus <- read_excel(paste0(PATH_occ, "/Other_Datasets/ConnerMuseum_WSU_FishCollection.xlsx"), sheet = "Catalog")
hb <- read.csv(paste0(PATH_occ, "/Other_Datasets/HB_1990_2019.csv"))
hb.key <- read.csv(paste0(PATH_occ, "/Other_Datasets/master_fish_species_Harneybasin.csv"))

#prep datasets
conner.mus <- conner.mus %>% 
  select(c(3,41,44,55,59:62)) %>%
  filter(COUNTRY == "United States") %>%
  mutate(dataset_ID = "ConnerMuseum") %>%
  rename(decimalLatitude = DecimalLatitude, 
         decimalLongitude = DecimalLongitude, 
         orig.species = TAXON_NAME_CURRENT, 
         year = Cat_year, 
         coordinateUncertaintyInMeters = CoordinateUncertaintyInMeters,
         family = Family,
         COUNTRY = countryCode) %>%
  filter(!is.na(decimalLatitude)) 
conner.mus[ conner.mus == "United States" ] <- "US" #fix common name mislabel

WA.EMAP <- WA.EMAP %>%
  select(c(4,5,9,11,12)) %>%
  mutate(dataset_ID = "WA_EMAP") %>%
  rename(decimalLatitude = Latitude,
         decimalLongitude = Longitude,
         year = Year,
         orig.species = Scientific_Name,
         individualCount = Count) %>%
  filter(!is.na(decimalLatitude))

vertnet18 <- vertnet18 %>%
  select(c(3,4,6:8)) %>%
  mutate(dataset_ID = "vertnet18") %>%
  rename(orig.species = GenusSpeci,
         decimalLatitude = decimallat,
         decimalLongitude = decimallon,
         coordinateUncertaintyInMeters = coordinate) %>%
  filter(!is.na(decimalLatitude))

columbia.basin <- columbia.basin %>%
  select(c(3,5,6,9,10)) %>%
  mutate(dataset_ID = "ColumbiaBasin") %>%
  rename(decimalLatitude = LAT,
         decimalLongitude = LON,
         orig.species = "Scientific Name",
         individualCount = Count) %>%
  mutate(year = str_extract(Date, "[0-9]{4}")) %>%
  select(-one_of("Date")) %>%
  filter(!is.na(decimalLatitude))
columbia.basin[ columbia.basin == "Cutthroat Trout" ] <- "Oncorhynchus clarkii" #fix common name mislabel
columbia.basin[ columbia.basin == "Paiute Sculpin" ] <- "Cottus beldingii" #fix common name mislabel

hb <- hb %>% 
  select(c(7,11:13,32,36:80)) %>%
  mutate(dataset_ID = "HarneyBasin") %>%
  filter(!is.na(Latitude)) %>% 
  pivot_longer(cols = AM:YP, names_to = "Code", values_to = "individualCount") %>%
  filter(individualCount != 0) %>%
  left_join(., hb.key) %>%
  rename(decimalLatitude = Latitude,
         decimalLongitude = Longitude,
         year = Year,
         elevation = Elev,
         orig.species = Species.scientific.name) %>%
  select(-c("Code", "Species.common.name")) %>%
  mutate(orig.species = gsub("\\s+", " ", orig.species))
  
#get species keys for matching names (since datasets have some non-focal species and to deal with possible synonyms)
##first find the unique specieskeys for species on the original list
pnw.gbif.taxa <- NULL
for(u in 1:nrow(species.list)){
  pnw.gbif.taxa1 <- name_backbone(name=species.list[u,2], 
                                  rank='species', kingdom='animals')
  pnw.gbif.taxa <- bind_rows(pnw.gbif.taxa1, pnw.gbif.taxa)
}
pnw.gbif.taxa <- pnw.gbif.taxa %>% filter(!duplicated(speciesKey))
  
comb.dat <- list(columbia.basin, conner.mus, vertnet18, WA.EMAP, hb) #list of occ data

##find species keys across all datasets, reduce datasets to only desired species
comb.dat.taxa <- lapply(comb.dat, function (x) {
  
  dat.species <- as.data.frame(unique(x$orig.species))
  taxa.key <- NULL
  
  for(i in 1:nrow(dat.species)){
    taxa.key1 <- name_backbone(name=dat.species[i,1], 
                                    rank='species', kingdom='animals')
    taxa.key <- bind_rows(taxa.key, taxa.key1)
  }
  
  taxa.key$orig.species <- dat.species[,1]
  x$speciesKey <- taxa.key$speciesKey[ match(x$orig.species, taxa.key$orig.species) ]
  
  x <- x[!is.na(x$speciesKey),] #removing records w/ no species key (this includes records only w/ family, genus, a hybrid)
  x$species <- taxa.key$species[ match(x$speciesKey, taxa.key$speciesKey) ] #add corrected species names to dataset
  
  x <- x %>% filter(speciesKey %in% pnw.gbif.taxa$speciesKey) #remove non-focal species
  x$year <- as.numeric(x$year)
  
  return(x)
  
} )

comb.no.GBIF <- comb.dat.taxa %>% reduce(full_join) #combine all non-GBIF datasets
write.csv(comb.no.GBIF, file = paste0(PATH_occ, "/CombinedOccDat_noGBIF_", date, ".csv")) #save step


##+ combine GBIF + all other data ----

pnw.gbif$dataset_ID <- "GBIF"
all.dat <- full_join(comb.no.GBIF, pnw.gbif)
all.dat <- all.dat[, c(5,8,7,1,2,6,11,4,3,10,9,12:53) ] #rearranging... 

#count occurrences in each dataset
dataset_summ <- all.dat %>% 
  group_by(species) %>%
  mutate(combined_ct = n()) %>%
  mutate(GBIF_ct = sum(dataset_ID == "GBIF")) %>%
  mutate(ColumbiaBasin_ct = sum(dataset_ID == "ColumbiaBasin")) %>%
  mutate(ConnerMuseum_ct = sum(dataset_ID == "ConnerMuseum")) %>%
  mutate(VertNet18_ct = sum(dataset_ID == "vertnet18")) %>%
  mutate(WA_EMAP_ct = sum(dataset_ID == "WA_EMAP")) %>%
  mutate(HarneyBasin_ct = sum(dataset_ID == "HarneyBasin")) %>%
  slice(1) %>% select(species, combined_ct, GBIF_ct, ColumbiaBasin_ct, ConnerMuseum_ct, VertNet18_ct, WA_EMAP_ct, HarneyBasin_ct)
write.csv(dataset_summ, file = paste0(PATH_occ, "/AllData_Occ_Summary_", date, ".csv"))


##+ filter all data together ----

###++ remove duplicates of lat long year ----
full_count <- all.dat %>%
  group_by(species) %>%
  mutate(combined_ct = n()) %>%
  distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE) %>%
  mutate(lat_long_norep = n()) %>%
  slice(1) %>% select(species, combined_ct, lat_long_norep)
all.dat <- all.dat %>% group_by(species) %>% distinct(decimalLatitude, decimalLongitude, year, .keep_all = TRUE)


###++ pnw extent only ----
#make data spatial
dat_sp <- st_as_sf(all.dat, coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo)
dat_alb <- dat_sp %>% st_transform(crs.albers) #N. American projected crs

#clip data within 10 km buffered PNW extent
occ <- dat_alb[pnw.buff,]  

#save counts
pnw_ct <- occ %>% group_by(species) %>% count(.) 
full_count <- merge(full_count, pnw_ct, by="species", all.x=TRUE)
colnames(full_count)[4] <- "pnw_only"

#take out hybrids!
occ <- occ %>% filter((orig.species != "Oncorhynchus mykiss/clarkii") %>% replace_na(TRUE))

#save filtered occurrences only in PNW extent
occ <- st_transform(occ, crs.geo)
st_write(occ, dsn = paste0(PATH_occ, "/PNWextent_ALLocc_", date, ".csv"), delete_layer = T, layer_options = "GEOMETRY=AS_XY")
occ <- read.csv(paste0(PATH_occ, "/PNWextent_ALLocc_2021-08-19.csv"))

###++ year filter: 1995 - 2015 ----
#save counts
ct_yr <- occ %>% 
  group_by(species) %>%
  filter(year >= 1995 & year <= 2015) %>% count(.)
full_count <- merge(full_count, ct_yr, by="species", all.x=TRUE)
colnames(full_count)[5] <- "ct_1995_2015"

#filter
occ <- occ %>% group_by(species) %>% filter(year >= 1995 & year <= 2015)

#save dataset
write.csv(occ, file = paste0(PATH_occ, "/1995_2015_ALLocc_", date, ".csv"), row.names = FALSE)
  
###++ year filter: 2004 - 2015 ----
#save counts
ct_yr <- occ %>% 
  group_by(species) %>%
  filter(year >= 2004 & year <= 2015) %>% count(.)
full_count <- merge(full_count, ct_yr, by="species", all.x=TRUE)
colnames(full_count)[6] <- "ct_2004_2015"

#filter
occ <- occ %>% group_by(species) %>% filter(year >= 2004 & year <= 2015)

#save dataset
write.csv(occ, file = paste0(PATH_occ, "/2004_2015_ALLocc_", date, ".csv"), row.names = FALSE)

#save filter count
full_count[is.na(full_count)] <- 0
write.csv(full_count, file = paste0(PATH, "/AllDat_filter_tracksheet_", date, ".csv"))


# visualize ----

occ <- read.csv(paste0(PATH_occ, "/PNWextent_ALLocc_2021-08-19.csv"), 
                colClasses=c("character","character",rep("numeric",6), rep("character",4),"numeric",rep("character",16),rep("numeric",5),rep("character", 18)))

par(mfrow=c(1,2))
hist(p95$ct_1995_2015, breaks = 10, xlab = "num. of records per species", main = "Combined record counts (1995-2015)", ylim = c(0,20))
hist(p04$ct_2004_2015, breaks = 10, xlab = "num. of records per species", main = "Combined record counts (2004-2015)")

latlong <- occ %>% select(c("species","decimalLatitude","decimalLongitude"))
latlong <- latlong[order(latlong$species),]
write.csv(latlong, file = paste0(PATH_occ, "/NWCASC_SpeciesCoordinates_", date, ".csv"), row.names = FALSE)
