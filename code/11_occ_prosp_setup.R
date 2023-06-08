##TITLE: Occurrence Data Set Up + QAQC
##AUTHOR: C. E. Moore
##Updated on 22 AUG 2022

#PURPOSE: Format and QAQC occurrence data for use in RCS calculation

#set up ----

library(tidyverse); library(sf); library(terra); library(exactextractr) #local

#paths
PATH <- "~/Documents/Projects/PNW_fishes" #local
# PATH <- "/home/chloe9mo/PNW_fishes" #ARC

PATH_prosp <- paste0(PATH, "/PROSPER")
PATH_nhd <- paste0(PATH, "/NHD")
PATH_occdat <- paste0(PATH, "/occurrence_data")

#spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

#read in combined nhd+prosper data
prosp.nhd <- st_read(dsn=paste0(PATH_prosp,"/combined_streamdata"), layer = "prosper_nhd_flowline") #%>% 
  st_transform(crs.albers)

#snap occ to flowlines ----

#read in occurrence data
occ.dat <- read.csv(paste0(PATH_occdat, "/1995_2015_ALLocc_2021-08-19.csv"))
#remove undesired salmon species
occ.dat <- subset(occ.dat, !(species %in% c("Oncorhynchus mykiss", "Oncorhynchus gorbuscha", "Oncorhynchus kisutch", "Oncorhynchus tshawytscha")))
dat_sp <- st_as_sf(occ.dat, coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo)
dat_alb <- dat_sp %>% st_transform(crs.albers)

#assign COMID to each point + calc distance to nearest line
nr.line <- st_nearest_feature(dat_alb, prosp.nhd, check_crs = TRUE)
nhd_near <- prosp.nhd[nr.line,]
occ.dat$COMID <- nhd_near$COMID
dist <- as.vector(st_distance(dat_alb, nhd_near, by_element = TRUE))
occ.dat$dist2strm_km <- dist

#save w/ COMID and distance added to occurrence data
write.csv(occ.dat, file = paste0(PATH_occdat, "/NHD_near_distance_occurrences_", Sys.Date(), ".csv"))

#apply expert QAQC ----
occ.dat <- read.csv(paste0(PATH_occdat, "/NHD_near_distance_occurrences_2021-11-02.csv")) 

#here maps (created with 13a_HTML_OccurrenceMaps.Rmd & 13b_PDF_OccurrenceMaps.Rmd) were sent to PNW fisheries experts for comment
#1X_occ_expert_QAQC.R was used to consolidate comments w/ help from VJ Catalan
#note: this is from 1995 to 2015 only, fitting PROSPER limits. However, max. distance from stream has not been applied yet

#read in comment actions
act <- read.csv(paste0(PATH_occdat, "/QAQC_CommentActions.csv"))

#start with the easy ones - removal of pts
r <- act %>% filter(Action == "remove")
occ.dat <- occ.dat %>% filter(!X %in% r$Occurrence.ID)

#change species
r <- act %>% filter(Action == "change to C. perplexus")
occ.dat[occ.dat$X %in% r$Occurrence.ID,]$species <- "Cottus perplexus"
r <- act %>% filter(Action == "change to Catostomus tsiltcoosensis")
occ.dat[occ.dat$X %in% r$Occurrence.ID,]$species <- "Catostomus tsiltcoosensis"
r <- act %>% filter(Action == "change to P. umpquae")
occ.dat[occ.dat$X %in% r$Occurrence.ID,]$species <- "Ptychocheilus umpquae"
r <- act %>% filter(Action == "Change to R. evermanni")
occ.dat[occ.dat$X %in% r$Occurrence.ID,]$species <- "Rhinichthys evermanni"
#remove S. confluentus
occ.dat <- occ.dat %>% filter(species != "Salvelinus confluentus")

#remove records based on lat/long
r <- act %>% filter(Action == "in ID, remove pts <-113")
occ.dat <- occ.dat %>%
  filter(case_when(species == unique(r$Species) ~
                     !(between(decimalLongitude, -117, -113.4) & decimalLatitude < 42.5),
                   T ~ T))
r <- act %>% filter(Action == "remove > -120.89") #note: all of these pts removed by another step already
#actions: "remove all points >-121.0", "remove points > 120.2 lon" also covered in this step, same spp
occ.dat <- occ.dat %>%
  filter(case_when(species == unique(r$Species) ~
                     !(decimalLongitude > -120.89),
                   T ~ T))
r <- act %>% filter(Action == "remove all points < -120.07")
#actions: "remove all points < 47.245", "remove points < -119.2" also covered in this step
occ.dat <- occ.dat %>%
  filter(case_when(species == unique(r$Species) ~
                     !(decimalLongitude < -120.07 | decimalLatitude < 47.245),
                   T ~ T))
r <- act %>% filter(Action == "remove pts < 121.0")
#actions: "remove pts > -115.0 (long)" also for this spp
occ.dat <- occ.dat %>%
  filter(case_when(species == unique(r$Species) ~
                     !(decimalLongitude < -121 | decimalLongitude > -115),
                   T ~ T))
r <- act %>% filter(Action == "remove pts > -114.4, < 42.594")
occ.dat <- occ.dat %>%
  filter(case_when(species == unique(r$Species) ~
                     !((decimalLongitude > -114.4 & decimalLatitude < 42.59) | decimalLongitude > -113.7),
                   T ~ T))
r <- act %>% filter(Action == "remove pts < 43.8")
#actions: "remove pts > -113" also for this spp
occ.dat <- occ.dat %>%
  filter(case_when(species == unique(r$Species) ~
                     !((decimalLongitude > -113 & decimalLatitude < 43.8) | decimalLatitude < 42.5),
                   T ~ T))
#save step
write.csv(occ.dat, file=paste0(PATH_occdat, "/Occ_PostExpQAQC.csv"), row.names = F)

# QAQC climate outliers ----

occ.dat <- read.csv(paste0(PATH_occdat, "/filter_tracking_datasets/Occ_PostExpQAQC.csv"))
occ.dat <- occ.dat %>% filter(dist2strm_km < 0.5) #500m distance from prosper stream

#get prism data
pr_st <- rast(list.files(paste0(PATH,"/PRISM"), pattern = "asc.asc$", full.names = T))
huc12 <- readRDS(paste0(PATH,"/HUC/PNW_huc12.rds"))
pnw_prism <- crop(pr_st, huc12)
rm(huc12, pr_st)

#break apart occ data by species
occ.list <- split(occ.dat, f = occ.dat$species)
# occ.list <- occ.list[c(1:2)]
# x <- occ.list[[2]]

#for each species ...
out.occ.list <- lapply(occ.list, function(x) {
  
  sp.pt <- st_as_sf(x, coords=c("decimalLongitude","decimalLatitude"), crs=crs.geo) %>% st_transform(crs=crs(pnw_prism))
  sp.pri <- terra::extract(pnw_prism, vect(sp.pt), method="simple", list=F) #extract values of prism data at each pt
  colnames(sp.pri) <- c("cell", "ppt", "tmax", "tmin")
  x <- cbind(x, sp.pri)
  sp.nhd <- prosp.nhd %>% right_join(., x) %>% #get occupied streams
    rename(ID = X) %>% #retain ID for occurrence data
    select(ID, yr_MEANspc, yr_MEANspp, yr_m_AT, yr_mn_Q, ppt, tmax, tmin) %>% #get prosper data
    st_drop_geometry()
  
  if (nrow(sp.nhd) == 1) { #can't have an outlier if there's only one data pt (these will get removed anyway)
    
  out.occ <- x %>% rename(ID = X) %>%
    mutate(out_q = "FALSE", out_sd = "FALSE")
    
  } else {
  
  sp.nhd$out_q <- rowSums(apply(sp.nhd[c('yr_MEANspp', 'yr_m_AT', 'yr_mn_Q', 'ppt', 'tmax', 'tmin')], 2, function(y) #for each var. calc, pts falling outside #SD of mean and flag
    y < quantile(y, 0.003, na.rm=T) | y > quantile(y, 0.997, na.rm=T))) > 0 #quantile
  sp.nhd$out_sd <- rowSums(apply(sp.nhd[c('yr_MEANspp', 'yr_m_AT', 'yr_mn_Q', 'ppt', 'tmax', 'tmin')], 2, function(y) #for each var. calc, pts falling outside #SD of mean and flag
    y < mean(y, na.rm=T)-3*sd(y, na.rm = T) | y > mean(y, na.rm=T)+3*sd(y, na.rm=T) )) > 0 #quantile
  sp.nhd <- sp.nhd %>% select(-all_of(c('ppt','tmax','tmin'))) #remove duplicates
  out.occ <- x %>% rename(ID = X) %>% 
    left_join(., sp.nhd, by='ID') %>% 
    select(-all_of(c('yr_MEANspc', 'yr_MEANspp', 'yr_m_AT', 'yr_mn_Q'))) %>%
    mutate(out_q = as.character(out_q),
           out_sd = as.character(out_sd))
  
  }
  
  return(out.occ)
  
})

#save flagged and filtered occurrences
occ.dat.out <- bind_rows(out.occ.list) #make it one big dataframe again
write.csv(occ.dat.out, file = paste0(PATH_occdat, "/occ_500mdist_outliers_flagged.csv"), row.names = F)
occ.dat.out <- occ.dat.out %>% filter(out_sd == FALSE)
write.csv(occ.dat.out, file = paste0(PATH_occdat, "/occ_500mdist_no_outliers.csv"), row.names = F)

#visualize ----

##+ check distance code worked ----
# library(mapview)
# 
# #read in occurrence data
# occ.dat <- read.csv(paste0(PATH_occdat, "/NHD_near_distance_occurrences_2021-11-02.csv"))
# #remove undesired salmon species
# occ.dat <- subset(occ.dat, !(species %in% c("Oncorhynchus mykiss", "Oncorhynchus gorbuscha", "Oncorhynchus kisutch", "Oncorhynchus tshawytscha")))
# dat_sp <- st_as_sf(occ.dat, coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo)
# dat_alb <- dat_sp %>% st_transform(crs.albers)
# 
# #subsetting to make sure it worked
# ## look at specific comids
# nhd_sub <- nhd_near %>% filter(COMID == 23156084) %>% unique()
# dat_sub <- dat_alb %>% filter(COMID == 23156084)
# ## look at specific distances
# dat_sub <- dat_alb %>% filter(dist2strm_km < 20 & dist2strm_km > 10)
# nhd_sub <- nhd_near %>% filter(COMID %in% dat_sub$COMID) %>% unique()
# 
# mapview(nhd_sub$geometry) +
#   mapview(dat_sub$geometry)

##+ plot occurrences' >< 5km distances ----
# library(maps)
# pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon"))
# p <- ggplot() + 
#   geom_sf(data = pnw, fill = "gray") +
#   geom_sf(data = dat_sp, aes(fill = dist2strm_km), shape=21, size=1.5, alpha=0.7) +
#   scale_fill_viridis_b(trans = "sqrt") +
#   theme_bw() +
#   ggtitle("Distance to nearest NHD stream (km)") + xlab("Longitude") + ylab("Latitude")
# ggsave("Occurrence_Strm_Distance.png")


##+ check spp points during filtering steps ----
# library(mapview)
# occ.dat %>% 
#   # filter(species == unique(r$Species)) %>% 
#   filter(species == "Cottus klamanthus") %>%
#   st_as_sf(coords=c("decimalLongitude","decimalLatitude"), crs=crs.geo) %>%
#   mapview()


##+ check outlier QAQC ----
# out.occ <- read.csv(paste0(PATH_occdat, "/Allocc_NHDdist_outliers_2021-12-01.csv"))
# table((out.occ %>% filter(outliers == "TRUE"))$species) #which species have outlier values?
# 
# out.sp <- out.occ %>% filter(species == "Rhinichthys osculus")
# sp.nhd <- prosp.nhd %>% right_join(., out.sp) %>% #get occupied streams
#   select(ID, yr_MEANspp, yr_m_AT, yr_mn_Q, ppt, tmax, tmin, outliers) %>% 
#   st_drop_geometry() %>%
#   pivot_longer(cols = c(starts_with("yr"), "ppt", "tmax", "tmin"), names_to = "env_var", values_to = "env_val")
# 
# ggplot(sp.nhd %>% group_by(env_var), aes(x=factor(1), y=env_val)) + geom_boxplot() +
#   geom_jitter(aes(color = outliers), width = 0.1) +
#   facet_wrap(~env_var, scales = 'free')
# 
# v <- "yr_mn_Q"
# ggplot() + 
#   geom_density(data = sp.nhd %>% filter(env_var == v), aes(x=env_val), fill="dodgerblue", alpha=0.5) +
#   geom_vline(xintercept= (sp.nhd %>% filter(env_var == v & outliers == "TRUE"))$env_val, color="red") +
#   # facet_wrap(~env_var, scales='free') +
#   theme_bw() +
#   theme(legend.position = "none") + labs(x="value")
# ggplot() +
#   geom_boxplot(data = sp.nhd %>% filter(env_var == v), aes(x=env_val)) +
#   geom_point(data = sp.nhd %>% filter(env_var == v & outliers == "TRUE"), aes(x=env_val, y = 0, color = "red")) +
#   theme_bw() +
#   theme(legend.position = "none") + labs(x="value")
