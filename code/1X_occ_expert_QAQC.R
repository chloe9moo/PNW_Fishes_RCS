#THIS SCRIPT WAS USED TO EASILY CHECK QAQC COMMENTS AND FIND PT IDS

#LOAD PACKAGES
#if installing:
#install.packages("tidyverse")
#install.packages("mapview")
#install.packages("sf")
library(tidyverse); library(mapview); library(sf)

#PATH TO FILES
PATH <- "~/Documents/Projects/PNW_fishes" #change path to match your computer
PATH_occdat <- paste0(PATH, "/occurrence_data")
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
huc12 <- readRDS(paste0(PATH, "/HUC/PNW_huc12.rds"))

#READ IN OCCURRENCE FILE
occ.dat <- read.csv(paste0(PATH_occdat, "/NHD_near_distance_occurrences_2021-11-02.csv"))

#SELECT SPECIES*** change only the words in quotes to match the species name
sp <- "Ptychocheilus oregonensis"

#FILTER OCCURRENCE DATA + MAKE IT SPATIAL
sp.occ <- occ.dat %>% 
  filter(species == sp) %>% 
  select(X, species, decimalLongitude, decimalLatitude) %>%
  st_as_sf(coords=c('decimalLongitude', 'decimalLatitude')) %>% st_set_crs(., crs.geo)

#MAP IT
# here record the number next to X in the google spreadsheet
mapview(sp.occ)
  # + mapview(st_as_sf(data.frame(long = -114.4011, lat = 42.5913), coords = c("long","lat")), col.regions="red")

# c.huc <- huc12 %>% select(huc12, name) %>%  #Cottus macroschelius
#   filter(grepl('[Uu]mpqua|[Ss]iuslaw|[Cc]oos|[Cc]oquille|[Ss]ixes', name)) %>%
#   st_union() %>%
#   st_buffer(dist=10000) #dist is meters
# c.huc <- huc12 %>% select(huc12, name) %>% 
#   filter(grepl('[Ll]ittle [Ss]almon|[Cc]lark [Ff]ork', name)) #Cottus cognatus Don Zaroban comment
# c.huc <- huc12 %>% select(huc12, name) %>% filter(grepl('[Ww]ood', name))
# c.huc <- huc12 %>% select(huc12, name) %>%  #Cottus rhotheus
#     filter(grepl('[Ss]almon [Rr]iver', name))
c.huc <- huc12 %>% select(huc12, name) %>%  #Ptychocheilus oregonensis
    filter(grepl('[Rr]ogue [Rr]iver', name))
# c.huc <- huc12 %>% select(huc12, name) %>% 
  # filter(grepl('[Pp]aulina|[Ee]ast [Ll]ake', name))
# c.f <- flowline %>% select(GNIS_NAME) %>%
  # filter(grepl('[Ss]nake [Rr]iver', GNIS_NAME)) %>%
  # st_union() %>%
  # st_buffer(5000) #5 km

mapview(r.pt, col.regions="red") + mapview(sp.occ) + mapview(c.huc, col.regions="red") + mapview(c.f)

r.pt <- sp.occ %>% 
  # st_filter(sf.pnw %>% filter(ID == "idaho")) %>%
  data.frame(st_coordinates(.)) %>%
  filter(Y < 47.245 | X.1 < -120.07) %>%
  st_as_sf()
r.pt <- sp.occ %>% st_filter(c.huc)
outside <- sp.occ %>% filter(!X %in% r.pt$X)
# outside <- sp.occ[!lengths(st_intersects(sp.occ, sf.pnw)), ] %>% data.frame(st_coordinates(.)) %>% 
#   filter(Y < 43.8) %>% 
#   st_as_sf()

#prep boundary + background
pnw <- maps::map('state', regions = c("Oregon", "Washington", "Idaho"), fill = TRUE, col="transparent") #load states of interest
IDs <- sapply(strsplit(pnw$names, ":"), function(x) x[1])
sf.pnw <- st_as_sf(pnw)
bb <- st_bbox(sf.pnw)

#hydrography load in
flowline <- st_read(paste0(PATH, "/NHD/"), layer = "PNW_highflow")
flowline <- flowline %>% st_set_crs(crs.geo) %>% st_make_valid(.) %>% st_intersection(., sf.pnw)
waterbody <- st_read(paste0(PATH,"/NHD/"), layer = "PNW_waterbody")
waterbody <- waterbody %>% st_set_crs(crs.geo) %>% st_make_valid(.) %>% st_intersection(., sf.pnw)

p2 <- ggplot() + xlim(c(bb[1], bb[3])) + ylim(c(bb[2], bb[4])) +
  geom_sf(data = sf.pnw, fill="gray89", color="black") +
  geom_sf(data = flowline, color="deepskyblue3", size=0.3) +
  geom_sf(data = waterbody, fill="deepskyblue3", color=NA) +
  geom_sf(data = sp.occ, fill=alpha("red", 0.5), shape=21, size=2) +
  # geom_point(data=data.frame(long = -114.4011, lat = 42.5913), aes(long, lat), color="red", size=3, shape=13) +
  # geom_sf(data = r.pt, fill=alpha("red", 0.5), shape=21, size=2) +
  # geom_sf(data = outside, fill=alpha("red", 0.5), shape=21, size=2) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black")) +
  labs(caption = paste0("Current records in red, total records =", nrow(sp.occ)))
  # labs(caption = "Records to be removed/changed in red")
ggsave(filename = paste0(PATH, "/occurrence_data/PDF_maps/", gsub(" ", "_", sp), "_mapcheck.png"), height = 4, width = 6)
