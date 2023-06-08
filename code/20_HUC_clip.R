##TITLE: Clip HUC to Extent
##AUTHOR: C. E. Moore
##Updated on 30 AUG 2021

#PURPOSE: Clip HUC data to desired extent, save as RDS for read in to other scripts
#note this takes up ~13G of RAM

# packages
library(tidyverse); library(sf)

# paths
PATH <- "~/Documents/Projects/PNW_fishes"
#in
PATH_HUC12 <- paste0(PATH, "/HUC/WBD_National_GDB.gdb")
#out
PATH_HUCout <- paste0(PATH, "/HUC")

# crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")

# read in huc12 to national level
huc12_raw <- st_read(dsn = PATH_HUC12, layer = "WBDHU12") %>%
  st_transform(crs.geo) %>%
  rename(WSAREA='areasqkm') %>% #rename for convenient fn use later
  filter((states %in% c("OR","WA","ID","BC","ID,OR,WA", 
                        "ID,WA","ID,OR","ID,NV","ID,UT","ID,WY","ID,MT","ID,MT,WY","ID,UT,WY","ID,NV,UT", "ID,NV,OR",
                        "OR,WA", "CA,OR","NV,OR", "CA,NV,OR")))
         # , hutype != 'F') #this filter was in original RCS (anurans) for removing 'frontal' hydrologic units, but I think nature of fish makes this removal a problem

# save as rds
saveRDS(huc12_raw, file = paste0(PATH_HUCout, "/PNW_huc12.rds"))
st_write(huc12_raw %>% select(huc12, WSAREA), dsn = paste0(PATH,"/HUC/PNW_huc12.shp"), delete_layer = T)

