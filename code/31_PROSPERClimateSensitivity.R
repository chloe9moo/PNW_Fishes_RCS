##TITLE: Climate Sensitivity Code -- PROSPER v.
##AUTHOR: C. E. Moore; modified from code by S. Silknetter + T. DuBose
##Updated on 22 AUG 2022

#PURPOSE: Using species occurrence records and environmental data, calculate the climate sensitivity (std dev of climate vars) for each species. 
##To be used to calculate RCS index with area of occupancy

#set up ----

#.libPaths(.libPaths()[3:1])
library(tidyverse); library(sf) #; library(raster); library(exactextractr)

#paths
PATH <- "~/Documents/Projects/PNW_fishes" #local
#PATH <- "/home/chloe9mo/PNW_fishes" #arc
# in
# PATH_HUC12 <- paste0(PATH, "/HUC/PNW_huc12.rds")
PATH_prosp <- paste0(PATH, "/PROSPER/combined_streamdata")
PATH_occ <- paste0(PATH, "/occurrence_data/occ_500mdist_no_outliers.csv")
PATH_huc_sp <- paste0(PATH, "/RCS_results/AOO/huc12_sp")
# out
PATH_CS_out <- paste0(PATH,"/RCS_results/CS")

#spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs") #contig us albers equal area

#read in data
#p <- readRDS(paste0(PATH_prosp, "/prosper_nhd_buffer.rds"))
nhd_p <- st_read(dsn=PATH_prosp, layer = "prosper_nhd_flowline") #%>% st_transform(crs.albers)

#get target species
species.list <- read.csv(paste0(PATH, "/occurrence_data/RCS_FocalSpecies.csv"))
species.list <- split(species.list$Species, seq(nrow(species.list)))


#watershed prosper CS calcs ----

#subset for testing
# x <- species.list[[3]]

all.huc.prosp <- lapply(species.list, function(x) {
  
  spp <- sub(" ", "_", x)
  huc.sp <- st_read(dsn = PATH_huc_sp, layer = spp) #call species HUC shp
  huc.sp <- st_transform(huc.sp, st_crs(nhd_p)) #match crs to prosper
  
  #extract values from streams overlapped by occupied hucs
  mn <- huc.sp %>%
    st_join(., nhd_p) %>% #add prosper values to overlapping watersheds
    #reduce to only unique streams in case some streams overlap multiple wtrsheds
    distinct(., COMID, .keep_all=TRUE) %>%
    filter(!is.na(yr_MEANspc)) %>%
    mutate(huc_mn_spc = mean(yr_MEANspc, na.rm=T),
           huc_mn_spp = mean(yr_MEANspp, na.rm=T),
           huc_mn_temp = mean(yr_m_AT, na.rm=T),
           huc_mn_Q = mean(yr_mn_Q, na.rm=T)) %>%
    st_drop_geometry() %>%
    dplyr::select(matches('huc_mn*')) %>%
    unique()
  
  sd <- huc.sp %>%
    st_join(., nhd_p) %>% #add prosper values to overlapping watersheds
    #reduce to only unique streams in case some streams overlap multiple wtrsheds
    distinct(., COMID, .keep_all=TRUE) %>%
    filter(!is.na(yr_MEANspc)) %>%
    mutate(huc_sd_spc = sd(yr_MEANspc, na.rm=T),
           huc_sd_spp = sd(yr_MEANspp, na.rm=T),
           huc_sd_temp = sd(yr_m_AT, na.rm=T),
           huc_sd_Q = sd(yr_mn_Q, na.rm=T)) %>%
    st_drop_geometry() %>%
    dplyr::select(matches('huc_sd*')) %>%
    unique()
  
  mn <- mn %>% 
    pivot_longer(cols=everything(), names_to = "value_type", values_to = "mean") %>%
    mutate(species=spp,
           value_type=case_when(grepl('huc_mn_spc', value_type) ~ 'strm_perm_class',
                                grepl('huc_mn_spp', value_type) ~ 'strm_perm_prob',
                                grepl('huc_mn_temp', value_type) ~ 'aug_temp',
                                grepl('huc_mn_Q', value_type) ~ 'baseflow'))
  sd <- sd %>% 
    pivot_longer(cols=everything(), names_to = "value_type", values_to = "sd") %>%
    mutate(species=spp,
           value_type=case_when(grepl('huc_sd_spc', value_type) ~ 'strm_perm_class',
                                grepl('huc_sd_spp', value_type) ~ 'strm_perm_prob',
                                grepl('huc_sd_temp', value_type) ~ 'aug_temp',
                                grepl('huc_sd_Q', value_type) ~ 'baseflow'))
  
  #combine + reorder
  huc.pr.val <- merge(mn, sd, by = c("species", "value_type"))
  
  return(huc.pr.val)
  
})

all.huc.prosp <- bind_rows(all.huc.prosp)
write.csv(all.huc.prosp, file = paste0(PATH_CS_out, "/CS_summ_huc_prosp_", Sys.Date(),".csv"), row.names = FALSE)

all.huc.prosp <- read.csv(paste0(PATH_CS_out, "/CS_summ_huc_prosp_2022-08-22.csv"))


#buffer prosper CS calcs ----

#read in occurrence data w/ nearest nhd stream comids
occ.dat <- read.csv(PATH_occ) %>% filter(species %in% cbind(species.list)) # using filtered data

#merge occurrence data w/ prosper by comid
prosp.occ <- occ.dat %>% 
  dplyr::select(-starts_with("yr_"), -contains("ppt"), -contains("tmax"), -contains("tmin")) %>%
  left_join(., nhd_p)

#find means
mn <- prosp.occ %>%
  group_by(species) %>%
  mutate(buf_mn_spc = mean(yr_MEANspc, na.rm=T),
         buf_mn_spp = mean(yr_MEANspp, na.rm=T),
         buf_mn_AT = mean(yr_m_AT, na.rm=T),
         buf_mn_Q = mean(yr_mn_Q, na.rm=T)) %>%
  dplyr::select(species, matches('buf_mn*')) %>%
  unique() %>% 
  pivot_longer(cols=contains('buf_mn'), names_to = "value_type", values_to = "mean") %>%
  mutate(value_type=case_when(grepl('buf_mn_spc', value_type) ~ 'strm_perm_class',
                              grepl('buf_mn_spp', value_type) ~ 'strm_perm_prob',
                              grepl('buf_mn_AT', value_type) ~ 'aug_temp',
                              grepl('buf_mn_Q', value_type) ~ 'baseflow'))


sd <- prosp.occ %>%
  group_by(species) %>%
  mutate(buf_sd_spc = sd(yr_MEANspc, na.rm=T),
         buf_sd_spp = sd(yr_MEANspp, na.rm=T),
         buf_sd_AT = sd(yr_m_AT, na.rm=T),
         buf_sd_Q = sd(yr_mn_Q, na.rm=T)) %>%
  dplyr::select(species, matches('buf_sd*')) %>%
  unique() %>% 
  pivot_longer(cols=contains('buf_sd'), names_to = "value_type", values_to = "sd") %>%
  mutate(value_type=case_when(grepl('buf_sd_spc', value_type) ~ 'strm_perm_class',
                              grepl('buf_sd_spp', value_type) ~ 'strm_perm_prob',
                              grepl('buf_sd_AT', value_type) ~ 'aug_temp',
                              grepl('buf_sd_Q', value_type) ~ 'baseflow'))

#combine + reorder
b.env.val <- merge(mn, sd, by = c("species", "value_type"))

write.csv(b.env.val, file = paste0(PATH_CS_out, "/CS_summ_buffer_prosp_", Sys.Date(),".csv"), row.names = FALSE)

#visuals ----

#visualize climate breadth across all taxa
# h <- b.env.val %>%
#   group_by(value_type) %>%
#   ggplot() +
#   geom_density(aes(x=mean, fill=value_type), alpha=.5) +
#   #geom_vline(aes(xintercept=mean(value), group=value_type))+
#   facet_wrap(~value_type, scales='free')+
#   theme_bw() + theme(#axis.text.x = element_text(angle=20),
#     axis.title.x=element_text(size=-1),
#     legend.position = 'none')
# ggsave(paste0(PATH_CS_out, "/buff_prosp_var_dist.jpg"), width=6, height=3)

# 
# #visualize occ pts and buffered NHD streams
# # library(ggplot2)
# b <- st_bbox(head(huc.sp, 1)) %>% st_as_sfc()
# c <- st_crop(buff.sp, b)
# h <- ggplot() +
#   geom_sf(data = c, fill = "grey", color = "black") +
#   geom_sf(data = nhd_p, color = "blue") +
#   coord_sf(xlim = st_coordinates(b)[c(1,2),1], # min & max of x values
#            ylim = st_coordinates(b)[c(2,3),2]) + # min & max of y values
#   theme_bw()
# ggsave(filename = paste0(PATH_prosp, "/nhd_buff_overlap.png"), plot = h)


