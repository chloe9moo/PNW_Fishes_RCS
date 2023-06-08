##TITLE: Climate Sensitivity Code -- PRISM
##AUTHOR: C. E. Moore; modified from code by S. Silknetter + T. DuBose
##Updated on 22 AUG 2022

#PURPOSE: Using species occurrence records and environmental data, calculate the climate sensitivity (std dev of climate vars) for each species. 
##To be used to calculate RCS index with area of occupancy

#set up ----

library(tidyverse); library(sf); library(terra); library(raster); library(exactextractr)
#note to self: https://tmieno2.github.io/R-as-GIS-for-Economists/extraction-speed-comparison.html#fnref85

#paths
PATH <- "~/Documents/Projects/PNW_fishes" #local
#PATH <- "/home/chloe9mo/PNW_fishes" #arc
# in
PATH_HUC12 <- paste0(PATH, "/HUC/PNW_huc12.rds")
PATH_env <- paste0(PATH, "/PRISM")
PATH_buff_sp <- paste0(PATH, "/RCS_results/AOO/buffer_sp")
PATH_huc_sp <- paste0(PATH, "/RCS_results/AOO/huc12_sp")
# out
PATH_CS_out <- paste0(PATH,"/RCS_results/CS")

#spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
#check the metadata - this is not exactly the same as what is given in the prosper metadata
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

#read in data, crop, remove
# climate -- note: reading in prism aug max temp, jan min temp, and annual precip
pr_st <- rast(list.files(PATH_env, pattern = "asc.asc$", full.names = T))
huc12 <- readRDS(PATH_HUC12)
pnw_prism <- crop(pr_st, huc12)
rm(huc12, pr_st)

#nhd_p <- st_read(dsn=PATH_prosp, layer = "prosper_nhd_flowline")

# aoo shapefiles
FILES_buffer <- list.files(path = PATH_buff_sp, pattern = ".shp")
FILES_huc <- list.files(path = PATH_huc_sp, pattern = ".shp")

#presence/absence matrix of spp in all watersheds ----

huc12.all <- readRDS(PATH_HUC12)
# st_crs(huc12.all)

hucIDlong <- NULL
for(i in 1:length(FILES_huc)){
  #read in the occupied watersheds
  shp <- st_read(dsn=PATH_huc_sp, layer=sub('.shp', '', FILES_huc[i]), quiet = T)
  
  #rearrange data to a long list of occupied watersheds
  hucIDlong1 <- data.frame(taxa=rep(sub(".shp", "", FILES_huc[i]), 
                                     length(shp$wtrshd_)),
                            watershed_ID=shp$wtrshd_) #%>%
  #mutate_if(is.factor, as.character)
  #add it all together into one long list
  hucIDlong <- rbind(hucIDlong, hucIDlong1)
}

#reducing huc spatial file for easier processing
huc12_red <- huc12.all %>% filter(huc12 %in% unique(hucIDlong$watershed_ID))
huc12_red <- huc12_red %>% rename(watershed.name = huc12) #make this easier for later

#convert this to a binary p/a matrix for easy multiplication later
hucIDbi <- hucIDlong %>%
  mutate(watershed_ID = as.character(paste(watershed_ID)),
         presence=1) %>%
  pivot_wider(names_from='taxa', values_from = 'presence', values_fill=list(presence = 0))
write.csv(hucIDbi, paste0(PATH_CS_out, "/HUC12_pa_matrix_", Sys.Date(), ".csv"), row.names = F)

#load back in, if necessary
hucIDbi <- read.csv(paste0(PATH_CS_out, "/HUC12_pa_matrix_2022-08-22.csv"))

#watershed climate mean calcs ----

#testing 
# hucs <- huc12_red[c(1:5),]

#extract the values of interest from the raster for each occupied polygon
##note: this function already weights by fraction of cell overlapped by polygon
r_pnw_prism <- raster::stack(pnw_prism)
huc12_red <- st_transform(huc12_red, crs = st_crs(pnw_prism))
env.values <- exact_extract(r_pnw_prism, huc12_red, fun = 'mean', force_df = TRUE, stack_apply = TRUE)

env.values$watershed.name <- huc12_red$watershed.name # add id column back to resulting dataframe
env.values <- env.values[,c(4,1:3)] #reorder resulting dataframe
  
write.csv(env.values, paste0(PATH_CS_out, "/cs_PRISMmn_all_huc12_", Sys.Date(), ".csv"), row.names = FALSE) #writing out a csv just in case
  
#check the watershed IDs still match (if character ->, often drop the leading 0)
# which(!(hucIDbi$watershed_ID %in% env.values$watershed.name))

#use matrix algebra to build a matrix to calculate standard of deviation
or.hucIDbi <- hucIDbi[order(hucIDbi$watershed_ID),] %>% dplyr::select(-watershed_ID)
or.env.values <- env.values[order(env.values$watershed.name),] %>% dplyr::select(-watershed.name)

all.mats <- sapply(or.env.values, function(x){
  
  cs.bimat <- x * or.hucIDbi
  cs.bimat <- cs.bimat %>%
    as_tibble() %>%
    bind_cols(watershed.name=env.values[order(env.values$watershed.name),1]) %>%
    pivot_longer(-watershed.name, names_to='species') %>%
    filter(value != 0) %>%
    mutate(watershed.name=as.character(watershed.name),
           value_origin='huc12')
  
  return(cs.bimat)
  
}, simplify = FALSE)

for (i in 1:length(all.mats)) { 
  
  all.mats[[i]] <- all.mats[[i]] %>% mutate(value_type=case_when(grepl('ppt', names(all.mats)[[i]]) ~ 'annual_ppt',
                                                                 grepl('tmax', names(all.mats)[[i]]) ~ 'aug_tmax',
                                                                 grepl('tmin', names(all.mats)[[i]]) ~ 'jan_tmin'))
  }

all.mats <- bind_rows(all.mats)
write.csv(all.mats, file = paste0(PATH_CS_out, "/CS_mn_prism_HUCvalues_", Sys.Date(), ".csv"), row.names = FALSE)

#testing it worked // w/ traci's version; note this is before the transformation starting at L112
# climate_df <- env.values[,c(1,2)]
# climate_df <- climate_df %>% filter(watershed.name %in% hucIDbi$watershed_ID)
# print(nrow(or.hucIDbi) == nrow(climate_df)) #check to make sure number of rows match
# climate_binmat <- climate_df[order(climate_df$watershed.name),2] * or.hucIDbi
# all.equal(climate_binmat, all.mats[[1]])

#summarize num. unique hucs occupied by species
View(
all.mats %>%
  group_by(species, value_origin) %>% tally() %>%
  mutate(n=n/3) %>%
  pivot_wider(names_from=value_origin, values_from=n)
)
#mean and standard deviation calculation
taxa.mn.sd.summ <- all.mats %>% 
  group_by(species, value_type) %>% 
  dplyr::summarise(mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
taxa.mn.sd.summ[is.na(taxa.mn.sd.summ)] <- 0

write.csv(taxa.mn.sd.summ, file = paste0(PATH_CS_out, "/CS_summ_huc_prism_", Sys.Date(), ".csv"))

#buffer climate mean calcs ----
all.buff <- lapply(FILES_buffer, function(x) {
  
  spp <- sub(".shp", "", x)
  buff.sp <- st_read(dsn = PATH_buff_sp, layer = spp) #call species buffer shp
  buff.sp <- st_transform(buff.sp, st_crs(pnw_prism)) #match crs to env. stack
  
  #extract values from env. layers overlapped by buffer shp
  ##note: this function already weights by fraction of cell overlapped by polygon
  mn <- exact_extract(r_pnw_prism, buff.sp, fun = 'mean', force_df = TRUE, stack_apply = TRUE)
  sd <- exact_extract(r_pnw_prism, buff.sp, fun = 'stdev', force_df = TRUE, stack_apply = TRUE)
  mn <- mn %>% 
    pivot_longer(cols=everything(), names_to = "value_type", values_to = "mean") %>%
    mutate(species=spp,
           value_type=case_when(grepl('ppt', value_type) ~ 'annual_ppt',
                                grepl('tmax', value_type) ~ 'aug_tmax',
                                grepl('tmin', value_type) ~ 'jan_tmin'))
  sd <- sd %>% 
    pivot_longer(cols=everything(), names_to = "value_type", values_to = "sd") %>%
    mutate(species=spp,
           value_type=case_when(grepl('ppt', value_type) ~ 'annual_ppt',
                                grepl('tmax', value_type) ~ 'aug_tmax',
                                grepl('tmin', value_type) ~ 'jan_tmin'))
  
  #combine + reorder
  b.env.val <- merge(mn, sd, by = c("species", "value_type"))
  
  return(b.env.val)
  
})

all.buff <- bind_rows(all.buff)
write.csv(all.buff, file = paste0(PATH_CS_out, "/CS_summ_buffer_prism_", Sys.Date(),".csv"), row.names = FALSE)

#summ. prism CS ----

# library(scales)
# all.buff$value_origin <- "buffer"
# taxa.mn.sd.summ$value_origin <- "huc12"
# cs <- all.buff
# 
# buff.summ <- all.buff %>%
#   group_by(species, value_type) %>%
#   #add ranks for quick comparison
#   mutate(rank_buff=rank(buff_area_sqkm),
#          rescale(aoo.out$mean_Rank, to = c(0,1))) %>%
#   rowwise() %>%
#   mutate(mean_rank=mean(c(rank_WS, rank_buff)),
#          sd_rank=sd(c(rank_WS, rank_buff)),
#          log10_buff_sqkm=log10(buff_area_sqkm),
#          log10_WS_sqkm=log10(huc_area_sqkm)) %>%
#   select(-c(huc_area_sqkm, buff_area_sqkm, total_pts, na_pts, na_ws))
# write.csv(AOOs_summarized, file = paste0(PATH_AOO_out, "/AOO_sumranks_", Sys.Date(), ".csv"), row.names = FALSE)


#visualize ----

##+ visualize climate breadth across all taxa ----
all.mats %>%
  group_by(value_type, value_origin) %>%
  ggplot() +
  geom_density(aes(x=value, fill=value_origin), alpha=.5) +
  #geom_vline(aes(xintercept=mean(value), group=value_type))+
  facet_wrap(~value_type, scales='free')+
  theme_bw() + theme(#axis.text.x = element_text(angle=20),
    axis.title.x=element_text(size=-1),
    legend.position = 'bottom')
# ggsave(paste0(PATH_CS_out, "/ws_climate_var_dens.jpg"), width=6, height=3)

##+ buffer vs. huc12 climate breadth density plot ----
ggplot() + geom_density(data = taxa.mn.sd.summ %>% group_by(value_type), aes(x=mean, fill="orange"), alpha=0.5) +
  geom_density(data = all.buff %>% group_by(value_type), aes(x=mean, fill="red"), alpha=0.5) +
  #geom_vline(aes(xintercept=mean(value), group=value_type))+
  facet_wrap(~value_type, scales='free') + 
  scale_fill_discrete(name = "value_origin", labels = c("huc12", "buffer")) +
  theme_bw() + theme(#axis.text.x = element_text(angle=20),
    axis.title.x=element_text(size=-1),
    legend.position = 'bottom')
# ggsave(paste0(PATH_CS_out, "/CS_distrib_prism.jpg"), width=6, height=3)

