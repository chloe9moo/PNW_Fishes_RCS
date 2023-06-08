##TITLE: RCS Index
##AUTHOR: C. E. Moore; modified from code by S. Silknetter + T. DuBose
##Updated on 22 AUG 2022

#PURPOSE: Calculate CS from std dev. values, merge prior AOO calculations into an RCS output, 
# and calculate Rarity and Climate Sensitivity Index (RCS) for each species + range metric

#set up ----

library(tidyverse)

#paths
PATH <- "~/Documents/Projects/PNW_fishes"
# in
PATH_AOO <- paste0(PATH, "/RCS_results/AOO")
PATH_CS <- paste0(PATH, "/RCS_results/CS")
# out
PATH_RCS <- paste0(PATH, "/RCS_results/")

#focal species
species.list <- read.csv(paste0(PATH,"/occurrence_data/RCS_FocalSpecies.csv")) %>%
  mutate(species = gsub(" ", "_", Species))

#import aoo and cs outputs (PROSPER)
AOO <- read.csv(paste0(PATH_AOO, "/AOO_sumranks_2022-08-22.csv"))
raw.AOO <- read.csv(paste0(PATH_AOO, "/AOO_2022-08-22.csv"))
CS.huc <- read.csv(paste0(PATH_CS, "/CS_summ_huc_prosp_2022-08-22.csv"))
CS.buff <- read.csv(paste0(PATH_CS, "/CS_summ_buffer_prosp_2022-08-22.csv"))

#import aoo and cs outputs (prism)
# AOO is the same
CS.huc.prsm <- read.csv(paste0(PATH_CS, "/CS_summ_huc_prism_2022-08-22.csv"))
CS.buff.prsm <- read.csv(paste0(PATH_CS, "/CS_summ_buffer_prism_2022-08-23.csv"))

#occ dataset (outliers and distance already filtered)
occ.dat <- read.csv(paste0(PATH, "/occurrence_data/occ_500mdist_no_outliers.csv")) %>%
  filter(species %in% species.list$Species)
length(unique(occ.dat$species)) #there should be 29 unique species


  
#PROSPER RCS ----

##+ scale climate breadth 0 to 1 ----

#huc
wtrshd.cs.index <- CS.huc %>%
  mutate(species = str_replace(species, "_", " ")) %>%
  filter(species %in% species.list$Species) %>%
  group_by(value_type) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd-min(sd))/ #calculating the index
           (max(sd)-min(sd))) %>% 
  pivot_wider(names_from='value_type', values_from='var_CS_index')%>% #switching from long to wide dataframe
  dplyr::select(-mean, -sd) %>%
  group_by(species) %>%
  dplyr::summarise(spc_huc.CS=mean(strm_perm_class, na.rm=T),
                   spp_huc.CS=mean(strm_perm_prob, na.rm=T),
                   temp_huc.CS=mean(aug_temp, na.rm=T),
                   q_huc.CS=mean(baseflow, na.rm=T)) %>%
  rowwise() %>%
  mutate(spc_huc_CSind=mean(c(spc_huc.CS, temp_huc.CS, q_huc.CS)),
         spp_huc_CSind=mean(c(spp_huc.CS, temp_huc.CS, q_huc.CS)))

#buffer
buff.cs.index <- CS.buff %>%
  filter(species %in% species.list$Species) %>%
  group_by(value_type) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd-min(sd))/ #calculating the index
           (max(sd)-min(sd))) %>% 
  pivot_wider(names_from='value_type', values_from='var_CS_index')%>% #switching from long to wide dataframe
  dplyr::select(-mean, -sd) %>%
  group_by(species) %>%
  dplyr::summarise(spc_buf.CS=mean(strm_perm_class, na.rm=T),
                   spp_buf.CS=mean(strm_perm_prob, na.rm=T),
                   temp_buf.CS=mean(aug_temp, na.rm=T),
                   q_buf.CS=mean(baseflow, na.rm=T)) %>%
  rowwise() %>%
  mutate(spc_buff_CSind=mean(c(spc_buf.CS, temp_buf.CS, q_buf.CS)),
         spp_buff_CSind=mean(c(spp_buf.CS, temp_buf.CS, q_buf.CS))) #calculate the mean climate breadth

##+ scale aoo 0 to 1 ----

#using log transformed here
AOOs_index <- AOO %>% 
  mutate(species = scientific_name) %>% #for easy joining
  filter(species %in% species.list$Species) %>%
  mutate(WS_AOO_ind = (log10_WS_sqkm - min(log10_WS_sqkm))/
           (max(log10_WS_sqkm) - min(log10_WS_sqkm)),
         buff_AOO_ind = (log10_buff_sqkm - min(log10_buff_sqkm))/
           (max(log10_buff_sqkm) - min(log10_buff_sqkm))) %>%
  dplyr::select(-scientific_name, -mean_rank, -sd_rank)

##+ create rcs table ----

RCS <- buff.cs.index %>%
  left_join(wtrshd.cs.index, by='species') %>%
  left_join(AOOs_index, by='species') %>%
  # Subtract the scaled AOO + CS values from 1 
  # (so that low values = commonness, high values = vulnerability.
  mutate(AOO_WS_adj = 1 - WS_AOO_ind,
         AOO_BUF_adj = 1 - buff_AOO_ind,
         CS_spc_buf_adj = 1 - spc_buff_CSind,
         CS_spc_huc_adj = 1 - spc_huc_CSind,
         CS_spp_buf_adj = 1 - spp_buff_CSind,
         CS_spp_huc_adj = 1 - spp_huc_CSind,
         # Calculate Relative Climate Sensitivity Index by taking average of "1-AOO" and "1-CS".
         #  Where large values of RCS indicate species with small AOO and low range of climate variables.
         RCS_WS_spc = (AOO_WS_adj + CS_spc_huc_adj)/2,
         RCS_buff_spc = (AOO_BUF_adj + CS_spc_buf_adj)/2,
         RCS_WS_spp = (AOO_WS_adj + CS_spp_huc_adj)/2,
         RCS_buff_spp = (AOO_BUF_adj + CS_spp_buf_adj)/2) %>%
  arrange(RCS_buff_spp) %>%
  mutate(SpFac = factor(species, levels = species[order(RCS_buff_spp)])) %>%
  dplyr::select(-rank_WS, -rank_buff) %>%
  relocate(species, RCS_WS_spc, RCS_buff_spc, RCS_WS_spp, RCS_buff_spp, 
           spc_huc_CSind, spp_huc_CSind, spc_huc.CS, spp_huc.CS, temp_huc.CS, q_huc.CS, 
           spc_buff_CSind, spp_buff_CSind, spc_buf.CS, spp_buf.CS, temp_buf.CS, q_buf.CS, 
           WS_AOO_ind, log10_WS_sqkm, buff_AOO_ind, log10_buff_sqkm)

# Exports RCS data table.
write.csv(RCS, file = paste0(PATH_RCS, "/RCS_table_prosp_", Sys.Date(), ".csv"), row.names = F)

#PRISM RCS ----

##+ scale climate breadth 0 to 1 (WC) ----
# CS.raw %>% filter(value_type=='tmax') %>% group_by(species, value_origin) %>% count() #num watersheds per spp

#huc
wtrshd.cs.index.prsm <- CS.huc.prsm %>%
  dplyr::select(-X) %>%
  mutate(species = str_replace(species, "_", " ")) %>%
  filter(species %in% species.list$Species) %>%
  group_by(value_type) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd-min(sd))/ #calculating the index
           (max(sd)-min(sd))) %>% 
  pivot_wider(names_from='value_type', values_from='var_CS_index')%>% #switching from long to wide dataframe
  dplyr::select(-mean, -sd) %>%
  group_by(species) %>%
  dplyr::summarise(ppt_huc.CS=mean(annual_ppt, na.rm=T),
                   tmax_huc.CS=mean(aug_tmax, na.rm=T),
                   tmin_huc.CS=mean(jan_tmin, na.rm=T)) %>%
  rowwise() %>%
  mutate(prism_huc_CSind=mean(c(ppt_huc.CS, tmax_huc.CS, tmin_huc.CS))) #calc mean climate breadth -- prism

#buffer
buff.cs.index.prsm <- CS.buff.prsm %>%
  mutate(species = str_replace(species, "_", " ")) %>%
  filter(species %in% species.list$Species) %>%
  group_by(value_type) %>% #group_by to scale the climate variables appropriately
  mutate(var_CS_index = (sd-min(sd))/ #calculating the index
           (max(sd)-min(sd))) %>% 
  pivot_wider(names_from='value_type', values_from='var_CS_index')%>% #switching from long to wide dataframe
  dplyr::select(-mean, -sd) %>%
  group_by(species) %>%
  dplyr::summarise(ppt_buff.CS=mean(annual_ppt, na.rm=T),
                   tmax_buff.CS=mean(aug_tmax, na.rm=T),
                   tmin_buff.CS=mean(jan_tmin, na.rm=T)) %>%
  rowwise() %>%
  mutate(prism_buff_CSind=mean(c(ppt_buff.CS, tmax_buff.CS, tmin_buff.CS))) #calc mean climate breadth -- prism

##+ scale aoo 0 to 1 ----

#using log transformed here
AOOs_index <- AOO %>% 
  mutate(species = scientific_name) %>% #for easy joining
  filter(species %in% species.list$Species) %>%
  mutate(WS_AOO_ind = (log10_WS_sqkm - min(log10_WS_sqkm))/
           (max(log10_WS_sqkm) - min(log10_WS_sqkm)),
         buff_AOO_ind = (log10_buff_sqkm - min(log10_buff_sqkm))/
           (max(log10_buff_sqkm) - min(log10_buff_sqkm))) %>%
  dplyr::select(-scientific_name, -mean_rank, -sd_rank)

##+ create rcs table ----

RCS.prsm <- buff.cs.index.prsm %>%
  left_join(wtrshd.cs.index.prsm, by='species') %>%
  left_join(AOOs_index, by='species') %>%
  # Subtract the scaled AOO + CS values from 1 
  # (so that low values = commonness, high values = vulnerability.
  mutate(AOO_WS_adj = 1 - WS_AOO_ind,
         AOO_BUF_adj = 1 - buff_AOO_ind,
         CS_pr_buf_adj = 1 - prism_buff_CSind,
         CS_pr_huc_adj = 1 - prism_huc_CSind,
         # Calculate Relative Climate Sensitivity Index by taking average of "1-AOO" and "1-CS".
         #  Where large values of RCS indicate species with small AOO and low range of climate variables.
         RCS_WS_prsm = (AOO_WS_adj + CS_pr_huc_adj)/2,
         RCS_buff_prsm = (AOO_BUF_adj + CS_pr_buf_adj)/2) %>%
  arrange(RCS_buff_prsm) %>%
  mutate(SpFac = factor(species, levels = species[order(RCS_buff_prsm)])) %>%
  dplyr::select(-rank_WS, -rank_buff) %>%
  relocate(species, RCS_WS_prsm, RCS_buff_prsm, 
           prism_huc_CSind, ppt_huc.CS, tmax_huc.CS, tmin_huc.CS,
           prism_buff_CSind, ppt_buff.CS, tmax_buff.CS, tmin_buff.CS, 
           WS_AOO_ind, log10_WS_sqkm, buff_AOO_ind, log10_buff_sqkm)
#Exports RCS data table.
write.csv(RCS.prsm, file = paste0(PATH, "/RCS_results/RCS_table_prism_", Sys.Date(), ".csv"), row.names = F)

#combine CS into one table
cs.ind.all <- wtrshd.cs.index %>%
  dplyr::select(species, spc_huc_CSind, spp_huc_CSind) %>%
  left_join(wtrshd.cs.index.prsm %>% dplyr::select(species, prism_huc_CSind)) %>%
  left_join(buff.cs.index %>% dplyr::select(species, spc_buff_CSind, spp_buff_CSind)) %>%
  left_join(buff.cs.index.prsm %>% dplyr::select(species, prism_buff_CSind)) %>%
  mutate(across(contains("CSind"), ~1-.x)) # switch for 1 = small climate breadth
write.csv(cs.ind.all, file = paste0(PATH_RCS, "/CS/CS_table_combined.csv"), row.names = F)
