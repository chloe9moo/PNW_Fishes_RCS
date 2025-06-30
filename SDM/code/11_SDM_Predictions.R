## 11. Model Predictions
## Project: PNW Fishes SDMs
## Author: CE Moore
## Date: 28 SEP 2022

## Purpose: Take MAXENT models and predict to NHD flowlines within species' ranges

#set up ----

# .libPaths(.libPaths()[3:1]) #for ARC only
library(tidyverse); library(sf); library(ENMeval)

#home
PATH <- "~/Documents/Projects/PNW_fishes" #local
# PATH <- "/data" #ARC

#get species run
species.list <- read.csv(paste0(PATH, "/occurrence_data/RCS_FocalSpecies.csv")) %>% 
  filter(N_final >= 20) %>% #species w/ >10 pts & >25% of their range in the study extent
  mutate(Species = gsub(" ", "_", Species))
#get occurrence data
occ <- read.csv(paste0(PATH, "/occurrence_data/occurrence_data_finalfilter.csv")) %>%
  filter(species %in% gsub("_", " ", species.list$Species))
#make list obj of species
species.list <- split(species.list$Species, seq(nrow(species.list)))

#subset for testing
# species.list <- species.list[1:2]

#set background type to predict
b.t <- "SpCont"

#bring in all env. variable flowlines
env.var <- read.csv(paste0(PATH, "/SDM/Input_Data/EnvPredictor_Table.csv")) %>% select(Variable.Name) %>% filter(!Variable.Name %in% c("smmr_fl", "wntr_fl"))
env <- st_read(dsn = paste0(PATH, "/SDM/Input_Data"), layer = "prosper_nhd_flowline")
env.swd <- env %>% 
  st_drop_geometry() %>%
  mutate(x = as.integer(COMID),
         y = as.integer(COMID)) %>%
  select(x, y, all_of(env.var$Variable.Name)) %>% #select only the variables we were interested in
  relocate(x, y) %>%
  replace(is.na(.), 0)

library(maps)
sf.pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon"))

#create range map folder key ----
r.map.key <- data.frame(folder = seq(1:length(list.dirs(paste0(PATH, "/range_maps"), recursive = F))),
                        species = NA)

for(i in 1:nrow(r.map.key)) {
  
  shp <- st_read(dsn=paste0(PATH, "/range_maps/", i), layer = "data_0") 
  sp <- unique(shp$BINOMIAL) #pull out species name
  r.map.key[i, 2] <- gsub(" ", "_", sp)
  
}

#prediction over range restricted NHD flowlines ----
## to double check IUCN range maps vs. convex hull clipping
species.list <- data.frame(spp = unlist(species.list)) %>% filter(spp %in% r.map.key$species)
species.list[nrow(species.list)+1, "spp"] <- "Acrocheilus_alutaceus"
species.list <- split(species.list$spp, seq(nrow(species.list)))

lapply(species.list, function(x) {
  
  mx <- readRDS(file=paste0(PATH, "/SDM/Output/Model_RDS/", x, "_ENMeval_", b.t, ".rds")) #get model
  res <- mx@results #get results
  
  ##+ select top model to predict ----
  #(this is only necessary if you forgot to save the top model from 10 >:( )
  if (any(res$cbi.val.avg > 0.6, na.rm = T)) {
    
    #if cbi values are > 0.7, select model based on highest CBI, with tie breaker of highest AUC
    cat("selecting based on CBI \n")
    opt.mod <- res %>%
      filter(cbi.val.avg == max(cbi.val.avg, na.rm = T)) %>%
      filter(auc.val.avg == max(auc.val.avg, na.rm = T))
    
  } else {
    
    #if cbi values are low, the model probably needs to be reworked, but pick the highest AUC
    cat("selecting based on AUC only \n")
    opt.mod <- res %>%
      filter(auc.val.avg == max(auc.val.avg))
    
  }
  
  opt.mod <- opt.mod[1,] #in case there's still a tie, just take one of them
  
  top.mod <- mx@models[[opt.mod$tune.args]] #pull out the top model
  
  ##+ get or create range shp ----
  if (any(r.map.key$species == x)) { 
    
    cat("Using IUCN range map for", x, "\n") 
    folder <- r.map.key[r.map.key$species == x, 1]
    shp <- st_read(dsn=paste0(PATH, "/range_maps/", folder), layer = "data_0", quiet = T) %>%
      filter(LEGEND == "Extant (resident)") %>%
      st_transform(st_crs(env)) %>% #match crs
      st_union() %>%
      st_buffer(10000) %>% #wiggle room
      st_make_valid()
    
    shp2 <- occ %>% filter(species == gsub("_", " ", x)) %>% #get species' occ data
      dplyr::select(species, decimalLongitude, decimalLatitude) %>%
      st_as_sf(coords=c('decimalLongitude', 'decimalLatitude')) %>% #convert to sf obj
      st_set_crs(., st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>% #convert to sf pts
      st_transform(st_crs(env)) %>% #match crs
      st_union() %>%
      st_convex_hull() %>%
      st_buffer(50000) %>% #do 50 km buffer on the hull
      st_sf() %>% st_make_valid() %>% filter(!st_is_empty(.))
    
  } else { 
      
    if (any((r.map.key %>% mutate(species = gsub(".*_", "", species)))$species == gsub(".*_", "", x))) {
      
      rkey2 <- r.map.key %>% mutate(species = gsub(".*_", "", species))
      folder <- rkey2[rkey2$species == gsub(".*_", "", x), 1]
      cat("Using IUCN range map for", x, "using name under", r.map.key[folder, 2], "\n") 
      
      shp <- st_read(dsn=paste0(PATH, "/range_maps/", folder), layer = "data_0", quiet = T) %>%
        filter(LEGEND == "Extant (resident)") %>%
        st_transform(st_crs(env)) %>% #match crs
        st_union() %>%
        st_buffer(10000) %>% #wiggle room
        st_make_valid()
      shp2 <- occ %>% filter(species == gsub("_", " ", x)) %>% #get species' occ data
        dplyr::select(species, decimalLongitude, decimalLatitude) %>%
        st_as_sf(coords=c('decimalLongitude', 'decimalLatitude')) %>% #convert to sf obj
        st_set_crs(., st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>% #convert to sf pts
        st_transform(st_crs(env)) %>% #match crs
        st_union() %>%
        st_convex_hull() %>%
        st_buffer(50000) %>% #do 50 km buffer on the hull
        st_sf() %>% st_make_valid() %>% filter(!st_is_empty(.))
      
      } else { 
        
        cat("Using convex hull for", x, "\n")
        
        shp <- occ %>% filter(species == gsub("_", " ", x)) %>% #get species' occ data
          dplyr::select(species, decimalLongitude, decimalLatitude) %>%
          st_as_sf(coords=c('decimalLongitude', 'decimalLatitude')) %>% #convert to sf obj
          st_set_crs(., st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>% #convert to sf pts
          st_transform(st_crs(env)) %>% #match crs
          st_union() %>%
          st_convex_hull() %>%
          st_buffer(50000) %>% #do 50 km buffer on the hull
          st_sf() %>% st_make_valid() %>% filter(!st_is_empty(.))
        
        }
    
    }
  
  ##+ clip env to range shp ----
  sp.env <- sf::st_filter(env, shp) 
  sp.env.swd <- env.swd %>% filter(y %in% (sp.env %>% mutate(COMID = as.integer(COMID)))$COMID)
  
  sp.env2 <- sf::st_filter(env, shp2) 
  sp.env.swd2 <- env.swd %>% filter(y %in% (sp.env2 %>% mutate(COMID = as.integer(COMID)))$COMID)
  
  ##+ predict ----
  pred <- dismo::predict(top.mod, sp.env.swd)
  sp.env.swd$pred <- pred
  sp.n.pred <- left_join(sp.env, sp.env.swd %>% select(y, pred) %>% mutate(y=as.character(y)), by=c("COMID"="y"))
  
  pred2 <- dismo::predict(top.mod, sp.env.swd2)
  sp.env.swd2$pred_hull <- pred2
  sp.n.pred2 <- left_join(sp.env2, sp.env.swd2 %>% select(y, pred_hull) %>% mutate(y=as.character(y)), by=c("COMID"="y"))
  
  #save
  st_write(sp.n.pred, dsn = paste0(PATH, "/SDM/Output/Top_Models/", x, "_", b.t, "_iucn.shp"), delete_layer = T)
  st_write(sp.n.pred2, dsn = paste0(PATH, "/SDM/Output/Top_Models/", x, "_", b.t, "_hull.shp"), delete_layer = T)
  
  ##+ plot ---- 

  p1 <- ggplot() +
    # geom_sf(data = shp, color="black") +
    geom_sf(data = env, color="lightgrey", size=0.3) +
    geom_sf(data = sp.n.pred, aes(color=pred), size=0.3) +
    geom_sf(data = (occ %>%
                      filter(species == gsub("_", " ", x)) %>%
                      st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"), 
                               crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))), 
            color="black", pch=4, alpha=0.5) +
    # geom_sf(data = st_read(dsn=paste0(PATH, "/SDM/Input_Data/Bg_Points"), layer=paste0(x, "_", b.t, "_BG")), color="black", pch=4) +
    geom_sf(data = sf.pnw, color="black", fill=NA) +
    scale_color_gradient(low = "#d7d8d9", high = "#FF0000", breaks = c(0, 0.25, 0.50, 0.75, 1.00), name="Predicted \nsuitability") +
    labs(x="Long", y="Lat", title = paste0(gsub("_", " ", x), " IUCN Extent")) +
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))
  ggsave(paste0(PATH, "/SDM/Output/Top_Models/Pred_Maps/", x, "_", b.t, "_iucn.png"), plot=p1, height = 8, width = 8)
  
  p2 <- ggplot() +
    # geom_sf(data = shp, color="black") +
    geom_sf(data = env, color="lightgrey", size=0.3) +
    geom_sf(data = sp.n.pred2, aes(color=pred_hull), size=0.3) +
    geom_sf(data = (occ %>%
                      filter(species == gsub("_", " ", x)) %>%
                      st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"), 
                               crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))), 
            color="black", pch=4, alpha=0.5) +
    # geom_sf(data = st_read(dsn=paste0(PATH, "/SDM/Input_Data/Bg_Points"), layer=paste0(x, "_", b.t, "_BG")), color="black", pch=4) +
    geom_sf(data = sf.pnw, color="black", fill=NA) +
    scale_color_gradient(low = "#d7d8d9", high = "#FF0000", breaks = c(0, 0.25, 0.50, 0.75, 1.00), name="Predicted \nsuitability") +
    labs(x="Long", y="Lat", title = paste0(gsub("_", " ", x), " Convex Hull Extent")) +
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))
  ggsave(paste0(PATH, "/SDM/Output/Top_Models/Pred_Maps/", x, "_", b.t, "_hull.png"), plot=p2, height = 8, width = 8)

          } #end of function
  )

#visualize ----
# nhd.pred <- st_read(dsn = paste0(PATH, "/SDM/Output/Top_Models/Catostomus_ardens_SpCont.shp"))

# #example data for testing
# p <- read.csv(paste0(PATH, "/ExData/gadbis.csv")) %>% select(-sp)
# b <- read.csv(paste0(PATH, "/ExData/background.csv")) %>% select(-sp)
# e <- read.csv(paste0(PATH, "/ExData/rivers.csv"))

# evalplot.stats(e = mx, stats = "or.mtp", color = "fc", x.var = "rm")
# 
# crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")

#+ avg. pred of presence flowlines? ----
# a.p <- nhd.pred %>% filter(COMID %in% sp.occ$x)
# hist(a.p$pred)
# #bg flowlines?
# b.p <- nhd.pred %>% filter(COMID %in% bg$x)
# hist(b.p$pred)
# ggplot() +
#   geom_density(data = b.p, aes(x=pred), fill="#a0b0bb", color="#a0b0bb", alpha=0.5) +
#   geom_density(data = a.p, aes(x=pred), fill="#FF0000", color="#FF0000", alpha=0.5) +
#   labs(x=paste0("Predicted Habitat Suitability for ", spp), y="Density") +
#   scale_x_continuous(expand = c(0.01, 0.01)) +
#   scale_y_continuous(expand = c(0.01, 0.01)) +
#   theme(panel.background = element_blank(),
#         axis.line = element_line(color="black"),
#         panel.border = element_blank(),
#         panel.grid.major = element_line(color="grey"),
#         legend.position = "right")

#+ comparison of prediction flowlines ----
for (i in 1:length(species.list)) {
  
  iucn <- st_read(paste0(PATH, "/SDM/Output/Top_Models/", species.list[[i]], "_", b.t, "_iucn.shp")) %>% st_drop_geometry()
  hull <- st_read(paste0(PATH, "/SDM/Output/Top_Models/", species.list[[i]], "_", b.t, "_hull.shp")) %>% st_drop_geometry()
  combo <- full_join(iucn, hull) #combine
  env.c <- env %>% left_join(combo %>% select(COMID, pred, pred_hull)) %>% 
    replace_na(list(pred = 0, pred_hull = 0)) %>% mutate(p.diff = pred - pred_hull)
  
  p3 <- ggplot() +
    geom_sf(data = env.c, aes(color=p.diff), size=0.3) +
    geom_sf(data = (occ %>%
                      filter(species == gsub("_", " ", species.list[[i]])) %>%
                      st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"),
                               crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))),
            color="black", pch=4, alpha=0.5) +
    # geom_sf(data = st_read(dsn=paste0(PATH, "/SDM/Input_Data/Bg_Points"), layer=paste0(x, "_", b.t, "_BG")), color="black", pch=4) +
    geom_sf(data = sf.pnw, color="black", fill=NA) +
    scale_color_gradient2(low = "#1860a8", mid = "#D4D4D4", high = "#f03030", name="IUCN - C.Hull") +
    labs(x="Long", y="Lat", title = gsub("_", " ", species.list[[i]])) +
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))
  ggsave(paste0(PATH, "/SDM/Output/Top_Models/Pred_Maps/", species.list[[i]], "_IUCNxCHull_comp.png"), plot=p3, height = 8, width = 8)
  
}


#pulling out high % pred lines ----
#going with convex hull predictions as of SEP 28 2022 - see email w/ J. Dunham
r.map.key[r.map.key$species == "Gila_alutaceus",]$species <- "Acrocheilus_alutaceus" #updating mismatched names

lapply(species.list, function(x) {
  
  f.l <- list.files(paste0(PATH, "/SDM/Output/Top_Models/"), pattern = paste0(x, ".*shp"))
  if (length(f.l) > 1) { 
    shp <- st_read(paste0(PATH, "/SDM/Output/Top_Models/", x, "_", b.t, "_hull.shp")) %>%
      rename(any_of(c(pred = "pred_hull", pred = "pred_iucn")))
  } else { 
    shp <- st_read(paste0(PATH, "/SDM/Output/Top_Models/", x, "_", b.t, ".shp"))
    }

  shp8 <- shp %>% 
    filter(pred > 0.8) %>%
    mutate(thresh = ifelse(pred >= 0.9, "NinetyP", "EightyP"))
  
  p4 <- ggplot() +
    geom_sf(data = sf.pnw, color="#3D3D3D", fill= "black") +
    geom_sf(data = env %>% filter(StrmOrd > 2), color="#5C5C5C", size=0.2) +
    geom_sf(data = shp8, aes(color=thresh), size=0.7) +
    # geom_sf(data = (occ %>%
    #                   filter(species == gsub("_", " ", species.list[[i]])) %>%
    #                   st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"),
    #                            crs=st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84"))),
    #         color="black", pch=4, alpha=0.5) +
    # geom_sf(data = st_read(dsn=paste0(PATH, "/SDM/Input_Data/Bg_Points"), layer=paste0(x, "_", b.t, "_BG")), color="black", pch=4) +
    # geom_sf(data = sf.pnw, color="black", fill=NA) +
    scale_color_manual(values = c("NinetyP" = "#75DBCD", "EightyP" = "#F15025"), 
                       name="Predicted\nSuitability", 
                       labels = c("> 90%", "> 80%"),
                       guide = guide_legend(override.aes = list(size = 3))) +
    labs(x="Long", y="Lat", title = gsub("_", " ", x)) +
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), legend.key = element_blank())
  # p4
  ggsave(paste0(PATH, "/SDM/Output/Top_Models/Pred_Maps/Threshold/", x, "_predthresh.png"), plot=p4, height = 8, width = 8)
  
  sp.com <- occ %>% filter(species == gsub("_", " ", x)) %>% 
    mutate(COMID = as.character(COMID)) %>%
    select(COMID) %>% left_join(., shp %>% select(COMID, pred)) %>%
    mutate(grp = case_when(pred >= 0.9 ~ "NinetyP",
                           pred >= 0.8 & pred < 0.9 ~ "EightyP",
                           pred < 0.8 & pred >= 0.5 ~ "P50-80",
                           TRUE ~ "LessThan50")) %>%
    mutate(grp = factor(grp, levels=c("LessThan50", "P50-80", "EightyP", "NinetyP")))
  p5 <- ggplot(sp.com, aes(grp)) +
    geom_bar() +
    scale_x_discrete(labels=c("<50%", "50-80%", "80-90%", ">90%")) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = "", title = paste0("Predicted suitability for streams\nwith an occurrence record for ", gsub("_", " ", x))) +
    theme_classic() +
    theme(plot.title = element_text(size = 10))
  ggsave(paste0(PATH, "/SDM/Output/Top_Models/Pred_Maps/Threshold/", x, "_bars.png"), plot=p5, height = 4, width = 4)
  
  shp8 %>% st_drop_geometry() %>%
    mutate(species = x) %>% 
    relocate(species) %>%
    write_csv(., file=paste0(PATH, "/SDM/Output/Top_Models/Filtered_CSVs/", x, "_filtered_output.csv"))
  
})





PATH <- getwd()
PATH <- sub("/code", "", PATH)

sdm <- st_read(paste0(PATH, "/SDM/DataRelease/SDMs_HULLextent/Acrocheilus_alutaceus_hull.shp"))
