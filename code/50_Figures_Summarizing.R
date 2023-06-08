##TITLE: PNW Fishes Figures and Summarizing
##AUTHOR: C. E. Moore
##Updated on 22 AUG 2022

#PURPOSE: Calculate CS from std dev. values, merge prior AOO calculations into an RCS output, 
# and calculate Rarity and Climate Sensitivity Index (RCS) for each species + range metric

#set up ----

library(tidyverse); library(cowplot); library(scales); library(ggpubr); library(sf)

#paths
PATH <- "~/Documents/Projects/PNW_fishes"
PATH_AOO <- paste0(PATH, "/RCS_results/AOO")
PATH_CS <- paste0(PATH, "/RCS_results/CS")
PATH_RCS <- paste0(PATH, "/RCS_results/")

#focal species
species.list <- read.csv(paste0(PATH,"/occurrence_data/RCS_FocalSpecies.csv")) %>%
  mutate(species = gsub(" ", "_", Species))

#AOO
AOO <- read.csv(paste0(PATH_AOO, "/AOO_sumranks_2022-08-22.csv"))
AOOs_index <- AOO %>% 
  mutate(species = scientific_name) %>% #for easy joining
  filter(species %in% species.list$Species) %>%
  mutate(WS_AOO_ind = (log10_WS_sqkm - min(log10_WS_sqkm))/
           (max(log10_WS_sqkm) - min(log10_WS_sqkm)),
         buff_AOO_ind = (log10_buff_sqkm - min(log10_buff_sqkm))/
           (max(log10_buff_sqkm) - min(log10_buff_sqkm))) %>%
  dplyr::select(-scientific_name, -mean_rank, -sd_rank)
#CS
CS.huc <- read.csv(paste0(PATH_CS, "/CS_summ_huc_prosp_2022-08-22.csv"))
CS.buff <- read.csv(paste0(PATH_CS, "/CS_summ_buffer_prosp_2022-08-22.csv"))
CS.huc.prsm <- read.csv(paste0(PATH_CS, "/CS_summ_huc_prism_2022-08-22.csv"))
CS.buff.prsm <- read.csv(paste0(PATH_CS, "/CS_summ_buffer_prism_2022-08-23.csv"))
cs.ind.all <- read.csv(paste0(PATH_CS, "/CS_table_combined.csv"))
#RCS
RCS <- read.csv("~/Documents/Projects/PNW_fishes/RCS_results/RCS_table_prosp_2022-08-23.csv")
RCS.prsm <- read.csv("~/Documents/Projects/PNW_fishes/RCS_results/RCS_table_prism_2022-08-23.csv")

#occ dataset (outliers and distance already filtered)
occ.dat <- read.csv(paste0(PATH, "/occurrence_data/occ_500mdist_no_outliers.csv")) %>%
  filter(species %in% species.list$Species)
length(unique(occ.dat$species)) #there should be 29 unique species


#Figure A: Study Extent ----
library(maps)

# spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

#load in elements
pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon")) %>%
  st_transform(crs.geo)
usa <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% st_transform(crs.geo)
nhd_p <- st_read(dsn=paste0(PATH, "/PROSPER/combined_streamdata"), layer = "prosper_nhd_flowline")
sf_use_s2(FALSE)
nhd_pnw <- nhd_p %>%
  st_transform(crs.geo) %>%
  st_filter(., pnw) %>%
  st_crop(., pnw)
# huc12 <- readRDS(paste0(PATH, "/HUC/PNW_huc12.rds"))
dat_sp <- read.csv(paste0(PATH, "/occurrence_data/occ_500mdist_no_outliers.csv")) %>%
  filter(species %in% species.list$Species) %>%
  st_as_sf(coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo)

p <- ggplot() +
  # geom_sf(data = usa, fill = "#FAFAFA") +
  geom_sf(data = pnw, fill = "gray") +
  # geom_sf(data = huc12, fill="gray") +
  geom_sf(data = nhd_pnw %>% filter(LENGTHK > 0.1), color="blue", size = 0.12) +
  geom_sf(data = pnw, fill = NA, color="black") +
  # geom_sf(data = dat_sp, pch=21, fill="white", alpha=0.6, size=1) +
  coord_sf(xlim=c(-125,-110), ylim = c(41.5, 49.5)) +
  labs(x="Longitude", y="Latitude") +
  theme(rect = element_rect(fill="transparent"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
  
p
ggsave(filename=paste0(PATH_RCS, "/RCS_Study_Extent_nopts.png"), plot = p, width=4, height = 4)

#full N. America inset
na <- st_as_sf(map("world", regions = c("canada", "usa(?!:Hawaii)", "mexico", "guatemala", "belize", "el salvador", "honduras", "nicaragua", "costa rica", "panama"), plot = FALSE, fill = TRUE)) %>%
  st_transform(st_crs("+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"))

p2 <- ggplot() +
  geom_sf(data = na, fill="grey", color="black") +
  geom_sf(data = st_as_sfc(st_bbox(pnw)), fill = "red", color="red", size=2, alpha=0.5) +
  labs(x="Longitude", y="Latitude") +
  #coord_sf(xlim=c(-180,-40), ylim = c(15, 85)) +
  theme(rect = element_rect(fill="transparent"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
p2
ggsave(filename=paste0(PATH_RCS, "/RCS_Study_Extent_inset.png"), width = 4, height = 6)

sf_use_s2(TRUE)


#Figure B: Grain size comparisons ----
library(terra)
g.dat_sp <- st_crop(dat_sp, y = st_bbox(c(xmin=-122.6, xmax=-122.4, ymin=46.2, ymax=46.4)))
g.buff <- st_buffer(g.dat_sp, dist = 1000) %>% st_union()

g.nhd_pnw <- st_crop(nhd_pnw, y = st_bbox(c(xmin=-122.8, xmax=-122.2, ymin=46, ymax=46.5)))
g.nhd_pnw.p <- g.nhd_pnw %>% filter(COMID %in% g.dat_sp$COMID)

g.huc <- st_filter(huc12, g.nhd_pnw)
prez <- st_join(g.dat_sp, huc12)  %>% # Store the watershed name as an attribute of the data.
  rename(watershed.name=all_of('huc12')) # renaming watershed name column
g.huc.p <- g.huc %>% filter(huc12 %in% prez$watershed.name)
prez <- st_join(g.huc.p, g.nhd_pnw)
g.nhd.h <- g.nhd_pnw %>% filter(COMID %in% prez$COMID)

pr_st <- rast(list.files(paste0(PATH, "/PRISM"), pattern = "asc.asc$", full.names = T)[[1]])
pr_st <- pr_st %>% crop(ext(g.huc)) %>% as.data.frame(., xy=T)
colnames(pr_st) <- c("x", "y", "value")
# pr_st <- pr_st %>% filter(between(x, -122.6, -122.4)) %>% filter(between(y, 46.2, 46.4))

b1 <- ggplot() +
  geom_tile(data=pr_st, aes(x=x, y=y, fill=value)) + 
  #geom_sf(data = g.nhd_pnw, color="blue", size = 0.2) +
  geom_sf(data = g.huc, color = "black", size = 0.5, fill=NA, alpha=0.3) +
  geom_sf(data = g.huc.p, color = "black", size = 0.65, fill="black", alpha=0.2) +
  geom_sf(data = g.buff, color="red", size=0.5, fill="red", alpha=0.2) +
  geom_sf(data = g.dat_sp, pch=21, fill="white", size=2) +
  scale_fill_viridis_c(guide = "none") +
  labs(x="Longitude", y="Latitude", title="PRISM") +
  coord_sf(xlim = c(-122.73, -122.31), ylim = c(46.15, 46.45)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color="black", fill=NA))
b1
ggsave(filename=paste0(PATH_RCS, "/GrainComp_PRISM.png"), plot=b1, width=4, height=5)

b2 <- ggplot() +
  # geom_tile(data=pr_st, aes(x=x, y=y, fill=value)) + 
  geom_sf(data = g.huc.p, color = "black", size = 0.4, fill="black", alpha=0.2) +
  geom_sf(data = g.nhd_pnw, color="blue", size = 0.4, alpha = 0.3) +
  geom_sf(data = g.nhd.h, color="blue", size = 0.5) +
  geom_sf(data = g.nhd_pnw.p, color="red", size = 0.5) +
  geom_sf(data = g.huc, color = "black", size = 0.4, fill=NA, alpha=0.3) +
  #geom_sf(data = g.buff, color="red", size=0.5, fill="red", alpha=0.2) +
  geom_sf(data = g.dat_sp, pch=21, fill="white", size=2) +
  scale_fill_viridis_c(guide = "none") +
  labs(x="Longitude", y="Latitude", title="PROSPER") +
  coord_sf(xlim = c(-122.73, -122.31), ylim = c(46.15, 46.45)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid = element_blank())
b2
ggsave(filename=paste0(PATH_RCS, "/GrainComp_PROSPER.png"), plot=b2, width=4, height=5)

#Figure C: PROSP RCS Dotplot code  ----
rcs.theme <- list(
  scale_y_discrete(name=""),
  theme_cowplot(),
  facet_wrap(~facet.group, scales = 'free'),
  xlim(1.01, -0.01),
  labs(x = "RCS Index"),
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        plot.background = element_rect(fill = 'white'),
        #plot.margin = unit(c(1,1,15,1), "mm"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.title.y=element_text(size=-1),
        axis.text.x = element_text(size=9, angle=30, hjust=.9),
        axis.title.x = element_text(size=10, face='bold'),
        strip.background = element_blank(),
        strip.text=element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        #legend.position=c(-0.01,-0.1), #this only plots correctly in the saved figure
        legend.box="vertical",
        legend.title = element_text(size=10),
        legend.text = element_text(size=9),
        #legend.margin=margin(r=-20),
        legend.box.margin=margin(0,0,0,0))
)

# RCS <- read.csv("~/Documents/Projects/PNW_fishes/RCS_results/RCS_table_prosp_2022-06-10.csv")
r <- RCS %>%
  dplyr::select(species, SpFac, RCS_WS_spc, RCS_buff_spc, RCS_WS_spp, RCS_buff_spp) %>%
  arrange(RCS_buff_spp) %>%
  mutate(SpFac = factor(species, levels = species[order(RCS_buff_spp)])) %>%
  pivot_longer(cols=c('RCS_WS_spc','RCS_buff_spc','RCS_WS_spp', 'RCS_buff_spp')) %>%
  mutate(grain.size=recode(name, RCS_buff_spc='Snapped', RCS_buff_spp='Snapped',
                           RCS_WS_spc='Watersheds', RCS_WS_spp='Watersheds'),
         prosper=recode(name, RCS_buff_spc='Perm. Class', RCS_buff_spp='Perm. Prob.',
                        RCS_WS_spc='Perm. Class', RCS_WS_spp='Perm. Prob.'),
         filt.rank=as.numeric(SpFac),
         # box.min=ifelse(filt.rank > 20, 1.01, .97), #to make the split axis half and half
         # box.max=ifelse(filt.rank > 20, 1.03, .92),
         # boy.min=ifelse(filt.rank > 20, filt.rank-0.4-20, filt.rank-0.4),
         # boy.max=ifelse(filt.rank > 20, filt.rank+0.4-20, filt.rank+0.4),
         facet.group=factor(case_when(filt.rank > 14 ~'vulnerable',
                                      T~'not as vulnerable'),
                            levels=c('vulnerable','not as vulnerable'))) %>%
  #left_join(con_plot_df, by=c('species'='scientific_name')) %>%
  filter(!duplicated(.)) %>% #View()
  group_by(facet.group) %>% #summarize(max(value))
  ggplot() +
  #geom_rect(aes(xmin=box.min, xmax=box.max, ymin=boy.min, ymax=boy.max))+ #used for conservation status
  geom_point(aes(x=value, y=SpFac, shape=prosper, fill=grain.size), size=3.5, alpha=0.6)+
  scale_shape_manual('PROSPER Type:', values = rep(c(21,23),2))+
  #scale_color_manual('Grain Size', values = c('black','black'))+
  scale_fill_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .7)(2)))+
  #scale_x_reverse(name="RCS Index", expand=c(0,0.01), breaks=c(0, .25, .5, .75, .9, 1))+ #to separate x axis along the plots
  rcs.theme +
  guides(fill=guide_legend(override.aes=list(shape=21), nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")),
         shape=guide_legend(nrow=1, direction = "horizontal",title.theme = element_text(size=10,face="bold")))

ggsave(filename=paste0(PATH, "/RCS_results/RCS_graincomp_prosp_", Sys.Date(), ".png"), plot=r, width = 8, height = 7) 

# Figure D: PRISM RCS Dotplot code  ----
# RCS.prsm <- read.csv("~/Documents/Projects/PNW_fishes/RCS_results/RCS_table_prism_2022-06-10.csv")
rp <- RCS.prsm %>%
  dplyr::select(species, SpFac, RCS_WS_prsm, RCS_buff_prsm) %>%
  arrange(RCS_buff_prsm) %>%
  mutate(SpFac = factor(species, levels = species[order(RCS_buff_prsm)])) %>%
  pivot_longer(cols=c('RCS_WS_prsm','RCS_buff_prsm')) %>%
  mutate(grain.size=recode(name, RCS_buff_prsm='1 km Buffer', RCS_WS_prsm='Watersheds'),
         filt.rank=as.numeric(SpFac),
         # box.min=ifelse(filt.rank > 20, 1.01, .97),
         # box.max=ifelse(filt.rank > 20, 1.03, .92),
         # boy.min=ifelse(filt.rank > 20, filt.rank-0.4-20, filt.rank-0.4),
         # boy.max=ifelse(filt.rank > 20, filt.rank+0.4-20, filt.rank+0.4),
         facet.group=factor(case_when(filt.rank > 14 ~'vulnerable',
                                      T~'not as vulnerable'),
                            levels=c('vulnerable','not as vulnerable'))) %>%
  #left_join(con_plot_df, by=c('species'='scientific_name')) %>%
  filter(!duplicated(.)) %>% #View()
  group_by(facet.group) %>% #summarize(max(value))
  ggplot() +
  #geom_rect(aes(xmin=box.min, xmax=box.max, ymin=boy.min, ymax=boy.max))+ #used for conservation status
  geom_point(aes(x=value, y=SpFac, fill=grain.size), shape=21, size=3.5, alpha=0.6)+
  scale_fill_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .7)(2)))+
  #scale_x_reverse(name="RCS Index", expand=c(0,0.01), breaks=c(0, .25, .5, .75, .9, 1))+
  guides(fill=guide_legend(override.aes=list(shape=21), nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold"))) +
  rcs.theme


ggsave(filename=paste0(PATH, "/RCS_results/RCS_graincomp_prism_", Sys.Date(), ".png"), plot=rp, width = 8, height = 7) 

# Figure: AOO + CS Dotplot code  ----
dot.theme <- list(
  scale_y_discrete(name=""),
  #facet_wrap(~facet.group, scales = 'free'),
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        panel.grid.major.x = element_blank(),
        plot.background = element_rect(fill = 'white'),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        #plot.margin = unit(c(1,1,15,1), "mm"),
        axis.text.y = element_text(size=8, face='italic'),
        axis.title.y=element_text(size=-1),
        axis.text.x = element_text(size=9, angle=30, hjust=.9),
        axis.title.x = element_text(size=10, face='bold'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        #legend.position=c(-0.01,-0.1), #this only plots correctly in the saved figure
        legend.box="vertical",
        legend.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.margin=margin(-2,-2,-2,-2),
        legend.box.margin=margin(-5,-5,-5,-5),
        legend.key = element_blank()),
  guides(fill=guide_legend(override.aes=list(shape=21), nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")),
         shape=guide_legend(nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")))
)

a <- AOOs_index %>%
  dplyr::select(species, WS_AOO_ind, buff_AOO_ind) %>%
  mutate(WS_AOO_ind = 1 - WS_AOO_ind, buff_AOO_ind = 1 - buff_AOO_ind) %>% #switch for 1 = small aoo
  mutate(SpFac=factor(species, levels = species[order(WS_AOO_ind)])) %>%
  pivot_longer(cols=c('WS_AOO_ind','buff_AOO_ind')) %>%
  mutate(grain.size=recode(name, WS_AOO_ind='Watershed',
                           buff_AOO_ind='1 km buffer'),
         filt.rank=as.numeric(SpFac),
         # box.min=ifelse(filt.rank > 45, 1.01, .97),
         # box.max=ifelse(filt.rank > 45, 1.03, .92),
         # boy.min=ifelse(filt.rank > 45, filt.rank-0.4-45, filt.rank-0.4),
         # boy.max=ifelse(filt.rank > 45, filt.rank+0.4-45, filt.rank+0.4),
         shape.group=grain.size) %>%
  #left_join(con_plot_df, by=c('species'='scientific_name')) %>%
  filter(!duplicated(.)) %>% #View()
  #group_by(facet.group)%>% summarize(max(value))
  ggplot()+
  #geom_rect(aes(xmin=box.min, xmax=box.max, ymin=boy.min, ymax=boy.max))+ #used for conservation status
  geom_point(aes(x=value, y=SpFac, shape=shape.group, color=shape.group), size=2, alpha=0.6)+
  scale_shape_manual('Grain Size', values=rep(c(16,0),2))+
  scale_color_manual('Grain Size', values = c('black','black','grey','grey'))+
  scale_fill_manual(values=c(viridis_pal(option='magma', begin=.3, end=.75)(3), 'lightgrey')) +
  scale_x_reverse(name="AOO Index", expand=c(0.01,0.01), breaks=c(0, .25, .5, .75, 1)) +
  dot.theme
a
ggsave(filename=paste0(PATH, "/RCS_results/AOO_graincomp_", Sys.Date(), ".png"), plot=a, width = 6, height = 9) 

prsm.cs <- cs.ind.all %>%
  arrange(prism_huc_CSind) %>%
  mutate(SpFac=factor(species, levels = species[order(prism_huc_CSind)])) %>%
  pivot_longer(cols=names(cs.ind.all)[-1]) %>%
  mutate(grain.size = case_when(grepl("huc", name) ~ "Watershed",
                                grepl("buff", name) ~ "Buffer/Snapped"),
         env.var = case_when(grepl("spp", name) ~ "PROSPER Prob.",
                             grepl("spc", name) ~ "PROSPER Class",
                             grepl("prism", name) ~ "PRISM"),
         # filt.rank=as.numeric(SpFac)
  ) %>%
  filter(!duplicated(.)) %>%
  #group_by(facet.group)%>% summarize(max(value))
  filter(env.var == "PRISM") %>%
  ggplot() +
  #geom_rect(aes(xmin=box.min, xmax=box.max, ymin=boy.min, ymax=boy.max))+ #used for conservation status
  geom_point(aes(x=value, y=SpFac, shape=grain.size), size=3, alpha=0.6)+
  scale_shape_manual('Grain Size', values=rep(c(16,0),2))+
  scale_color_manual('Grain Size', values = c('black','black','grey','grey'))+
  scale_fill_manual(values=c(viridis_pal(option='magma', begin=.3, end=.75)(3), 'lightgrey')) +
  scale_x_reverse(name="CS Index", expand=c(0.01,0.01), breaks=c(0, .25, .5, .75, 1), limits = c(1,0)) +
  dot.theme
prsm.cs
ggsave(filename=paste0(PATH, "/RCS_results/CS_prism_graincomp_", Sys.Date(), ".png"), plot=prsm.cs, width = 6, height = 10) 

prsp.cs <- cs.ind.all %>%
  arrange(spp_huc_CSind) %>%
  mutate(SpFac=factor(species, levels = species[order(spp_huc_CSind)])) %>%
  pivot_longer(cols=names(cs.ind.all)[-1]) %>%
  mutate(grain.size = case_when(grepl("huc", name) ~ "Watershed",
                                grepl("buff", name) ~ "Buffer/Snapped"),
         env.var = case_when(grepl("spp", name) ~ "PROSPER Prob.",
                             grepl("spc", name) ~ "PROSPER Class",
                             grepl("prism", name) ~ "PRISM"),
         # filt.rank=as.numeric(SpFac)
  ) %>%
  filter(!duplicated(.)) %>%
  filter(env.var %in% c("PROSPER Prob.", "PROSPER Class")) %>%
  ggplot() +
  #geom_rect(aes(xmin=box.min, xmax=box.max, ymin=boy.min, ymax=boy.max))+ #used for conservation status
  geom_point(aes(x=value, y=SpFac, shape=grain.size, fill=env.var), size=3, alpha=0.6)+
  scale_shape_manual('Grain Size', values=rep(c(21,23),2))+
  scale_fill_manual('PROSPER Type:', values=c(viridis_pal(begin = .3, end = .7)(2)))+
  # scale_color_manual('PROSPER Type:', values = c('black','black','grey','grey'))+
  # scale_fill_manual(values=c(viridis_pal(option='magma', begin=.3, end=.75)(3), 'lightgrey')) +
  scale_x_reverse(name="CS Index", expand=c(0.01,0.01), breaks=c(0, .25, .5, .75, 1), limits = c(1,0)) +
  dot.theme
prsp.cs
ggsave(filename=paste0(PATH, "/RCS_results/CS_prosper_graincomp_", Sys.Date(), ".png"), plot=prsp.cs, width = 6, height = 10) 



#Figure G: Distribution of buff vs huc climate breadth vals ----
##+ PROSPER comparison ----
ggplot() + 
  geom_density(data = CS.huc %>% group_by(value_type), aes(x=sd, fill="#8ed9aa"), alpha=0.5) +
  geom_density(data = CS.buff %>% group_by(value_type), aes(x=sd, fill="#859fba"), alpha=0.5) +
  facet_wrap(~value_type, scales='free', labeller = as_labeller(c(aug_temp = "August Stream Temp.", baseflow = "Summer Baseflow", strm_perm_class = "Predicted Permanence Class", strm_perm_prob = "Predicted Permanence Prob."))) +
  scale_fill_discrete(type = c("#8ed9aa", "#859fba"), name = "Grain Size", labels = c("Watershed", "Snapped")) +
  theme_bw() + 
  labs(title = "PROSPER", y = "Standard Deviation Density") +
  theme(axis.title.x=element_text(size=-1),
        legend.position = 'bottom',
        plot.margin = margin(5, 10, 5, 5, unit = "pt"))
ggsave(paste0(PATH_RCS, "/CS_distrib_prosper.jpg"), width=6, height=5)


##+ PRISM comparison ----
ggplot() + 
  geom_density(data = CS.huc.prsm %>% group_by(value_type), aes(x=sd, fill="#8ed9aa"), alpha=0.5) +
  geom_density(data = CS.buff.prsm %>% group_by(value_type), aes(x=sd, fill="#859fba"), alpha=0.5) +
  facet_wrap(~value_type, scales='free', labeller = as_labeller(c(annual_ppt="Annual Precipitation", aug_tmax="August Max. Temp.", jan_tmin="January Min. Temp."))) +
  scale_fill_discrete(type = c("#8ed9aa", "#859fba"), name = "Grain Size", labels = c("Watershed", "1 km Buffer")) +
  theme_bw() + 
  labs(title = "PRISM", y = "Standard Deviation Density") +
  theme(axis.title.x=element_text(size=-1),
        legend.position = 'bottom')
ggsave(paste0(PATH_RCS, "/CS_distrib_prism.jpg"), width=6.5, height=4)


# Figure: Map of RCS values ----

# add RCS index
occ.dat <- occ.dat %>% left_join(., RCS, by="species") %>% select(species, decimalLatitude, decimalLongitude, dist2strm_km, RCS_WS_spp, RCS_buff_spp)

# spatial crs
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

dat_sp <- st_as_sf(occ.dat, coords=c("decimalLongitude","decimalLatitude"), crs = crs.geo)

#read in HUCs
huc12 <- readRDS(paste0(PATH, "/HUC/PNW_huc12.rds"))
occ.huc <- st_join(huc12, dat_sp, join = st_contains)
occ.huc <- occ.huc %>% 
  group_by(huc12) %>%
  mutate(avg.RCS.WS = mean(RCS_WS_spp),
         avg.RCS.buff = mean(RCS_buff_spp)) %>%
  select(huc12, avg.RCS.WS, avg.RCS.buff) %>%
  unique()

library(maps)
pnw <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>% filter(ID %in% c("washington", "idaho", "oregon"))
p <- ggplot() +
  geom_sf(data = pnw, fill = "gray") +
  geom_sf(data = occ.huc, aes(fill = avg.RCS.WS), lwd = 0) +
  # geom_sf(data = dat_sp, aes(fill = RCS_WS_spp), shape = 21, size = 2, alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Spatial distribution of RCS") + xlab("Longitude") + ylab("Latitude")
p
ggsave(filename="Avg_RCS_HUC12.png", plot = p)

# Figure F: n. occurence points vs. RCS / AOO----
# RCS.prsm <- read.csv(paste0(PATH_RCS, "/RCS_table_prism_2022-06-10.csv"))
# RCS <- read.csv(paste0(PATH_RCS, "/RCS_table_prosp_2022-06-10.csv"))

RCS.p <- RCS %>%
  dplyr::select(species, contains("RCS")) %>%
  left_join(RCS.prsm %>% dplyr::select(species, contains("RCS"))) %>%
  left_join(species.list %>% dplyr::select(Species, N_final), by=c("species"="Species"))

o.v.rcs <- RCS.p %>%
  pivot_longer(contains("RCS")) %>%
  mutate(grain.size = case_when(grepl("buff", name)~"Buffer/Snapped",
                                grepl("WS", name)~"Watershed"),
         env.var = case_when(grepl("spc", name)~"PROSPER Class",
                             grepl("spp", name)~"PROSPER Prob.",
                             grepl("prsm", name)~"PRISM")) %>%
  ggplot() +
  # geom_smooth(method = "gam", aes(x=N_final, y=value, color=grain.size), se=F) +
  geom_smooth(method="lm", formula = y~x, aes(x=N_final, y=value, color=grain.size), se=F) +
  geom_point(aes(x = N_final, y = value, shape=env.var, fill=grain.size), alpha=0.6) + 
  stat_cor(aes(x = N_final, y = value, color=grain.size), p.accuracy = 0.0001, size=4, show.legend = FALSE, label.x.npc = 0.6) + 
  scale_shape_manual('Climate Data:', values = c(22, 21, 23)) +
  scale_fill_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .8)(2))) +
  scale_color_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .8)(2))) +
  labs(x="Number of Occurrence Points", y="RCS Index") +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.line = element_line(color="black"),
        legend.key = element_blank()) +
  guides(fill=guide_legend(override.aes=list(shape=21), nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")),
         shape=guide_legend(nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")))
o.v.rcs

ggsave(filename = paste0(PATH_RCS,"/RCS_vs_occurrencePts.png"), plot = o.v.rcs, height = 6, width = 5)

AOO.p <- AOOs_index %>%
  select(species, contains("ind")) %>%
  mutate(across(contains("ind"), ~1-.x)) %>%
  left_join(species.list %>% select(Species, N_final), by=c("species"="Species")) %>%
  pivot_longer(contains("ind")) %>%
  mutate(grain.size = case_when(grepl("buff", name)~"Buffer",
                                grepl("WS", name)~"Watershed")) %>%
  ggplot() +
  geom_smooth(method = "lm", formula='y~x', aes(x=log(N_final), y=value, color=grain.size), se=F) +
  geom_point(aes(x = log(N_final), y = value, shape = grain.size, fill = grain.size), alpha = 0.6) +
  scale_shape_manual('Grain Size:', values = c(21,23)) +
  scale_fill_manual('Grain Size:', values = c("#D3D3D3", "#3D3D3D")) +
  scale_color_manual('Grain Size:', values = c("#D3D3D3", "#3D3D3D")) +
  labs(x="log(# of Occurrence Points)", y="AOO Index") +
  ylim(c(0,1)) +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.line = element_line(color="black"),
        legend.key = element_blank())
ggsave(filename = paste0(PATH_RCS,"/AOO_vs_occurrencePts.png"), plot = AOO.p, height = 6, width = 6)

# Figure: PROSPER boxplot ----

prob <- CS.buff %>% 
  filter(species %in% species.list$Species) %>%
  mutate(grain.size = "snapped") %>%
  select(-sd) %>%
  pivot_wider(names_from = value_type, values_from = mean) %>%
  select(species, grain.size, strm_perm_prob) %>%
  mutate(prob_bin = cut_number(strm_perm_prob, 4)) %>%
  mutate(prob_bin = case_when(prob_bin == "[0.323,0.554]" ~ "32% - 55%",
                              prob_bin == "(0.554,0.597]" ~ "56% - 59%",
                              prob_bin == "(0.597,0.653]" ~ "60% - 65%",
                              prob_bin == "(0.653,0.726]" ~ "66% - 73%"),
         prob_bin2 = case_when(strm_perm_prob < 0.45 ~ "< 45%",
                               strm_perm_prob >= 0.45 & strm_perm_prob < 0.55 ~ "45% - 55%",
                               strm_perm_prob >= 0.55 & strm_perm_prob < 0.65 ~ "55% - 65%",
                               strm_perm_prob >= 0.65 ~ " > 65%")) %>%
  mutate(prob_bin = as.factor(prob_bin),
         prob_bin2 = as.factor(prob_bin2)) %>%
  left_join(., RCS %>% select(species, RCS_buff_spp))
prob$prob_bin2 <- fct_relevel(prob$prob_bin2, "< 45%", "45% - 55%", "55% - 65%", " > 65%")
  
class <- CS.buff %>%
  filter(species %in% species.list$Species) %>%
  mutate(grain.size = "snapped") %>%
  select(-sd) %>%
  pivot_wider(names_from = value_type, values_from = mean) %>%
  select(species, grain.size, strm_perm_class) %>%
  mutate(class_bin = cut_number(strm_perm_class, 4)) %>%
  mutate(class_bin = case_when(class_bin == "[-3.44,0.672]" ~ "-3.00 - 0.67",
                              class_bin == "(0.672,1.17]" ~ "0.63 - 1.16",
                              class_bin == "(1.17,2.1]" ~ "1.17 - 2.09",
                              class_bin == "(2.1,3.99]" ~ "2.10 - 4.00"),
         class_bin2 = case_when(strm_perm_class < 0 ~ "< 0",
                                strm_perm_class >= 0 & strm_perm_class < 1 ~ "0 - 1",
                                strm_perm_class >= 1 & strm_perm_class < 2 ~ "1 - 2",
                                strm_perm_class >= 2 ~ "> 2")) %>%
  mutate(class_bin = as.factor(class_bin),
         class_bin2 = as.factor(class_bin2)) %>%
  left_join(., RCS %>% select(species, RCS_buff_spc))
  # rbind(., (CS.huc %>% mutate(grain.size = "Watershed")))
class$class_bin2 <- fct_relevel(class$class_bin2, "< 0", "0 - 1", "1 - 2", "> 2")

box.theme <- list(
  scale_fill_manual(values=c(viridis_pal(begin = .3, end = .7)(4))),
  ylim(c(0.0,1.0)),
  labs(x="Predicted Stream Permanence Probability", y="RCS Value"),
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        panel.grid.major = element_line(color = "#D3D3D3"),
        legend.position = "none")
)

p1 <- ggplot() +
  geom_boxplot(data=prob, aes(x=prob_bin, y=RCS_buff_spp, fill=prob_bin)) +
  geom_jitter(data=prob, aes(x=prob_bin, y=RCS_buff_spp), color="#808080", shape=16, position = position_jitter(0.2)) +
  box.theme
p2 <- ggplot() +
  geom_boxplot(data=prob, aes(x=prob_bin2, y=RCS_buff_spp, fill=prob_bin2)) +
  geom_jitter(data=prob, aes(x=prob_bin2, y=RCS_buff_spp), color="#808080", shape=16, position = position_jitter(0.2)) +
  box.theme
ggsave(filename = paste0(PATH_RCS, "/Perm_Prob_RCS_boxplot_equalgroups.png"), plot = p1, width = 5, height = 5)
ggsave(filename = paste0(PATH_RCS, "/Perm_Prob_RCS_boxplot_evenbreaks.png"), plot = p2, width = 5, height = 5)
## OR
p3 <- ggplot() +
  geom_boxplot(data=class, aes(x=class_bin, y=RCS_buff_spc, fill=class_bin)) +
  geom_jitter(data=class, aes(x=class_bin, y=RCS_buff_spc), color="#808080", shape=16, position = position_jitter(0.2)) +
  box.theme +
  xlab("Predicted Stream Permanence Class")
p4 <- ggplot() +
  geom_boxplot(data=class, aes(x=class_bin2, y=RCS_buff_spc, fill=class_bin2)) +
  geom_jitter(data=class, aes(x=class_bin2, y=RCS_buff_spc), color="#808080", shape=16, position = position_jitter(0.2)) +
  box.theme +
  xlab("Predicted Stream Permanence Class")
ggsave(filename = paste0(PATH_RCS, "/Perm_Class_RCS_boxplot_equalgroups.png"), plot = p3, width = 5, height = 5)
ggsave(filename = paste0(PATH_RCS, "/Perm_Class_RCS_boxplot_evenbreaks.png"), plot = p4, width = 5, height = 5)

# compare CS values at diff. grain sizes ----
cor.test(cs.ind.all$spc_buff_CSind, cs.ind.all$spp_buff_CSind, method="pearson") #0.925 rho, p = 6.95e-13
cor.test(cs.ind.all$spc_huc_CSind, cs.ind.all$spp_huc_CSind, method="pearson") #0.971 rho, p < 2.2e-16
cor.test(cs.ind.all$spc_buff_CSind, cs.ind.all$spc_huc_CSind, method="pearson") #0.833 rho, p = 2.08e-8
cor.test(cs.ind.all$spp_huc_CSind, cs.ind.all$spp_buff_CSind, method="pearson") #0.883 rho, p = 2.22e-10


#Figure E: compare prism vs prosp rcs ----

##+ find RCS deltas ----
RCS.p <- RCS.p %>%
  mutate(dClim_Prob_WS = (RCS_WS_prsm - RCS_WS_spp) / RCS_WS_spp * 100,
         dClim_Prob_B = (RCS_buff_prsm - RCS_buff_spp) / RCS_buff_spp * 100,
         mn_prob_d = (dClim_Prob_WS + dClim_Prob_B) / 2,
         dClim_Clas_WS = (RCS_WS_prsm - RCS_WS_spc) / RCS_WS_spc * 100,
         dClim_Clas_B = (RCS_buff_prsm - RCS_buff_spc) / RCS_buff_spc * 100,
         mn_Clas_d = (dClim_Clas_WS + dClim_Clas_B) / 2)

crcs <- RCS.p %>%
  pivot_longer(contains("RCS")) %>%
  mutate(grain.size = case_when(grepl("buff", name)~"Buffer/Snapped",
                                grepl("WS", name)~"Watershed"),
         env.var = case_when(grepl("spc", name)~"PROSPER Class",
                             grepl("spp", name)~"PROSPER Prob.",
                             grepl("prsm", name)~"PRISM")) %>%
  dplyr::select(-name, -N_final) %>%
  pivot_wider(names_from = env.var, values_from = value)


##+ prism vs. class ----
p.c1 <- ggplot() +
  geom_smooth(data=crcs, method = "lm", formula='y~x', aes(x = `PROSPER Class`, y = PRISM, color=grain.size), se=F) +
  #this can show change between huc and buffer:
  # geom_line(data=crcs %>% mutate(ad = abs(mn_Clas_d)) %>% slice_max(ad, n=10), 
  #           aes(x = `PROSPER Class`, y = PRISM, group = species), 
  #           linetype = "twodash", size=0.7, color = "red") +
  #change point size based on percent change
  geom_point(data=crcs %>% mutate(ad = abs(mn_Clas_d)), aes(x = `PROSPER Class`, y = PRISM, fill=grain.size, size=ad), shape=21, alpha=0.6) + 
  #highlight species with biggest difference between env. var types:
  # geom_point(data=crcs %>% mutate(ad = abs(mn_Clas_d)) %>% slice_max(ad, n=10), 
  #            aes(x = `PROSPER Class`, y = PRISM), 
  #            shape=1, color="red", size=2) +
  stat_cor(data=crcs, aes(x = `PROSPER Class`, y = PRISM, color=grain.size), p.accuracy = 0.0001, size=4, show.legend = FALSE) + 
  scale_fill_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .8)(2))) +
  scale_color_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .8)(2))) +
  scale_x_continuous(limits = c(0.15,1), breaks = c(0.25,0.50,0.75,1.00)) +
  scale_size_continuous('Percent Change:') +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.line = element_line(color="black"),
        legend.key = element_blank()) 
  # guides(fill=guide_legend(override.aes=list(shape=21), nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")),
  #        shape=guide_legend(nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")))
p.c1
ggsave(filename = paste0(PATH_RCS,"/RCSprospclass_vs_RCSprism.png"), plot = p.c1, height = 5.5, width = 5)


##+ prism vs. % ----
p.c2 <- ggplot() +
  geom_smooth(data=crcs, method = "lm", formula='y~x', aes(x = `PROSPER Prob.`, y = PRISM, color=grain.size), se=F) +
  geom_point(data=crcs %>% mutate(ad = abs(mn_prob_d)), aes(x = `PROSPER Prob.`, y = PRISM, fill=grain.size, size=ad), shape=21, alpha=0.6) + 
  #no size scaling:
  # geom_point(data=crcs, aes(x = `PROSPER Prob.`, y = PRISM, fill=grain.size), shape=21, alpha=0.6) + 
  stat_cor(data=crcs, aes(x = `PROSPER Prob.`, y = PRISM, color=grain.size), p.accuracy = 0.0001, size=4, show.legend = FALSE) + 
  scale_fill_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .8)(2))) +
  scale_color_manual('Grain Size:', values=c(viridis_pal(begin = .3, end = .8)(2))) +
  scale_x_continuous(limits = c(0.15,1), breaks = c(0.25,0.50,0.75,1.00)) +
  scale_size_continuous('Percent Change:') +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.line = element_line(color="black"),
        legend.key = element_blank()) 
# guides(fill=guide_legend(override.aes=list(shape=21), nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")),
#        shape=guide_legend(nrow=1, direction = "horizontal", title.theme = element_text(size=10,face="bold")))
p.c2
ggsave(filename = paste0(PATH_RCS,"/RCSprospprob_vs_RCSprism.png"), plot = p.c2, height = 5.5, width = 5)


#Table B: Summary Table ----
#species, point count, RCS scores, % diff in rankings
tabB <- RCS.p %>%
  mutate(across(.cols=where(is.numeric), ~round(.x, digits = 2))) %>%
  relocate(species, N_final, 
            #RCS scores
            RCS_WS_spp, RCS_WS_spc, RCS_WS_prsm, 
            RCS_buff_spp, RCS_buff_spc, RCS_buff_prsm, 
            #Percent Change
            dClim_Prob_WS, dClim_Prob_B, mn_prob_d, 
            dClim_Clas_WS, dClim_Clas_B, mn_Clas_d) %>%
  arrange(species) %>%
  rename(Species = species, 
         Num_Occurrence_Points = N_final,
         RCS_ProspProb_HUC = RCS_WS_spp, 
         RCS_ProspClass_HUC = RCS_WS_spc, 
         RCS_Prism_HUC = RCS_WS_prsm,
         RCS_ProspProb_Buff = RCS_buff_spp, 
         RCS_ProspClass_Buff = RCS_buff_spc,
         RCS_Prism_Buff = RCS_buff_prsm,
         PctD_Prob_HUC = dClim_Prob_WS,
         PctD_Prob_Buff = dClim_Prob_B,
         Mean_PctD_Prob = mn_prob_d,
         PctD_Class_HUC = dClim_Clas_WS, 
         PctD_Class_Buff = dClim_Clas_B, 
         Mean_PctD_Class = mn_Clas_d)
write.csv(tabB, file = paste0(PATH_RCS, "/TableB_SpeciesRCS_Summary.csv"), row.names = F)


#Appendix Table: All results ----

AOO <- read.csv(paste0(PATH_RCS, "/AOO/AOO_2022-08-22.csv")) %>% dplyr::select(scientific_name, huc_area_sqkm, buff_area_sqkm)
CS <- CS.buff.prsm %>% 
  mutate(species = gsub("_", " ", species)) %>%
  filter(species %in% species.list$Species) %>%
  pivot_wider(names_from = value_type, values_from = c(mean, sd)) %>%
  left_join(., CS.buff %>%
              pivot_wider(names_from = value_type, values_from = c(mean, sd))) %>%
  rename_with(.cols = where(is.numeric), ~ paste0(.x, "_buf")) %>%
  left_join(., CS.huc.prsm %>% 
              dplyr::select(-X) %>%
              mutate(species = gsub("_", " ", species)) %>%
              filter(species %in% species.list$Species) %>%
              pivot_wider(names_from = value_type, values_from = c(mean, sd)) %>%
              rename_with(.cols = where(is.numeric), ~ paste0(.x, "_huc"))) %>%
  left_join(., CS.huc %>% 
              mutate(species = gsub("_", " ", species)) %>%
              pivot_wider(names_from = value_type, values_from = c(mean, sd)) %>%
              rename_with(.cols = where(is.numeric), ~ paste0(.x, "_huc"))) %>%
  relocate(species,
           mean_annual_ppt_buf, sd_annual_ppt_buf, mean_annual_ppt_huc, sd_annual_ppt_huc,
           mean_aug_tmax_buf, sd_aug_tmax_buf, mean_aug_tmax_huc, sd_aug_tmax_huc,
           mean_jan_tmin_buf, sd_jan_tmin_buf, mean_jan_tmin_huc, sd_jan_tmin_huc,
           mean_aug_temp_buf, sd_aug_temp_buf, mean_aug_temp_huc, sd_aug_temp_huc,
           mean_baseflow_buf, sd_baseflow_buf, mean_baseflow_huc, sd_baseflow_huc,
           mean_strm_perm_class_buf, sd_strm_perm_class_buf, mean_strm_perm_class_huc, sd_strm_perm_class_huc,
           mean_strm_perm_prob_buf, sd_strm_perm_prob_buf, mean_strm_perm_prob_huc, sd_strm_perm_prob_huc)
apx <- RCS.p %>%
  select(species, N_final, contains("RCS_")) %>%
  relocate(species, N_final, contains("buff"), contains("WS")) %>%
  relocate(species, N_final, contains("spp"), contains("spc"), contains("prsm")) %>%
  left_join(., AOO %>% relocate(contains("buff"), contains("huc")), by=c("species"="scientific_name")) %>%
  left_join(., CS) %>%
  mutate(across(.cols=where(is.numeric), ~ round(.x, digits = 3))) %>%
  arrange(species)
write.csv(apx, file=paste0(PATH_RCS, "/Appendix_Summary.csv"), row.names = F)

#RXV: Comparing spc vs spp; WS vs buff ----
png(paste0(PATH_RCS, "/spc_spp_huc_buff_comparisons.png"))
par(mfrow=c(2,2))
plot(RCS$RCS_WS_spc, RCS$RCS_WS_spp, main = "Watershed SPC vs SPP")
abline(lm(RCS$RCS_WS_spc ~ RCS$RCS_WS_spp), col = "red")
plot(RCS$RCS_buff_spc, RCS$RCS_buff_spp, main = "Snapped SPC vs SPP")
abline(lm(RCS$RCS_buff_spc ~ RCS$RCS_buff_spp), col = "red")
plot(RCS$RCS_WS_spc, RCS$RCS_buff_spc, main="SPC watershed vs snapped")
abline(lm(RCS$RCS_WS_spc ~ RCS$RCS_buff_spc), col = "red")
plot(RCS$RCS_WS_spp, RCS$RCS_buff_spp, main="SPP watershed vs snapped")
abline(lm(RCS$RCS_WS_spp ~ RCS$RCS_buff_spp), col = "red")
dev.off()