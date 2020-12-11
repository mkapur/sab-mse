## Climate-rec paper
## Figure 1

library(sf)
library(tibble)
require(ggplot2)
require(dplyr)
require(patchwork)
require(ggsidekick)
require(rgdal)
library(marmap)
require(rerddap)
require(here)

## political boundaries ----
load(here("input","downloads","eez_nepac_regions.RDA"))
regions <- eez_nepac_regions
usa <- map_data("world") %>% filter(region %in% c("USA","Canada"))
usa$long  <- usa$long+360

## LME Clips ----
# BS
BS <- tibble(X = seq(-200 + 360, -155 + 360, length.out = 100),
             Y = seq(55, 65, length.out = 100)) %>% as.matrix()
n50BS_poly <- st_linestring(BS)
n50BS_geom <- st_sfc(list(n50BS_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n50BS_shape <- st_sf(n50BS_geom)
bering <- st_crop(regions %>% filter(Region_Name == 'Alaska'), 
                  n50BS_shape)

##AI
AI <- tibble(X = seq(-200 + 360, -170 + 360, length.out = 100),
             Y = seq(45, 65, length.out = 100)) %>% as.matrix()
n50AI_poly <- st_linestring(AI)
n50AI_geom <- st_sfc(list(n50AI_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n50AI_shape <- st_sf(n50AI_geom)
aleut <- st_crop(regions %>% filter(Region_Name == 'Alaska'), 
                 n50AI_shape)

## Gulf
GOA <- tibble(X = seq(-170 + 360, -130 + 360, length.out = 100),
              Y = seq(45, 60, length.out = 100)) %>% as.matrix()
n50GOA_poly <- st_linestring(GOA)
n50GOA_geom <- st_sfc(list(n50GOA_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n50GOA_shape <- st_sf(n50GOA_geom)
gulf <- st_crop(regions %>% filter(Region_Name == 'Alaska'), 
                n50GOA_shape)

## S Pt Conception
Conc <- tibble(X = seq(-135 + 360, -115 + 360, length.out = 100),
               Y = seq(28, 36, length.out = 100)) %>% as.matrix()
n50Conc_poly <- st_linestring(Conc)
n50Conc_geom <- st_sfc(list(n50Conc_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n50Conc_shape <- st_sf(n50Conc_geom)
ptcon <- st_crop(regions %>% filter(Region_Name == 'US West Coast'), 
                 n50Conc_shape)

## climate layers ----
#* ocean currents ----
shape <- readOGR(dsn = "C:/Users/MKapur/Dropbox/UW/sab-growth/raw_data/Major_Ocean_Currents_arrowPolys_30m_8", 
                 layer = 'Major_Ocean_Currents_arrowPolys_30m_8')

shapefile_df <- fortify(shape)
shapefile_df$long <- shapefile_df$long+360


#* ocean bathy -- slowish ----
# NEPbath0 <- getNOAA.bathy(lon1 = -180, lon2 = -120,
#                           lat1 = 25, lat2 = 65, 
#                           resolution = 4)
# row.names(NEPbath0) <- as.numeric(row.names(NEPbath0))+360
# 
# NEPbath <- fortify(NEPbath0)
# NEPbath$x <- NEPbath$x+360


#* ocean SST ----
# https://rmendels.github.io/pices2017.nb.html
sstInfo <- info('jplMURSST41')
# get latest daily sst
murSST <- griddap(sstInfo, latitude = c(22., 75), 
                  longitude = c(-179.99, -115), time = c('last','last'), 
                  fields = 'analysed_sst', stride = 2)
murSST$data$lon= murSST$data$lon+360

##AI Chunk - SHOULD RUN AT 179 ONCE CORRECTED
murSSTB <- griddap(sstInfo, latitude = c(22., 75), 
                  longitude = c(170,180), time = c('last','last'), 
                  fields = 'analysed_sst', stride = 2)
# murSSTB$data$lon= (murSSTB$data$lon-360)*-1
summary(murSSTB$data$lon)
# save(list(c(murSST,murSSTB)), file = here("input","downloads","murSST.RDA"))
# load(here("input","downloads","murSST.RDA"))
mycolor <- colors$temperature
## use this to check the coords are correct
# ggplot(data = murSSTB$data, aes(x = lon, y = lat, fill = analysed_sst)) +
#   # geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +
#   geom_raster(interpolate = FALSE) +
#   scale_fill_gradientn(colours = mycolor, na.value = NA) +
#   theme_bw() + ylab("latitude") + xlab("longitude") +
#   # coord_fixed(1.3, xlim = c(-140+360, -105+360),  ylim = c(22., 51.)) +
#   ggtitle("Latest MUR SST")

## settings
land.col <- "cornsilk1"
EEZ.border.col <- 'black'
EEZ.fill.col <- c('grey22','grey44','grey66')
currents.col <- c("#000000", "#009E73", "#0072B2","#e79f00","#e79f00")
# c("gold", "dodgerblue","#ABDDA4", "#ABDDA4" ,"grey22")

## manual legend placement
x1 <- rep(180, 4); x1end <- rep(190, 4)
y1 <- rev(seq(30,35,length.out = 4))

## Fig 1: SST  ----
ocsst <- ggplot(data = regions) + 
  theme_sleek() +
  labs(x ="",y="", fill = 'SST (mean *C)') +
  geom_raster(data = murSST$data, aes(x = lon, y = lat, 
                                      fill = analysed_sst),
              interpolate = TRUE) +
  
  geom_raster(data = murSSTB$data, aes(x = lon, y = lat, 
                                      fill = analysed_sst),
              interpolate = TRUE) +
  scale_fill_gradientn(colours = colors$temperature, na.value = NA) +
  ## EEZ
  geom_sf(lwd = 1, col = 'black', fill = EEZ.fill.col, alpha = 0.2) +
  ## land
  # geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = NA) +
  coord_sf(xlim = c(175, 240), ylim = c(25.5, 60)) 


# Fig 2: LMEs ----
lmes <- ggplot(data = regions) + 
  theme_sleek() +
  ## EEZ and land
  geom_sf(lwd = 1, col = EEZ.border.col, fill = EEZ.fill.col, alpha = 0.2) +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = land.col) +

 
  ## manual legend
  annotate("segment", x = x1, xend = x1end, y = y1, yend = y1,
           colour = currents.col[1:4], size = 1.1) +
  annotate("text", x = x1+20, y = y1,
           label =  c('Alaskan Current',
                      'N. Pacific Current', 
                      'S. California Bight',
                      'California Current'),
           size = 2) +
  ##* add LMEs ----
  geom_sf(data = bering, fill = 'blue', alpha = 0.9, color = NA) +
  geom_sf(data = aleut, fill = 'blue', alpha = 0.9, color = NA) +
  geom_sf(data = gulf, fill = 'pink', alpha = 0.9, color = NA) +
  geom_sf(data = ptcon, fill = 'grey22', alpha = 0.9, color = NA) +
  labs(x ="",y="") +
  ## add currents
  geom_polygon(data = shapefile_df,
               aes(x = long, y = lat, group = group, fill = id),
               fill = rep(currents.col,
                          length(shapefile_df$id)/5), size = 0.2) +
  geom_label(aes(x = x1[1]+5, y = 55, label = "BS/AI"), fill = "grey88",
             size = 2) +
  geom_label(aes(x = x1[1]+30, y = 55, label = "GOA"), fill = "grey88",
             size = 2) +
  geom_label(aes(x = x1[1]+50, y = 36, label = "S. Pt Conception"), fill = "grey88",
             size = 2) + 
  coord_sf(xlim = c(175, 240), ylim = c(25.5, 65)) 


## ocsst takes 5 mins to render, 
ggsave(lmes   | ocsst,
       file = here("figs","lmes_sst.png"),
       width = 10, height = 8, dpi = 420)

## option 2: with bathy ----
autoplot(NEPbath0, geom=c("r", "c"), colour="white", size=0.05) +
  theme_sleek() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = land.col) +
  # geom_polygon(data = fortify(regions), aes(x = long, y = lat), fill = land.col) +
  
  coord_sf(xlim = c(180, 240), ylim = c(25.5, 65)) +
  labs(x ="",y="") +theme(legend.position = 'none')

ggsave(last_plot(),
       file = here("figs","climrec_fig1_bathy.png"),
       width = 6, height = 6, dpi = 420)
