
# Hi Maia, here's how I added the breakpoint lines, by building each one up in steps
# - First as a tibble
# - Then as a st_linestring
# - Then as a st_sfc
# - Then as a sf object (via st_sf)
# Let me know if you find a simpler way that works! I was in a hurry, so I just
# used this.

# Make breakpoints shapefiles
# 36 N
n36 <- tibble(X = seq(-120 + 360, -129 + 360, length.out = 100),
							Y = rep(36, 100)) %>% as.matrix()
n36_poly <- st_linestring(n36)
n36_geom <- st_sfc(list(n36_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n36_shape <- st_sf(n36_geom)
save(n36_shape, file = here("input","downloads","n36_shape.rda"))


# 50 N
n50 <- tibble(X = seq(-123 + 360, -138 + 360, length.out = 100),
							Y = rep(50, 100)) %>% as.matrix()
n50_poly <- st_linestring(n50)
n50_geom <- st_sfc(list(n50_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n50_shape <- st_sf(n50_geom)
save(n50_shape, file = here("input","downloads","n50_shape.rda"))

# 50 BC
n50BC <- tibble(X = seq(-125 + 360, -138 + 360, length.out = 100),
              Y = rep(50, 100)) %>% as.matrix()
n50BC_poly <- st_linestring(n50BC)
n50BC_geom <- st_sfc(list(n50BC_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
n50BC_shape <- st_sf(n50BC_geom)

# 145 W
w145 <- tibble(X = rep(-145 + 360, 100),
							 Y = seq(55, 62, length.out = 100)) %>% as.matrix()
w145_poly <- st_linestring(w145)
w145_geom <- st_sfc(list(w145_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
w145_shape <- st_sf(w145_geom)

# 130 W
w130 <- tibble(X = rep(-130 + 360, 100),
							 Y = seq(45, 55, length.out = 100)) %>% as.matrix()
w130_poly <- st_linestring(w130)
w130_geom <- st_sfc(list(w130_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
w130_shape <- st_sf(w130_geom)

# BC test angle
BCreg <- regions %>% filter(Region_Name == 'British Columbia')

splitLat <- function(sample_x = 25, xmin, xmax, ymin, ymax){
  ydiff = ymax - ymin
  xdiff = xmax - xmin
  int = ymin - ydiff/xdiff * xmin
  # split_y <- (ydiff/xdiff)*sample_x + int
  return(list(ydiff, xdiff, int))
}
xmin = -140 + 360
xmax = -125 + 360
ymin = 45
ymax = 55
xdomain <- seq(xmin,xmax, length.out = 100)

outerLower <- matrix(c(xmin, ymin,xmax,ymin,xmax,ymax,xmin, ymin),
                ncol=2, byrow=TRUE)
pts = list(outerLower)
pl1 = st_polygon(pts)
BCangle_geom <- st_sfc(list(pl1), crs = "+proj=longlat +datum=WGS84 +no_defs")
BCangle_shape <- st_sf(BCangle_geom)
B2 <- st_intersection(BCreg,BCangle_shape)

outerUpper <- matrix(c(xmin, ymin,xmin, ymax,xmax,ymax,xmin, ymin),
                     ncol=2, byrow=TRUE)
pts = list(outerUpper)
p12 = st_polygon(pts)
BCangle_geom <- st_sfc(list(p12), crs = "+proj=longlat +datum=WGS84 +no_defs")
BCangle_shape <- st_sf(BCangle_geom)
B1 <- st_intersection(BCreg,BCangle_shape)
save(B1,file =  "./input/cleaned/sub_area_clips_B1.Rdata")
save(B2,file =  "./input/cleaned/sub_area_clips_B2.Rdata")

par(mfrow = c(1,2))
plot(BCreg$geometry)
plot(pl1, col = 'gold')
plot(p12, add = TRUE, col = 'blue')

BCbath <- st_read("./_writing/figures/BC_EEZ_bathy.gdb")
# install_github("ericpante/marmap") ## to use with GGplot
library(marmap)
xmin = -140
xmax = -125 
ydif = splitLat(sample_x = 25, xmin, xmax, ymin, ymax)[[1]]
xdif = splitLat(sample_x = 25, xmin, xmax, ymin, ymax)[[2]]
int = splitLat(sample_x = 25, xmin, xmax, ymin, ymax)[[3]]
blues <- colorRampPalette(c("red","purple","blue",
                            "cadetblue1","white"))
BCbath <- getNOAA.bathy(lon1 = -140, lon2 = -120,
                        lat1 = 45, lat2 = 65, resolution = 4)

# https://stackoverflow.com/questions/27214282/add-bathymetry-lines-to-ggplot-using-marmap-package-and-getnoaa-bathy

# BCangle <- tibble(X = xdomain,
#                Y = splitLat(sample_x = xdomain, xmin, xmax, 
#                             ymin = 50, ymax = 54)) %>% as.matrix()
# BCangle_poly <- st_linestring(BCangle)
# BCangle_geom <- st_sfc(list(BCangle_poly), crs = "+proj=longlat +datum=WGS84 +no_defs")
# BCangle_shape <- st_sf(BCangle_geom)
# #first, get the boundaries of the shapefile
# bbox <- BCangle_shape %>% 
#   st_bbox() %>% st_as_sfc( crs = "+proj=longlat +datum=WGS84" )
# 
# #      xmin      ymin      xmax      ymax 
# # 122.29929  11.71056 124.42607  14.50061 
# 
# spdf <- regions %>% filter(Region_Name == 'British Columbia') 
# 
# BCreg <- spdf %>% #check if a line the spatial df intersecgt with the defined boundary-box
#   mutate( passes_through_box = as.numeric( st_intersects(spdf, bbox) ) ) %>%
#   group_by( SN ) %>%
#   mutate( passes_through_box_anyime = ifelse( any( passes_through_box == 1),
#                                               "yes", "no" ) )
# 
# 
# # https://stackoverflow.com/questions/32233576/how-to-clip-a-polygon-shapefile-by-another-polygon-shapefile-in-r
# gClip <- function(shp, bb){
#   if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
#   else b_poly <- as(extent(bb), "SpatialPolygons")
#   gIntersection(shp, b_poly, byid = T)
# }
# 
# 
# zones_clipped_w <- gClip(BCreg, BCangle_shape)
## clip EEZ to lims
# https://gis.stackexchange.com/questions/282524/cropping-sf-object-in-r
## xlims are negvals plus 360

## these are the defaults
# regLims <- data.frame(ymax = c(75,75,55,55,50,49,36),
#                       ymin = c(45,50,50,50,45,36,30),
#                       xmax = c(-145,-132, -130, -122,-124,-121,-115),
#                       xmin = c(-195,-145, -145, rep(-130,4)) )

## dropping 130
regLims <- data.frame(ymax = c(75,75,55,50,50,36),
                      ymin = c(45,50,50,45,36,30),
                      xmax = c(-145,-132, -122,-124,-121,-115),
                      xmin = c(-195,-145, -145, rep(-145,3)) )

reglims$Region <- 7:1
reglims$Region <- ifelse(regLims$Region == 7, 4,
                         ifelse(regLims$Region < 7 &regLims$Region > 4 , 3))


regLims[-2,]

regLims$xmax <- regLims$xmax + 360
regLims$xmin <- regLims$xmin + 360

clips <- list()
for(i in 1:nrow(regLims)){
  if(i < 3){
  clips[[i]] <- st_crop(regions %>% filter(Region_Name == 'Alaska'), 
                        xmin = regLims$xmin[i], 
                        xmax = regLims$xmax[i], 
                        ymin = regLims$ymin[i], 
                        ymax = regLims$ymax[i])
  } else if (i > 2 & i < 5){
    clips[[i]] <- st_crop(regions %>% filter(Region_Name == 'British Columbia'), 
                          xmin = regLims$xmin[i], 
                          xmax = regLims$xmax[i], 
                          ymin = regLims$ymin[i], 
                          ymax = regLims$ymax[i])
  } else{
    clips[[i]] <- st_crop(regions %>% filter(Region_Name == 'US West Coast'), 
                          xmin = regLims$xmin[i], 
                          xmax = regLims$xmax[i], 
                          ymin = regLims$ymin[i], 
                          ymax = regLims$ymax[i])
  }
}


save(clips,file =  "./input/cleaned/sub_area_clips_50N.Rdata")
BCtest <- st_crop(regions %>% filter(Region_Name == 'British Columbia'), 
                  xmin = 215, xmax = 228, ymin = 53.5, ymax = 75)

# install.packages("PBSmapping")
# require(PBSmapping)
# # BCclip <- regions[w130_shape,]
# df <- importShapefile("./_writing/figures/World_EEZ_v11_20191118_HR_0_360/eez_v11_0_360.shp")
# df_sub <- clipLines(df, xlim = c(50 , 55) , ylim = c(130 , 115), keepExtra = TRUE )
# dfSL <- PolySet2SpatialLines( df_sub )
# # Plot
# ggplot(data = dfSL) +
# 	geom_sf() 
# 	geom_sf(data = n36_shape, col = "red") +
# 	geom_sf(data = n50_shape, col = "red") +
# 	geom_sf(data = w145_shape, col = "red") +
# 	geom_sf(data = w130_shape, col = "red") +
# 	# coord_sf(crs = "+proj=laea +lat_0=48 +lon_0=210 +ellps=GRS80 +units=m +no_defs") +
# 	coord_sf(xlim = c(165, 245), ylim = c(30, 65)) +
# 	theme_classic()
# Save
# ggsave(here::here("figs", "map-sub-areas-maia.pdf"), width = 6, height = 4)
# ggsave(here::here("figs", "map-sub-areas-maia.png"), width = 6, height = 4)
