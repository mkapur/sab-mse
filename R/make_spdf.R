# https://coastalmap.marine.usgs.gov/GISdata/basemaps/boundaries/eez/NOAA/useez_noaa.htm
# and sent from Lisa Lacko 31 Mar 2020

# http://www.marineregions.org/downloads.php
# v3 2020 03 17 marine and lone zones union eez

# worldEZ
worldReg0 <- rgdal::readOGR( "./figures/World_EEZ_v11_20191118_HR_0_360",
                          layer = "eez_boundaries_v11_0_360")

worldReg1 <- spTransform(worldReg0, CRS("+proj=longlat")) # reproject

spdf_fortified_world0 <- broom::tidy(worldReg1, region = "TERRITORY1") 

spdf_fortified_world1 <- spdf_fortified_world0 %>%
  filter(long < 245 & long > 160  ) %>%
  filter( lat < 75 & lat > 30) %>%
mutate(long = ifelse(long > 180, -1*(360-long), long)) 


ids <- spdf_fortified_world1 %>%
  group_by(id) %>%
  summarise(mal = max(lat), mil = min(long)) %>%
  filter(mal < 35 & mil < -125)  %>%
  select(id) %>%
  data.frame()

spdf_fortified_world <- spdf_fortified_world1 %>%
  filter(!(id %in% c(ids$id)) )

save(spdf_fortified_world,file = paste0("./figures/spdf_fortified_world.Rdata"))

bcReg0 <- rgdal::readOGR( "./figures/BCEEZ",
                          layer = "canadaEEZ")
bcReg1 <- spTransform(bcReg0, CRS("+proj=longlat")) # reproject

spdf_fortified_BC <- broom::tidy(bcReg1, region = "NAME") %>%
  filter(long < -124 & lat < 56 )



save(spdf_fortified_BC,file = paste0("./figures/spdf_fortified_BC.Rdata"))

usReg0 <- rgdal::readOGR( "./figures/useez",
                          layer = "useez")
usReg1 <- spTransform(usReg0, CRS("+proj=longlat")) # reproject

spdf_fortified_US0 <- broom::tidy(usReg1, region = "NAME") %>%
  filter(long < -115 & lat > 30 ) 
## drop hawaii
ids <- spdf_fortified_US0 %>%
  group_by(id) %>%
  summarise(mal = max(lat), mil = min(long)) %>% 
  filter(mal < 35 & mil < -125) %>%
  select(id) %>% 
  data.frame()

spdf_fortified_US <- spdf_fortified_US0 %>%
  filter(!(id %in% c(ids$id)) & order > 5 )
save(spdf_fortified_US,file = paste0("./figures/spdf_fortified_US.Rdata"))

## sanity chex
ggplot() +
  # geom_polygon(data = spdf_fortified_US, aes( x = long, y = lat, group = group),
  #              fill=NA, color="red") +
  # geom_polygon(data = spdf_fortified_BC, aes( x = long, y = lat, group = group),
  #              fill="red", color="black") +
  geom_polygon(data = spdf_fortified_world, aes( x = long, y = lat, group = group),
               fill="blue", color="blue") +

  theme_void() 
