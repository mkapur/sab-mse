
# Hi Maia, This should allow you to plot the regions
# - The sf package has good online support and interfaces with ggplot2
# - Extra points / lines usually need to be converted to sf objects first
# - Plot extra sf lines / points with more geom_sf() calls (e.g. below)



# Plot Regions
# Load data
load("eez_nepac_regions.rda")
# Plot equal area
map_regions <- ggplot(data = eez_nepac_regions) +
	geom_sf() +
	# geom_sf(data = , aes(x = , y = )) + # Point data should be an sf object
	# This one gives an "equal area" projection:
	# coord_sf(crs = "+proj=laea +lat_0=48 +lon_0=210 +ellps=GRS80 +units=m +no_defs") +
	# This gives a mercator projection:
	coord_sf() +
	theme_classic()
# Save
ggsave(here::here("figs", "map-regions.pdf"), width = 6, height = 4)
ggsave(here::here("figs", "map-regions.png"), width = 6, height = 4)
