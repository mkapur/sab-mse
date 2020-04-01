require(ggplot2)
require(dplyr)
require(patchwork)

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3",
                "#0072B2", "#D55E00", "#CC79A7", "navy", "#F0E442" )
cbbPalette <- c("grey22", "seagreen2", "goldenrod", "skyblue",
                "blue", "brown", "brown", "pink" )
## Figure 1 map of strata ----
usa <- map_data("world") 
load("C:/Users/mkapur/Dropbox/UW/sab-idx/runs/2020-01-23_nx=500_Triennial_WCGBTS_BCs_BCo_AK_DOM_LL_GOA_baseQ=AK_DOM_LL1980_2018/Data_Geostat.Rdata")

survLims <-  Data_Geostat %>% 
  group_by(Region) %>% 
  summarise(ymin = min(Lat),
            ymax = max(Lat), 
            xmin = min(Lon),
            xmax = max(Lon))

Data_Geostat %>% filter(Region == 'BC') %>%
  group_by(Survey) %>% 
  summarise(lat_min = min(Lat),
            lat_max = max(Lat), 
            lon_min = min(Lon),
            lon_max = max(Lon))

Data_Geostat %>% filter(Region == 'BC') %>%
  # group_by(Survey) %>% 
  filter(Lon > -126) 
  summarise(round(max(Lon),10))

## load polygons  
load(paste0("./figures/spdf_fortified_BC.Rdata"))
load(paste0("./figures/spdf_fortified_US.Rdata"))

## clockwise from A1; two for A2 
regLims <- data.frame(ymax = c(65,65,65,55,65,50,47,36),
                      ymin = c(50,50,55,50,50,47,36,30), 
                      xmax = c(-145,-138, -130, -130, rep(-120,3),-115), 
                      xmin = c(-180,-145, -138, -138, rep(-130,4)) )

mgmtLims <- data.frame(ymax = c(65, 49),
                      ymin = c(49, 30), 
                      xmax = c(-180, -115), 
                      xmin = c(-132.2, -132.2))

demoLims <- data.frame(ymax = c(65,65,65,50,36),
                       ymin = c(50,50,50,36,30), 
                       xmin = c(-180,-145, -130,-130,-130), 
                       xmax = c(-145,-130, -115, -115, -115))

ggplot() + geom_polygon(data = usa, aes(x = long, y = lat, group = group), 
                        fill = 'grey22') +
  kaputils::theme_mk(base_size = 16) + 
  theme(axis.title =element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
  scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,60,10), 
                     labels =  paste(seq(30,60,10), "°N"))  +
  
  ## CURRENT MGMT BOUNDARIES
  # geom_polygon(data = spdf_fortified_US, aes( x = long, y = lat, group = group),
  #              fill= 'red', color="red") +
  # geom_polygon(data = spdf_fortified_BC, aes( x = long, y = lat, group = group),
  #              fill= 'red', color="red")+
  geom_polygon(data = spdf_fortified_world, aes( x = long, y = lat, group = group),
               fill= NA, color="red") +
  # geom_rect(data = mgmtLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
            # fill =NA, size = 1, colour = 'red') +
  
  # Complexity: actual survey boundaries
  geom_rect(data = survLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  fill = NA , size = 1, colour = 'blue',linetype = 'dotted', alpha = 0.2) +

  
  # DEMOGRAPHIC BOUNDARIES
  geom_rect(data = demoLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  fill = NA , size = 1, colour = 'black',linetype = 'dashed', alpha = 0.2) +

  # ## OM SUB AREAS
  geom_rect(data = regLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
            fill = rev(cbbPalette[1:8]), size = 1, colour = NA, alpha = 0.3) +
  geom_label(aes(
      x = c(rep(-125, 4), -132, -140, -155),
      y = c(33, 40, 49.5, rep(53, 4)),
      label = c("W1", "W2", "B1", "B2", "A3", "A2", "A1")
    ),
    size = 5,
    fill = cbbPalette[c(1:6,8)],
    color = c("grey88", rep('black', 3), 'grey88', rep('black',2))
  ) +
  labs(x = "", y = "") +
  coord_quickmap()  



ggsave(plot = last_plot(),
       file = paste0("./figures/Fig1_strata_mapsC_WithSurvey.png"),
       width = 10, height = 8, units = 'in', dpi = 420)


## Figure 2 panel of OM Indices (made by me with VAST) ----

## load tables_for_ss3.csv
OM1File <- "C:/Users/mkapur/Dropbox/UW/sab-idx/runs/2020-01-23_nx=500_Triennial_WCGBTS_BCs_BCo_AK_DOM_LL_GOA_baseQ=AK_DOM_LL1980_2018/Table_for_SS3.csv"
vastc <- read.csv(OM1File) %>%
  mutate(TYPE = 'Abundance', Source = 'VAST',
         lci = Estimate_metric_tons-SD_mt,
         uci = Estimate_metric_tons+SD_mt) %>%
  select(Year, Fleet, Estimate_metric_tons, SD_log, TYPE, Source, uci, lci ) 
for(i in 1:nrow(vastc)){
  vastc$Fleet2[i] <- ifelse(vastc$Fleet[i] == "California_current",
                            "WC", 
                            ifelse(vastc$Fleet[i] == "British_Columbia", "BC",
                                   "AK"))
  if(vastc$Fleet[i] == 'AllAreas')   vastc$Fleet2[i] <- "ALL"
}

p1 <- vastc %>%
  filter(Fleet != 'Eastern_Bering_Sea' & Fleet != 'AllAreas') %>%
  ggplot(.,   aes(x = Year, y = Estimate_metric_tons, col = Fleet)) +
  kaputils::theme_mk(base_size = 16) + theme(legend.position = c(0.75,0.75))+
  scale_color_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  labs(x = 'Year', y = 'Estimate (mt)', title = paste0('OM1 Indices by Region'), 
       fill = "", color = "") +
  geom_line(lwd = 0.9)+
  # geom_point(pch = 1, cex = 3) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci, fill = Fleet),
              alpha = 0.2,
              show.legend = FALSE)


p2 <- vastc %>%
  filter(Fleet == 'AllAreas') %>%
  ggplot(.,   aes(x = Year, y = Estimate_metric_tons, col = Fleet)) +
  kaputils::theme_mk(base_size = 16) + theme(legend.position = c(0.75,0.75)) +
scale_color_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  labs(x = 'Year', y = 'Estimate (mt)', title = paste0('OM2 PLACEHOLDER'), 
       fill = "", color = "") +
  geom_line(lwd = 0.9)+
  # geom_point(pch = 1, cex = 3) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci, fill = Fleet),
              alpha = 0.2,
              show.legend = FALSE)

## p3 is allareas from same run as p1
p3 <- vastc %>%
  filter(Fleet == 'AllAreas') %>%
  ggplot(.,   aes(x = Year, y = Estimate_metric_tons, col = Fleet)) +
  kaputils::theme_mk(base_size = 16) + theme(legend.position = c(0.75,0.75))+
  scale_color_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  labs(x = 'Year', y = 'Estimate (mt)', title = paste0('OM3 Panmictic Index'), 
fill = "", color = "") +
  geom_line(lwd = 0.9)+
  # geom_point(pch = 1, cex = 3) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci, fill = Fleet),
              alpha = 0.2,
              show.legend = FALSE)




ggsave(plot = (p1  | p2  | p3),
       file = paste0("./figures/Fig2_OM_indices.png"),
       width = 17, height = 10, units = 'in', dpi = 720)



## Figure X input growth curves by OM ----

## OM2 uses values from kapur et al. 2019



# https://github.com/mkapur/sab-growth/blob/master/SAB_plot_master.R L112

# ypreds <- read.csv("C:/Users/mkapur/Downloads/SAB_predicts_2019-10-04_phase2.csv")
# ypreds$gamREG <- ypreds$gamREG
# levels(ypreds$REG) <- c('Alaska','British Columbia','US West Coast')
# levels(ypreds$Sex) <- c('Females','Males')
# levels(ypreds$Period) <- c('pre-2010','2010-2018','All Years')
# for(i in 1:nrow(ypreds)){
#   ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1,
#                               'All Years', paste(ypreds$Period[i]))
# }
# fd_summary_gamREG <- ypreds %>%
#   # filter(Age < 75) %>%
#   group_by(Age, Sex, gamREG,Period) %>%
#   dplyr::summarise(meanL = mean(Length_cm), sdmeanL = sd(Length_cm), meanPred = mean(Predicted))
# saveRDS(fd_summary_gamREG, file = "./OM2_RegionalCurves.Rdata")


p2 <- ggplot(readRDS("./OM2_RegionalCurves.Rdata"), aes(x = Age, col = gamREG, group = gamREG)) +
  kaputils::theme_mk(base_size = 16) +
  scale_color_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +  
  scale_alpha(guide = 'none') +
  scale_y_continuous(limits = c(0,110)) +
  scale_x_continuous(limits = c(0,65)) +
  geom_line(aes(y = meanPred, col = gamREG), lwd = 1.1)+
  labs(y = 'Length (cm)', x= 'Age (years)', col = "") +
  facet_wrap(~Sex +Period, ncol = 4)

## use ONEREG from dropbox
# ypreds <- read.csv("C:/Users/mkapur/Downloads/ONEREG_SAB_predicts_2019-10-29.csv")
# ypreds$gamREG <- ypreds$gamREG
# levels(ypreds$REG) <- c('Alaska','British Columbia','US West Coast')
# levels(ypreds$Sex) <- c('Females','Males')
# levels(ypreds$Period) <- c('pre-2010','2010-2018','All Years')
# for(i in 1:nrow(ypreds)){
#   ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1,
#                               'All Years', paste(ypreds$Period[i]))
# }
# fd_summary_gamREG <- ypreds %>%
#   # filter(Age < 75) %>%
#   group_by(Age, Sex, gamREG,Period) %>%
#   dplyr::summarise(meanL = mean(Length_cm), sdmeanL = sd(Length_cm), meanPred = mean(Predicted))
# 
# saveRDS(fd_summary_gamREG, file = "./OM3_RangewideCurve.Rdata")


p3 <- ggplot(readRDS("./OM3_RangewideCurve.Rdata"), aes(x = Age)) +
  kaputils::theme_mk(base_size = 16) +theme(legend.position = 'none') +
  scale_color_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +  
  scale_alpha(guide = 'none') +
  scale_y_continuous(limits = c(0,110)) +
  scale_x_continuous(limits = c(0,65)) +
  geom_line(aes(y = meanPred, col = gamREG), lwd = 1.1)+
  labs(y = 'Length (cm)', x= 'Age (years)', col = "") +
  facet_wrap(~Sex , ncol = 4)
  
  
  ggsave(plot = (p2  / p3),
         file = paste0("./figures/FigX_GrowthCurves.png"),
         width = 17, height = 10, units = 'in', dpi = 720)
  