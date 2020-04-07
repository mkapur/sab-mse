require(ggplot2)
require(dplyr)
require(patchwork)

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3",
                "#0072B2", "#D55E00", "#CC79A7", "navy", "#F0E442" )
# cbbPalette <- c("grey22", "seagreen2", "goldenrod", "skyblue",
#                 "blue", "brown", "pink" )
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
regLims <- data.frame(ymax = c(65,65,65,65,50,49,36),
                      ymin = c(50,50,50,50,49,36,30), 
                      xmax = c(-145,-132, -130, rep(-120,3),-115), 
                      xmin = c(-180,-145, -132, rep(-130,4)) )

mgmtLims <- data.frame(ymax = c(65, 65, 49),
                      ymin = c(49, 49, 30), 
                      xmax = c(-180, -115, -115), 
                      xmin = c(-132, -132, -132))

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
  # geom_polygon(data = spdf_fortified_world, aes( x = long, y = lat, group = group),
  #              fill= NA, color="red") +
  geom_rect(data = mgmtLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  fill =NA, size = 1, colour = 'red') +
  
  # Complexity: actual survey boundaries
  # geom_rect(data = survLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  # fill = NA , size = 1, colour = 'blue',linetype = 'dotted', alpha = 0.2) +

  
  # DEMOGRAPHIC BOUNDARIES
  geom_rect(data = demoLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  fill = NA , size = 1, colour = 'black',linetype = 'dashed', alpha = 0.2) +

  # ## OM SUB AREAS
  geom_rect(data = regLims, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
            fill = rev(cbbPalette[1:7]), size = 1, colour = NA, alpha = 0.3) +
  geom_label(aes(
      x = c(rep(-125, 4), -131.5, -140, -155),
      y = c(33, 40, 49.5, rep(53, 4)),
      label = c("C1", "C2", "B1", "B2", "A3", "A2", "A1")
    ),
    size = 5,
    fill = cbbPalette[c(1:7)],
    color = c("grey88", rep('black', 3), 'grey88', rep('black',2))
  ) +

  labs(x = "", y = "") +
  coord_quickmap()  

ggsave(plot = last_plot(),
       file = paste0("./figures/Fig1_strata_mapsC.png"),
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
vastc <- vastc %>%
  filter(Fleet != 'Eastern_Bering_Sea' & Fleet != 'AllAreas') 

vastc$SubArea[vastc$Fleet == "British_Columbia"] <- "B1, B2, A3 (BC)"
vastc$SubArea[vastc$Fleet == "California_current"] <- "C1, C2 (CC)"
vastc$SubArea[vastc$Fleet == "Gulf_of_Alaska"] <- "A1, A2 (AK)"
vastc$SubArea[vastc$Fleet == "Aleutian_Islands"] <- "A1 (AK)"

genPal <- c('brown','dodgerblue','goldenrod','grey22')
  ggplot(vastc,   aes(x = Year, y = Estimate_metric_tons, col = SubArea)) +
  kaputils::theme_mk(base_size = 16) +
    theme(legend.position = c(0.75,0.75))+
  scale_color_manual(values = genPal) +
  scale_fill_manual(values = genPal) +
  labs(x = 'Year', y = 'Estimate (mt)', 
       fill = "Sub Area(s) (Mgmt. Area)", color = "Sub Area(s) (Mgmt. Area)") +
  geom_line(lwd = 0.9)+
  # geom_point(pch = 1, cex = 3) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci, fill = SubArea),
              alpha = 0.2,
              show.legend = TRUE)
  ggsave(plot = last_plot(),
         file = paste0("./figures/Fig2_OM_indices.png"),
         width = 10, height = 7, units = 'in', dpi = 720)

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

## Figure X pseudo datplot
FleetNames <- c("WCGBTS", "Triennial", "StdTrap","StRs","AK_DOM_LL", "GOA_Trawl",
                "WC_HKL_POT","WC_TWL","BC_Trap","BC_Hook","BC_Trawl", "AK_LL_POT","AK_Trawl")

Fleets0 <- data.frame(FleetNames) %>% mutate(MR = NA)
Fleets0$MR[c(1:2,7:8)] <- "CC"
Fleets0$MR[c(3:4,9:11)] <- "BC"
Fleets0$MR[c(5:6,12:13)] <- "AK"
Fleets0$Type <- rep(NA, 13)
Fleets0$Type[1:6] <- "Survey"
Fleets0$Type[7:13] <- "Fishery"

Fleets <- rbind(Fleets0,Fleets0)
Fleets$Type[14:nrow(Fleets)] <- NA
## Add Age Comps
## CC: All fleets and have ages
Fleets$Type[Fleets$MR=='CC' & is.na(Fleets$Type)] <- "AgeComp"
## BC: Ages mostly from trap and both surveys
Fleets$Type[Fleets$MR=='BC' & is.na(Fleets$Type)] <- "AgeComp"
## AK: Ages from AKPOT and LL Surv
Fleets$Type[Fleets$MR =='AK' & is.na(Fleets$Type) & FleetNames %in% FleetNames[c(5,12)]] <- "AgeComp"
## Add LenComps
Fleets <- rbind(Fleets,Fleets0)
Fleets$Type[27:nrow(Fleets)] <- NA

## CC: from WCGBTS
Fleets$Type[Fleets$MR=='CC' & FleetNames == FleetNames[1] & is.na(Fleets$Type)] <- "LenComp"
## BC: commercial trawl
Fleets$Type[Fleets$MR=='BC' &  FleetNames == FleetNames[11] & is.na(Fleets$Type)] <- "LenComp"
## AK: fixed gear  and trawl fishery
Fleets$Type[Fleets$MR =='AK' & is.na(Fleets$Type) & FleetNames %in% FleetNames[c(12,13)]] <- "LenComp"

Fleets <- Fleets[!is.na(Fleets$Type),]

## basic plot (no years)
ggplot(Fleets, aes(x = FleetNames, y = Type, fill = Type)) +
  kaputils::theme_mk(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90), legend.position = 'none') +
  geom_tile(color = 'white') +
  scale_fill_manual(values = PNWColors::pnw_palette('Starfish')) +
  facet_wrap(.~MR, drop = TRUE, 
             scales = 'free_x') +
  labs(x = "Fleet Names", y = "Data Type")

ggsave(plot = last_plot(),
       file = paste0("./figures/Data_basic.png"),
       width = 10, height = 7, units = 'in', dpi = 720)
## expand grid into YEAR and drop rows if not applicable


dat <- expand.grid(Fleets,
            "Year" = 1970:2019 ) %>%
  mutate(FleetNames)
  mutate( LengthComp = NA,
          AgeComp = NA,
          CAAL = NA,
          Catch1  = NA,
          Catch2  = NA,
          Abund1 = NA,
          Abund2 = NA)


## Fill in survey years
dat$Abund1[dat$MR == 'CC' & dat$Year  %in% 1980:2004] <- FleetNames[2]
dat$Abund2[dat$MR == 'CC' & dat$Year >= 2003] <- FleetNames[1]

dat$Abund1[dat$MR == 'BC' & dat$Year %in% 1991:2009] <- FleetNames[3]
dat$Abund2[dat$MR == 'BC' & dat$Year  %in% 2003:2019] <- FleetNames[4]

dat$Abund1[dat$MR == 'AK' & dat$Year >= 1979 ] <- FleetNames[5]
dat$Abund2[dat$MR == 'AK' & dat$Year  >= 1980] <- FleetNames[6]


## Fill in Comps
dat$AgeComp[dat$MR == 'CC' & dat$Year  %in% 1980:2004] <- FleetNames[2]

dat$AgeComp[dat$MR == 'BC' & dat$Year %in% 1991:2009] <- FleetNames[3]

dat$AgeComp[dat$MR == 'AK' & dat$Year >= 1979 ] <- FleetNames[5]

 data.frame("Year" = 1960:2019,
                  LengthComp =NA,
                  AgeComp = NA,
                  CAAL = NA,
                  Catches = NA,
                  Abund = NA)



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
  