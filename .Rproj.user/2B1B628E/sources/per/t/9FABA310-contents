require(ggplot2)
require(dplyr)
require(patchwork)

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "navy", "#F0E442" )

## Figure 1 map of strata ----
usa <- map_data("world") 

p1 <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = 'grey22') +
  coord_quickmap() + 
  kaputils::theme_mk(base_size = 16) +  theme(axis.title =element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
  scale_y_continuous(expand = c(0,0), limits = c(30,71), breaks = seq(30,70,10), 
                     labels =  paste(seq(30,75,10), "°N"))  +
  geom_hline(yintercept = 49, lwd = 1.1, col = 'seagreen') +
  geom_vline(xintercept = -132, lwd = 1.1, col = 'seagreen') +
  geom_rect(fill = 'white', aes(xmin = -180, xmax = -132.5,
                                ymin = 30, ymax = 48.5)) +
  geom_label(aes( x = c(rep(-122,2),-150), 
                  y = c(40,rep(53,2)),
                 label = c("California Current (CC)",
                           "British Columbia (BC)", "Alaska (AK)")),
             size = 5, fill = 'seagreen') +
  labs(title = 'Spatial Stratification OM1 & EM1 (Current Management)')

p2 <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = 'grey22') +
  coord_quickmap() + 
  kaputils::theme_mk(base_size = 16) +  theme(axis.title =element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(-180,-119), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
  scale_y_continuous(expand = c(0,0), limits = c(30,71), breaks = seq(30,70,10), 
                     labels =  paste(seq(30,75,10), "°N"))  +
  geom_hline(yintercept = c(36,50), lwd = 1.1, col = 'goldenrod') +
  geom_vline(xintercept = c(-135,-145), lwd = 1.1, col = 'goldenrod') +
  geom_rect(fill = 'white', aes(xmin = -180, xmax = -135.5,
                                ymin = 30, ymax = 49.8)) +
  geom_label(aes( x = c(rep(-125,3),-140,-155), 
                  y = c(33,40, rep(53,3)),
                  label =paste0('Region ',1:5)), 
             size = 5, fill = 'goldenrod') +
  labs(title = 'Spatial Stratification OM2 & EM2 (demographic breaks)')

ggsave(plot = (p1  | p2),
       file = paste0("./figures/Fig1_strata_maps.png"),
       width = 17, height = 11, units = 'in', dpi = 720)


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
  