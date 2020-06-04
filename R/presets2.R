## Presets.R
## M S KAPUR 
## Check & load necessary packages and system opts

## CRAN pckgs ----
if(!require("here"))   install.packages("here")
if(!require("tidyverse"))   install.packages("tidyverse") 
if(!require("lubridate"))   install.packages("lubridate") 
if(!require("patchwork"))   install.packages("patchwork") 
if(!require("padr"))   install.packages("padr") 
if(!require("dplyr"))   install.packages("dplyr") 
if(!require("ggplot2"))   install.packages("ggplot2") 
if(!require("reshape2"))   install.packages("reshape2") 
if(!require("TMB"))   install.packages("TMB") 
if(!require("TMBhelper"))   install.packages("TMBhelper") 

## Git pckgs ----
if(!require("PNWColors"))   remotes::install_github("jakelawlor/PNWColors")
if(!require("hrbrthemes"))   remotes::install_github("hrbrmstr/hrbrthemes")
if(!require("r4ss"))   remotes::install_github("r4ss/r4ss")
if(!require("kaputils"))   remotes::install_github("mkapur/kaputils")

options(scipen = 999)
verbose = TRUE 
setwd(here())