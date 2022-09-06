###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 2-add-bioclimatic-vars.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@udsa.gov)
##
## This R script is in support of 'Climate change induced shifts in suitable
## habitat projected for PNW tree species with spatial-Bayesian models' by 
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - add biologically meaningful climate varibles
## - bioclim variables (abrv 'bio**') based on:
##   (a) WorldClim: https://worldclim.org/data/bioclim.html
##   (b) O'Donnell & Ignizio (2012): https://pubs.usgs.gov/ds/691/ds691.pdf
## - growing degree days based on McMaster and Wilhelm (1997)
##
## About - output:
## - ch2-2-data.Rdata
##   - location: path_data (USB)
##   - same contents as output from 1...R, but with subsetted climate vars
##
###############################################################################
library(dplyr)        # for ...E V E R Y T H I N G...
library(stringr)      # for string help (e.g. str_detect in functions)
library(reshape2)     # melt() and dcast()
library(purrr)        # working with lists
library(tidyr)

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")
path_data   <- paste0(path_top, "data/")

# == Load data =========================================================== ####
load(paste0(path_data, "ch2-1-all-vars.Rdata"))
spp <- names(FIA_spp)

# == bioclim - base on USGS calculations ================================= ####
# - (1) subset to climate & add rolling 3-month 'quarter' summaries ------ ####
# identify names to keep
keep_vars <- c("plot_id", "tmeanann", "tmaxann", "tminann", "pptann")
clim_vars <- lapply(c("tmean", "tmax", "tmin", "ppt"), 
                    function(i) {
                      grep(i, names(FIA %>% select(-keep_vars)), 
                           value = TRUE)
                    }) %>% unlist

df <- FIA %>% 
  select(keep_vars, clim_vars) %>%
  group_by(plot_id) %>%
  # abbreviate with r before center month
  mutate(pptr01 = mean(c(ppt12, ppt01, ppt02)),
         pptr02 = mean(c(ppt01, ppt02, ppt03)),
         pptr03 = mean(c(ppt02, ppt03, ppt04)),
         pptr04 = mean(c(ppt03, ppt04, ppt05)),
         pptr05 = mean(c(ppt04, ppt05, ppt06)),
         pptr06 = mean(c(ppt05, ppt06, ppt07)),
         pptr07 = mean(c(ppt06, ppt07, ppt08)),
         pptr08 = mean(c(ppt07, ppt08, ppt09)),
         pptr09 = mean(c(ppt08, ppt09, ppt10)),
         pptr10 = mean(c(ppt09, ppt10, ppt11)),
         pptr11 = mean(c(ppt10, ppt11, ppt12)),
         pptr12 = mean(c(ppt11, ppt12, ppt01))) %>%
  mutate(tmeanr01 = mean(c(tmean12, tmean01, tmean02)),
         tmeanr02 = mean(c(tmean01, tmean02, tmean03)),
         tmeanr03 = mean(c(tmean02, tmean03, tmean04)),
         tmeanr04 = mean(c(tmean03, tmean04, tmean05)),
         tmeanr05 = mean(c(tmean04, tmean05, tmean06)),
         tmeanr06 = mean(c(tmean05, tmean06, tmean07)),
         tmeanr07 = mean(c(tmean06, tmean07, tmean08)),
         tmeanr08 = mean(c(tmean07, tmean08, tmean09)),
         tmeanr09 = mean(c(tmean08, tmean09, tmean10)),
         tmeanr10 = mean(c(tmean09, tmean10, tmean11)),
         tmeanr11 = mean(c(tmean10, tmean11, tmean12)),
         tmeanr12 = mean(c(tmean11, tmean12, tmean01))) %>%
  mutate(tmaxr01 = mean(c(tmax12, tmax01, tmax02)),
         tmaxr02 = mean(c(tmax01, tmax02, tmax03)),
         tmaxr03 = mean(c(tmax02, tmax03, tmax04)),
         tmaxr04 = mean(c(tmax03, tmax04, tmax05)),
         tmaxr05 = mean(c(tmax04, tmax05, tmax06)),
         tmaxr06 = mean(c(tmax05, tmax06, tmax07)),
         tmaxr07 = mean(c(tmax06, tmax07, tmax08)),
         tmaxr08 = mean(c(tmax07, tmax08, tmax09)),
         tmaxr09 = mean(c(tmax08, tmax09, tmax10)),
         tmaxr10 = mean(c(tmax09, tmax10, tmax11)),
         tmaxr11 = mean(c(tmax10, tmax11, tmax12)),
         tmaxr12 = mean(c(tmax11, tmax12, tmax01))) %>%
  mutate(tminr01 = mean(c(tmin12, tmin01, tmin02)),
         tminr02 = mean(c(tmin01, tmin02, tmin03)),
         tminr03 = mean(c(tmin02, tmin03, tmin04)),
         tminr04 = mean(c(tmin03, tmin04, tmin05)),
         tminr05 = mean(c(tmin04, tmin05, tmin06)),
         tminr06 = mean(c(tmin05, tmin06, tmin07)),
         tminr07 = mean(c(tmin06, tmin07, tmin08)),
         tminr08 = mean(c(tmin07, tmin08, tmin09)),
         tminr09 = mean(c(tmin08, tmin09, tmin10)),
         tminr10 = mean(c(tmin09, tmin10, tmin11)),
         tminr11 = mean(c(tmin10, tmin11, tmin12)),
         tminr12 = mean(c(tmin11, tmin12, tmin01))) %>%
  ungroup()

# - (2) reformat data & summarize bioclim variables ---------------------- ####
# Notes:
# - following the common notation of 'bio**' etc.
# - if a tie exists in seasonal or monthly summaries, we're following the USGS
#   methods that take the first chronologically occurring stat

# grab the low-hanging fruits first...
bioclim_vars <- df %>% select(plot_id, tmeanann, pptann) %>%
  rename(bio1 = tmeanann, # mean annual temp
         bio12 = pptann) # total annual precip

# reformat data for ease of calculation...
df <- df %>% 
  melt(id.vars = "plot_id") %>%
  mutate(var2 = variable) %>%
  separate(col = var2, into = c("stats", "month_season"), sep = -2) %>%
  dcast(plot_id + month_season ~ stats) %>%
  group_by(plot_id)

# get the rest of the variables!
bioclim_vars <- df %>% 
  summarize(bio2 = mean(tmax - tmin, na.rm = TRUE), # mean diurnal range
            bio4 = sd(tmean, na.rm = TRUE), # temp seasonality
            bio5 = max(tmax, na.rm = TRUE), # max temp of warmest month
            bio6 = min(tmin, na.rm = TRUE), # min temp of coldest month
            bio15 = sd(ppt, na.rm = TRUE) / # precipitation seasonality
              mean(ppt, na.rm = TRUE)) %>%
  mutate(bio7 = bio5 - bio6) %>% # temp annual range
  inner_join(bioclim_vars, by = "plot_id")

bioclim_vars <- df %>%
  filter(pptr == max(pptr, na.rm = TRUE)) %>%
  summarize(bio8 = tmeanr[1], # mean temp of wettest quarter
            bio16 = pptr[1]) %>% # precipitation in wettest quarter
  inner_join(bioclim_vars, by = "plot_id")

bioclim_vars <- df %>%
  filter(pptr == min(pptr, na.rm = TRUE)) %>%
  summarize(bio9 = tmeanr[1], # mean temp of driest quarter
            bio17 = pptr[1]) %>% # precipitation in driest quarter
  inner_join(bioclim_vars, by = "plot_id")

bioclim_vars <- df %>%
  filter(tmeanr == max(tmeanr, na.rm = TRUE)) %>%
  summarize(bio10 = tmeanr[1], # mean temp of warmest quarter
            bio18 = pptr[1]) %>% # precipitation in warmest quarter
  inner_join(bioclim_vars, by = "plot_id")

bioclim_vars <- df %>%
  filter(tmeanr == min(tmeanr, na.rm = TRUE)) %>%
  summarize(bio11 = tmeanr[1], # mean temp of coldest quarter
            bio19 = pptr[1]) %>% # precipitation in coldest quarter
  inner_join(bioclim_vars, by = "plot_id")

bioclim_vars <- df %>%
  filter(ppt == max(ppt, na.rm = TRUE)) %>%
  summarize(bio13 = ppt[1]) %>% # precipitation in wettest month
  inner_join(bioclim_vars, by = "plot_id")

bioclim_vars <- df %>%
  filter(ppt == min(ppt, na.rm = TRUE)) %>%
  summarize(bio14 = ppt[1]) %>% # precipitation in driest month
  inner_join(bioclim_vars, by = "plot_id")

# - (3) FIA & FIA_spp: add bioclim variables ----------------------------- ####
FIA <- FIA %>% inner_join(bioclim_vars, by = "plot_id")
FIA_spp <- FIA_spp %>% map(~.x %>% inner_join(bioclim_vars, by = "plot_id"))

# clean-up space
rm(df, keep_vars, bioclim_vars)

# == Growing deg days ==================================================== ####
# Sum the month-days where the tmean exceeds 5C (frost-free days)
# - (1) key for no. days / month ----------------------------------------- ####
day_key <- data.frame(month = paste0("tmean", 
                                     c(rep(0, 9), rep(1, 3)),
                                     c(1:9, 0:2)),
                      n_days = c(31, 28, 31, 30, 31, 30, 
                                 31, 31, 30, 31, 30, 31))
# - (2) reformat data & summarize gdd5 ----------------------------------- ####
df <- FIA %>% 
  select(plot_id, 
         grep("tmean", names(FIA %>% select(-tmeanann)), value = TRUE)) %>%
  melt(id.vars = "plot_id",
       variable.name = "month") %>%
  inner_join(day_key, by = "month") %>%
  group_by(plot_id) %>%
  filter(value >= 5) %>%
  summarize(gdd5 = sum(n_days)) 

# - (3) FIA & FIA_spp: add growing degree-day (5C) variable -------------- ####
FIA <- FIA %>% inner_join(df, by = "plot_id")
FIA_spp <- FIA_spp %>% map(~.x %>% inner_join(df, by = "plot_id"))

# clean up space 
rm(df, day_key)
  

# == heat-moisture index ================================================= ####
# add annual heat-moisture index variables
# - form of the heat-moisture index motivated by ClimateNA variable
FIA <- FIA %>%
  group_by(plot_id) %>%
  mutate(hmi = c(tmeanann + 10)/(pptann/1000)) %>%
  ungroup
FIA_spp <- FIA_spp %>% map(function(sp.x) {
  sp.x %>% inner_join(FIA %>% select(plot_id, hmi), 
                      by = "plot_id")
})

# == FIA, FIA_spp: keep only pot. modeling climate variables ============= ####
keep_vars <- c(grep("bio", names(FIA), value = TRUE), "gdd5", "hmi")
keep_vars <- c("plot_id", "fiahex_id", "LAT_ACTUAL", "LON_ACTUAL",
               "hex_centroid.x_UTM11N_ACTUAL", "hex_centroid.y_UTM11N_ACTUAL",
               "plot_obsv", "fiahex_obsv", "n_plots_hex", "PLT_CN",
               "ES", keep_vars)

FIA <- FIA %>% select(keep_vars, spp)
FIA_spp <- FIA_spp %>% map(~.x %>% select(keep_vars, p_a))

# == Export as Rdata ===================================================== ####
save(sppcd, spp_names, ES_key, plot_key,  # reference dfs
     FIA, FIA_spp,                        # data (Y, r.e., W, U)
     file = paste0(path_data, "ch2-2-data.Rdata"))

