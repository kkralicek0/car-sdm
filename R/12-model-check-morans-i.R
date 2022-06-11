###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 12-model-check-morans-i.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Climate change induced shifts in suitable
## habitat for five tree species in the Pacific Northwest projected with 
## spatial-Bayesian hierarchical models' by Karin Kralicek, Jay Ver Hoef, 
## Tara Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - A type of Bayesian model-check.. evaluating model performance here with
##   Moran's I at different distance classes... these summaries are part of 
##   Table 3 in the manuscript (code for AUC & ESS are in script 9...R)
## - Concerning memory, this is a CPU-light / memory-heavy script, so...
##   - only estimating Moran's I for 100 of the 2,000 MCMC samples 
##   - calc with group matrices & tapply (motivated by Wright et al. (2019))
##   - Will only have enough memory for 3 species at a time... so run for 3
##     and then run again for the other 2, saving by species as we go...
##
## About - output:
## - ch2-12-current-residuals-location.rds; ch-12-obsv-dist-mat-true.rds
##   - location: path_output
## - ch2-12-moransi_d_avg.rds
##   - location: path_output
##   - average distance in each class, this object also appears in the rdata
##     file ch2-12-moransi-d_avg.rdata
## - ch2-12-moransi_*.rds, where * is a spp abbreviation (abpr, psmem, etc.)
##   - location: path_output
##   - this is the Moran's I stat for 100 of the MCMC samples (every 20th) for
##     that species; all five rds's will be subsequently (re)loaded and 
##     combined together in this script to be saved as an object in the rdata 
##     file ch2-12-moransi-d_avg.rdata
##   - once all of the ch2-12-moransi_*.rds files exists, will basically just
##     create a list of those objects for ease of storage... so once that obj
##     exists, these species-specific rds objects can be deleted.
## - ch2-12-moransi-d_avg.rdata
##   - location: path_output
##   - all the sp-specific moransi rds and d_avg info
###############################################################################
# data manipulation
library(magrittr)    # for set_rownames()
library(dplyr)       # for ... E V E R Y T H I N G ...
library(purrr)       # working with lists
library(reshape2)    # melt() and dcast()
library(stringr)     # for string help (e.g. str_detect in functions)

# plotting
library(ggplot2)
library(sf)

# parallel
library(parallel)

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_output <- paste0(path_top, "output/")
path_modRDS <- paste0(path_output, "tune fits - rds/")

# path to functions script
path_fn <- paste0(path_top, "scripts - real data/")

# == (Create|Load) residuals and dist mat ================================ ####
# Residuals & distance matrix (for observed FIA plots only)
# - here loading residuals for all MCMC samples... subsetting to 100 out of 
#   2000 samples is done below...
if (file.exists(paste0(path_output, "ch2-12-current-residuals-location.rds")) & 
    file.exists(paste0(path_output, "ch2-12-obsv-dist-mat-true.rds"))) {
  resid_current <- readRDS(paste0(path_output, "ch2-12-current-residuals-location.rds"))
  d_lt <- readRDS(paste0(path_output, "ch2-12-obsv-dist-mat-true.rds"))
} else {
  # residuals for observed FIA plots, for each MCMC (& keep lat/lon for d_lt)
  resid_current <- readRDS(paste0(path_output, "ch2-7-preds-PRISM.rds"))$preds_PRISM %>%
    map(~.x[c("newdat.x", "maps_1k")]) %>%
    map(function(sp.x) {
      # extract observed data
      obsv <- sp.x$newdat.x %>% 
        filter(y_index == 1) %>% # only observed FIA plots for residuals
        select(plot_id, p_a, LAT_ACTUAL, LON_ACTUAL) %>%
        arrange(plot_id) 
      
      # calculate the residuals for each of the 2k draws from the posterior
      # - this also ordered like obsv (i.e. by plot_id)
      # - will return a list of vectors
      resid_mcmc <- sp.x$maps_1k %>% map(function(pred.k) {
          obsv$p_a - (pred.k %>% 
                        filter(plot_id %in% obsv$plot_id) %>%
                        arrange(plot_id) %>%
                        pull(mu_hat))
        })
      
      # return a list of residuals by MCMC sample & a df with plot locations
      list(resid_mcmc = resid_mcmc,
           resid_loc  = obsv %>% select(plot_id, LAT_ACTUAL, LON_ACTUAL))
    })
  
  d_lt <- resid_current %>% 
    map(~.x$resid_loc) %>%
    map(function(sp.x) {
      sp.x <- sp.x %>% 
        # convert to sf object for plot points
        st_as_sf(coords = c("LON_ACTUAL", "LAT_ACTUAL"), crs = "NAD83") %>%
        # distance preserving proj: UTM Zone 11N (EPSG: 32611)
        st_transform(crs = 32611) %>%
        # should still be sorted by plot_id, but be safe...
        arrange(plot_id) %>%
        # big matrix with all the distances... (diag already set to zero)
        st_distance(., .) %>%
        units::set_units("km") %>%
        units::drop_units()
      
      # return the lower triangle of this distance matrix
      sp.x[lower.tri(sp.x)]
    })
  
  resid_current <- resid_current %>% map(~.x$resid_mcmc)
  
  saveRDS(resid_current, 
          file = paste0(path_output, "ch2-12-current-residuals-location.rds"))
  saveRDS(d_lt, 
          file = paste0(path_output, "ch2-12-obsv-dist-mat-true.rds"))
}

# == Moran's I: save for each spp ======================================== ####
# - set up --------------------------------------------------------------- ####
# MCMC samples to calc Moran's I on 
samp_n <- (1:100)*20

# vector of spp names
spp <- c("abpr", "psmem", "qudo", "quga4", "quke")

# list to hold results
moransi <- vector("list", length = length(spp)) %>% setNames(spp)

# function for Moran's I
fn.moranI <- function(resid_mcmc, dist_grp.lt) {
  resid_mcmc %>% map(function(r.k) {
    tapply(outer(r.k - mean(r.k), 
                 r.k - mean(r.k))[lower.tri(outer(r.k - mean(r.k), 
                                                  r.k - mean(r.k)))],
           dist_grp.lt,
           mean,
           na.rm = TRUE) * 
      (length(r.k) / (var(r.k) * (length(r.k) - 1)))
  })
}

# - calc & save by spp (also save d_avg) --------------------------------- ####
# to do this section, need to save by species because this takes up a lot of
# memory to calculate Moran's I... after saving files by species, load them
# back in again and resave as a list of these objects (after this, the original
# files that are named by species can be deleted as they will be redundant...

# Steps 1 & 2:
# - do in two batches, closing R in between... first batch of spp is 
#   {abpr, psmem, qudo} and then do {quga4, quke}. Code below is set to 
#   batch 2 atm. This will require some tweaking of the code below.
if (!file.exists(paste0(path_output, "ch2-12-moransi_quke.rds"))) {
  # - make groups for the 8 distance classes ------------------------ ####
  d_grp <- d_lt %>% map(function(d_lt.x) {
    # make group identifier for each neighbors pair based on the distance class
    # they belong to
    # - note: d_lt.x is the lower triangle of the distance matrix
    # - the distance groups are:
    #   1: (0-5] km
    #   2: (5-10] km
    #   3: (10-15] km
    #   4: (15-25] km
    #   5: (25-50] km
    #   6: (50-100] km
    #   7: (100-200] km
    #   8: (200, ...] km, wherever this maxs out
    
    x <- ceiling(d_lt.x / 5) # make 5 km groups
    x[x == 4 | x == 5]  <- 4 # d_class (15, 25]
    x[x > 5  & x <= 10] <- 5 # d_class (25, 50]
    x[x > 10 & x <= 20] <- 6 # d_class (50, 100]
    x[x > 20 & x <= 40] <- 7 # d_class (100, 200]
    x[x > 40] <- 8 # d_class (100, 200]
    
    # return this lower triangle vector of group membership...
    return(x)
    
  })
  
  # - get avg distance in each class -------------------------------- ####
  d_avg <- list(d_lt, d_grp) %>% pmap(~tapply(.x, .y, mean, na.rm = TRUE))
  
  saveRDS(d_avg, paste0(path_output, "ch2-12-moransi_d_avg.rds"))
  
  # - calculate & save Moran's I, by spp in for-loop ---------------- ####
  for(sp.i in spp) {
    start_time <- Sys.time()
    moransi[[sp.i]] <- fn.moranI(resid_current[[sp.i]][samp_n],
                                 d_grp[[sp.i]])
    print(paste0(sp.i, ": ", Sys.time() - start_time))
    
    # save along the way, just in case
    saveRDS(moransi[[sp.i]],
            file = paste0(path_output, "ch2-12-moransi_", sp.i, ".rds"))
  }
}

# Step 3: 
# - load sp-specific Moran's I rds into a list 
moransi <- spp %>% lapply(function(sp.i) {
  if (file.exists(paste0(path_output, "ch2-12-moransi_", sp.i, ".rds"))) {
    readRDS(paste0(path_output, "ch2-12-moransi_", sp.i, ".rds"))
  }
}) %>% setNames(spp)

# Step 4:
# - save moransi & d_avg for easy manipulation later...
save(moransi, d_avg,
     file = paste0(path_output, "ch2-12-moransi-d_avg.rdata"))

# == Summarize & plot ==================================================== ####
# MCMC samples to calc Moran's I on 
samp_n <- (1:100)*20

# load moransi & d_avg
load(paste0(path_output, "ch2-12-moransi-d_avg.rdata"))

# mashing into a plot-able format... this needs work still
mi_df <- list(moransi, d_avg) %>%
  pmap(function(mi.x, d_avg.x) {
    # atm, mi.x is a list of 100 where each element is a named vector of 
    # length 8 (i.e. length(d_avg) & in the same order as d_avg)... 
    
    data.frame(moransi = mi.x %>% unlist, # length 8 x 100 = 800
               # repeat "20" eight times, then "40" eight times, etc.
               MCMC_sample = rep(samp_n, each = length(d_avg.x)),
               # repeat the average distance (groups 1:8) 100 times
               dist_avg = rep(d_avg.x %>% unname, length(samp_n)))
  })

(mi_df %>% 
   bind_rows(.id = "species") %>% 
   mutate(species = toupper(species))) %>%
  ggplot() +
  geom_line(aes(x = dist_avg, y = moransi, group = MCMC_sample), 
            alpha = .1) +
  facet_grid(species ~ ., scales = "free_y") +
  xlab("distance (thousand km)") +
  ylab("Moran's I") +
  ggtitle("Moran's I by species",
          subtitle = "Based on 100 MCMC samples (every 20th of 2,000)") +
  theme_bw()

# summarize values of Moran's I
mi_df %>% 
  map(function(sp.x) {
  sp.x %>% 
    group_by(dist_avg) %>%
    summarize(mean = mean(moransi),
              max  = max(moransi),
              abs_max = max(abs(moransi))) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}) 

# distance class with greatest absolute mean Moran's I 
mi_df %>% 
  map(function(sp.x) {
    sp.x %>% 
      group_by(dist_avg) %>%
      summarize(mean = mean(moransi),
                max  = max(moransi),
                abs_max = max(abs(moransi))) %>%
      filter(abs(mean) == max(abs(mean))) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }) %>%
  bind_rows(.id = "spp")

