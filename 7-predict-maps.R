###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 7-predict-maps.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Climate change induced shifts in suitable
## habitat projected for PNW tree species with spatial-Bayesian models' by 
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - generate prediction maps for different scenarios.
## - for each prediction data set, 2k maps will be produced that correspond to 
##   2k MCMC samples kept after burn-in. The prediction data sets are:
##   - current climate, via Norm81m PRISM
##   - future climate, via ClimateNA with AOGCM + pathway combinations:
##     - CCSM4 + {RCP4.5, RCP8.5}
##     - HADGEM2-ES + {RCP4.5, RCP8.5}
##
## ASIDE - ch3 & ch2 relation:
## - ch3 work will also be based on these results, but (atm) will add another 
##   future downscaled source (MACAv2-LIVNEH, in addition to ClimateNA) using 
##   the same (AOGCM + pathway + period) combinations above
##
## About - output:
## - ch2-7-predictions-PRISM.rds
##   - location: path_output
##   - current predictions (based on climate data use to fit the models);
##   - both 2k maps and mean-prediction map
## - ch2-7-predictions-ClimateNA.rds
##   - location: path_output
##   - future predictions (based on ClimateNA data)
##     - (2 aogcm) x (2 pathways) x (3 time periods)
##   - both 2k maps & mean-prediction
## - ch2-7-drop_plot_ids.rds
##   - location: path_ouput
##   - key identifying 29 plots that CNA dataset doesn't cover
###############################################################################
library(magrittr)    # for set_rownames()
library(dplyr)       # for ... E V E R Y T H I N G ...
library(purrr)       # working with lists
library(reshape2)    # melt() and dcast()
library(stringr)     # for string help (e.g. str_detect in functions)

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_output <- paste0(path_top, "output/")
path_modRDS <- paste0(path_output, "tune fits - rds/")

# path to functions script
path_fn <- paste0(path_top, "scripts - real data/")

# future climate data - ClimateNA
path_prodataCNA <- paste0(path_USB,
                          "/Karin phd thesis files/r - processed master data/",
                          "future climate data/processed data - ClimateNA/")

# == Load - final MCMC run & fns ========================================= ####
# load final MCMC samples/run
run_name <- "final_MCMC_samples"
mod_spp <- readRDS(paste0(path_output, "ch2-6-final-MCMC-samples.rds")) 

# - extract spp names
spp <- names(mod_spp)

# fns
# - this also in 0-...R but no need to load all the other fns...
fn.inv_logit <- function(x) exp(x) / (1 + exp(x))

# == Load - PRISM (current climate) ====================================== ####
# formatted & standardized (i.e. ready to predict)
dat_ls <- lapply(spp, function(x) {
  readRDS(paste0(path_output, "ch2-5-sampler-ready-data_", x, ".rds"))
}) %>% setNames(spp)

# unformated, contains info on how to standardize future climate data
load(paste0(path_data, "ch2-4-modeling-data.Rdata"))

# == Load - ClimateNA (future climate) =================================== ####
# future scenarios of interest
future_scenarios <- c("CCSM4_rcp45_m2085","CCSM4_rcp85_m2085",
                      "HadGEM2_ES_rcp45_m2085", "HadGEM2_ES_rcp85_m2085")

# read in processed ClimateNA data & add label for scenario to data
CNA_bio <- lapply(future_scenarios, function(x) {
  # read in each {aogcm x pathway x period} data set...
  readRDS(paste0(path_prodataCNA, "ch2-ClimateNA-processed_", x, ".rds")) %>%
    # add a label for this scenario to the data.frame
    mutate(scenario = x)
}) %>% setNames(future_scenarios)

# == Create key for 29 plots to drop {qudo, quga4} ======================= ####
# About:
# - ClimateNA says 29 plots are outside of its modeling extent and therefore 
#   we don't have future climate data for these sites 
# - These 29 occur off the coast of CA (islands) and are all non-forest (only
#   impacts qudo & quga4 datasets)

drop_plot_ids <- CNA_bio[[1]] # same across all CNA scenarios
drop_plot_ids <- drop_plot_ids[!complete.cases(drop_plot_ids),]$plot_id

# save key for easy reference
saveRDS(drop_plot_ids,
        file = paste0(path_output, "ch2-7-drop_plot_ids.rds"))

# == Predict on 2k MCMC samples (current & future) ======================= ####
# About pps_pa:
# - posterior predicted simulations (pps) of presence & absence (pa)
#   instead of thresholding, we're using our prediction of mu_hat to 
#   simulate species presence or absence at plot locations. pps_pa can then
#   be used to calculate range size (via BKDE or Voronoi-polygon methods), 
#   and ANN dist metrics... 
#
# Notes to self on previous versions of this script:
# - my objects might have maps_1k instead of being labeled 2k, b/c just 
#   missed correcting this originally... below is corrected.
# - in my objects, ignore se_lp. Originally, (prior to 03/13/21) instead of 
#   back-transforming se_me (standard error on the scale of the mean (prob; 
#   mu_hat)) to the scale of the linear predictor, had incorrectly calculated 
#   se_lp as the SE of lp_hat -- this is incorrect b/c not a linear operation. 
#   Take away: don't use se_lp in my output from this script (although 
#   corrected so that this isn't an issue in the code below)
#
# - Set seed ------------------------------------------------------------- ####
# doing this for for pps_pa
set.seed(42)

# - predict & save: PRISM ------------------------------------------------ ####
# get predictions
start_time <- Sys.time()
preds_PRISM <- list(mod_spp, dat_ls) %>% 
  # get predictions...
  pmap(function(mod.x, dat.x) {
  # = Grab data ========================================================== ####
  newdat.x <- dat.x$moddat
  
  # = Predict: 2k maps (and drop plots) ================================== ####
  # predictions for 2k maps (use parameter samples from a given iteration)
  maps_2k <- lapply(1:nrow(mod.x$MCMC_theta), function(iter) {
    # vector with correct z repeated (lines up fiahex_id's with plot_id's)
    # - re-assemble z vector from z_ob- and z_un-blocks; order matches z_index...
    z <- mod.x$MCMC_z[iter, ] %>% map(~.x %>% unlist) %>% unlist
    Uz <- as.vector(dat.x$U %*% z)
    
    # vector with correct v repeated (lines up ES's with plot_id's)
    # - again, order matches Q's cols...
    v  <- mod.x$MCMC_v[iter, ] %>% map(~.x %>% unlist) %>% unlist
    Qv <- as.vector(dat.x$Q %*% v)
    
    # betas - note beta_0 always first...
    parms <- mod.x$MCMC_theta[iter, ]
    p <- sum(grepl("beta_", names(parms)))
    XB <- parms$beta_0 + as.matrix(newdat.x[,1:(p-1)]) %*% unlist(parms[2:p])
    
    # predictions
    lp_hat <- XB + Qv + Uz
    mu_hat <- fn.inv_logit(lp_hat)
    
    # format data to return (& drop 29 plots)
    df_out <- data.frame(plot_id = newdat.x$plot_id,
                         lp_hat = lp_hat,
                         mu_hat = mu_hat,
                         scenario = "current",
                         stringsAsFactors = FALSE)
    df_out <- df_out %>% filter(!plot_id %in% drop_plot_ids)
      
    return(df_out)
  })
  
  # = Predict: 'mean-pred' map (and drop plots) ========================== ####
  # for the scenario, take the mean of the 2k predictions
  # - get 'mean' of samples ---------------------------------------------- ####
  mean_theta <- mod.x$MCMC_theta %>% 
    apply(2, mean) %>%
    bind_rows %>%
    melt(variable.name = "theta") %>% 
    filter(theta != "n_sim")
  
  mean_z <- mod.x$MCMC_z %>% 
    apply(2, mean) %>% 
    bind_rows %>% 
    melt(variable.name = "fiahex_id") %>%
    mutate(fiahex_id = gsub("fiahex_id", "", fiahex_id))
  
  mean_v <- mod.x$MCMC_v %>% 
    apply(2, mean) %>% 
    bind_rows %>% 
    melt(variable.name = "ES") %>%
    mutate(ES = gsub("ES", "", ES))
  
  # - get mean-predictions ----------------------------------------------- ####
  # vector with correct z repeated (line up fiahex_id's with plot_id's)
  # - re-assemble z vector from z_ob- and z_un-blocks; order matches z_index...
  z <- mean_z$value
  Uz <- as.vector(dat.x$U %*% z)
  
  # vector with correct v repeated (line up ES's with plot_id's)
  # - again, order matches Q's cols...
  v  <- mean_v$value
  Qv <- as.vector(dat.x$Q %*% v)
  
  # betas - note beta_0 always first...
  parms <- mean_theta
  p <- sum(grepl("beta_", parms$theta))
  XB <- parms$value[1] + as.matrix(newdat.x[,1:(p-1)]) %*% unlist(parms$value[2:p])
  
  # predictions
  lp_hat <- XB + Qv + Uz
  mu_hat <- fn.inv_logit(lp_hat)
  
  # format data to return
  maps_mean <- data.frame(plot_id = newdat.x$plot_id,
                          lp_hat = lp_hat,
                          mu_hat = mu_hat,
                          scenario = "current",
                          stringsAsFactors = FALSE)

    # clean up space
  rm(z, Uz, v, Qv, parms, p, XB, lp_hat, mu_hat)
  
  # - drop the 29 plots (28 really) -------------------------------------- ####
  maps_mean <- maps_mean %>% filter(!plot_id %in% drop_plot_ids)
  newdat.x  <- newdat.x %>% filter(!plot_id %in% drop_plot_ids)
  
  # Unless keep plots associated with ES and fiahex_ids, drop those associated
  # with these drop plots
  mean_z <- mean_z %>% filter(fiahex_id %in% newdat.x$fiahex_id)
  mean_v <- mean_v %>% filter(ES %in% newdat.x$ES)
  
  # = return predictions ================================================= ####
  return(list(newdat.x = newdat.x %>% filter(!plot_id %in% drop_plot_ids),
              maps_2k = maps_2k, 
              maps_mean = maps_mean,
              maps_mean.etc = list(mean_theta = mean_theta,
                                   mean_z = mean_z,
                                   mean_v = mean_v)))
}) %>%
  # add pps_pa & se vars...
  map(function(preds.x) {
    # n obsv for pps_pa's rbinom call
    n.x <- nrow(preds.x$maps_mean) 
    
    # 'mean-pred' map
    # - add pps_pa
    preds.x$maps_mean$pps_pa <- rbinom(n.x, 1, preds.x$maps_mean$mu_hat)
    # - add se_mu based on 2k maps to 'mean-pred' map df
    #   (kk: this is where in previous versions I had incorrectly calc 
    #    se_lp as sd(lp_hat) -- so removed!)
    preds.x$maps_mean <- preds.x$maps_mean %>%
      inner_join(preds.x$maps_2k %>%
                   bind_rows %>%
                   group_by(plot_id) %>%
                   summarize(se_mu = sd(mu_hat)) %>% 
                   ungroup,
                 by = "plot_id")
    
    # 2k maps (MCMC) 
    # - add pps_pa to each run
    preds.x$maps_2k <- preds.x$maps_2k %>% map(function(df.x) {
      df.x$pps_pa <- rbinom(n.x, 1, df.x$mu_hat)
      return(df.x)
    })
    # - add MCMC_samp_n
    preds.x$maps_2k <- list(preds.x$maps_2k, c(1:length(preds.x$maps_2k))) %>%
      pmap(function(df.x, x) {
        df.x <- df.x %>% mutate(MCMC_samp_n = x)
      })
    
    return(preds.x)
  })
run_time <- Sys.time() - start_time
run_time

# save predictions
saveRDS(list(preds_PRISM = preds_PRISM,
             run_name = run_name,
             run_time = run_time), 
        file = paste0(path_output, "ch2-7-preds-PRISM.rds"))

# - Set seed ------------------------------------------------------------- ####
# doing this for for pps_pa
set.seed(396)

# - predict & save: ClimateNA -------------------------------------------- ####
# get predictions & reformat so org at top level by spp, then scenario
start_time <- Sys.time()
preds_CNA <- CNA_bio %>% map(function(newdat.x) {
  # work with an individual scenario ({aogcm x pathway x period} combo)
  list(mod_spp, dat_ls, data_spp.std) %>% 
    # get predictions...
    pmap(function(mod.x, dat.x, std.x) {
      # = Order covars (X) for this species ============================== ####
      # - subset to plots & vars of interest & std vars ------------------ ####
      clim.x <- std.x$mean_sd$clim_var %>% as.vector # kk8/17: had to add this `as.vector()` too...
      scenario.x <- newdat.x$scenario[1]
      newdat.x <- newdat.x %>% 
        filter(plot_id %in% dat.x$moddat$plot_id) %>% 
        select(plot_id, one_of(clim.x)) # kk8/17: older dplyr did not have `any_of()`
      
      # standardize prediction data by same standardization as training data
      # - (i.e. use the same mean & sd used to std current climate data)
      newdat.x <- newdat.x %>%
        melt(id.vars = "plot_id", 
             variable.name = "clim_var") %>%
        inner_join(std.x$mean_sd, by = "clim_var") %>%
        mutate(value_std = (value - std_mean) / std_sd) %>%
        select(plot_id, clim_var, value_std)  %>% 
        dcast(plot_id ~ clim_var, value.var = "value_std")
      
      # - get correct variable forms (lin/quad/interactions) -------- ####
      # (1) identify variables to manipulate ----
      # variables to square
      var_quad <- grep("bio", clim.x, value = TRUE)
      var_quad.name <- paste0(var_quad, "_2")
      
      # variables to interact
      var_inter <- list(("bio10" %in% clim.x) & ("bio18" %in% clim.x),
                        ("bio8" %in% clim.x)  & ("bio16" %in% clim.x),
                        ("bio8" %in% clim.x)  & ("bio13" %in% clim.x))
      var_inter <- list(switch(var_inter[[1]] + 1, NA, c("bio10", "bio18")),
                        switch(var_inter[[2]] + 1, NA, c("bio8", "bio16")),
                        switch(var_inter[[3]] + 1, NA, c("bio8", "bio13")))
      var_inter.name <- list("bio10_bio18",
                             "bio8_bio16",
                             "bio8_bio13")[!is.na(var_inter)]
      var_inter <- var_inter[!is.na(var_inter)]
      
      # (2) format data: explicitly add new variable forms ----
      # add quadratic terms
      for(i in 1:length(var_quad)) {
        varname <- var_quad.name[i]
        var_lin <- var_quad[i]
        
        newdat.x <- newdat.x %>% 
          mutate(!!varname := (!!as.name(var_lin))^2)
      }
      
      # add interaction terms
      for(i in 1:length(var_inter)) {
        varname  <- var_inter.name[[i]]
        var_lin1 <- var_inter[[i]][1]
        var_lin2 <- var_inter[[i]][2]
        
        newdat.x <- newdat.x %>% 
          mutate(!!varname := (!!as.name(var_lin1)) * (!!as.name(var_lin2)))
      }
      
      # (3) order variables in df ----
      # doing this to help with getting correct starting parameter values assoc.
      # with the right variables...
      # kk8/17: had to change this line too.. was `select(starts_with(c("gdd5", "hmi", "bio")), everything())`
      newdat.x <- suppressWarnings(newdat.x %>% 
                                     select(one_of(c("gdd5", "hmi")), 
                                            starts_with("bio"), 
                                            everything())) 
      
      # = Predict: 2k maps (and drop plots) ============================== ####
      # predictions for 2k maps (use parameter samples from a given iteration)
      maps_2k <- lapply(1:nrow(mod.x$MCMC_theta), function(iter) {
        # vector with correct z repeated (line up fiahex_id's with plot_id's)
        # - re-assemble z vector from z_ob- and z_un-blocks; order matches z_index...
        z <- mod.x$MCMC_z[iter, ] %>% map(~.x %>% unlist) %>% unlist
        Uz <- as.vector(dat.x$U %*% z)
        
        # vector with correct v repeated (line up ES's with plot_id's)
        # - again, order matches Q's cols...
        v  <- mod.x$MCMC_v[iter, ] %>% map(~.x %>% unlist) %>% unlist
        Qv <- as.vector(dat.x$Q %*% v)
        
        # betas - note beta_0 always first...
        parms <- mod.x$MCMC_theta[iter, ]
        p <- sum(grepl("beta_", names(parms)))
        XB <- parms$beta_0 + as.matrix(newdat.x[,1:(p-1)]) %*% unlist(parms[2:p])
        
        # predictions
        lp_hat <- XB + Qv + Uz
        mu_hat <- fn.inv_logit(lp_hat)
        
        # format data to return (& drop 29 plots)
        df_out <- data.frame(plot_id = newdat.x$plot_id,
                             lp_hat = lp_hat,
                             mu_hat = mu_hat,
                             scenario = scenario.x,
                             stringsAsFactors = FALSE)
        df_out <- df_out %>% filter(!plot_id %in% drop_plot_ids)
        
        return(df_out)
        
      })
      
      # = Predict: 'mean-pred' map (and drop plots) ======================= ####
      # predictions for the 'mean-pred' scenario (use mean of samples)
      # - get 'mean' of samples ------------------------------------------ ####
      mean_theta <- mod.x$MCMC_theta %>% 
        apply(2, mean) %>%
        bind_rows %>%
        melt(variable.name = "theta") %>% 
        filter(theta != "n_sim")
      
      mean_z <- mod.x$MCMC_z %>% 
        apply(2, mean) %>% 
        bind_rows %>% 
        melt(variable.name = "fiahex_id") %>%
        mutate(fiahex_id = gsub("fiahex_id", "", fiahex_id))
      
      mean_v <- mod.x$MCMC_v %>% 
        apply(2, mean) %>% 
        bind_rows %>% 
        melt(variable.name = "ES") %>%
        mutate(ES = gsub("ES", "", ES))
      
      # - get mean-predictions ------------------------------------------- ####
      # vector with correct z repeated (line up fiahex_id's with plot_id's)
      # - re-assemble z vector from z_ob- and z_un-blocks; order matches z_index...
      z <- mean_z$value
      Uz <- as.vector(dat.x$U %*% z)
      
      # vector with correct v repeated (line up ES's with plot_id's)
      # - again, order matches Q's cols...
      v  <- mean_v$value
      Qv <- as.vector(dat.x$Q %*% v)
      
      # betas - note beta_0 always first...
      parms <- mean_theta
      p <- sum(grepl("beta_", parms$theta))
      XB <- parms$value[1] + as.matrix(newdat.x[,1:(p-1)]) %*% unlist(parms$value[2:p])
      
      # predictions
      lp_hat <- XB + Qv + Uz
      mu_hat <- fn.inv_logit(lp_hat)
      
      # format data to return
      maps_mean <- data.frame(plot_id = newdat.x$plot_id,
                              lp_hat = lp_hat,
                              mu_hat = mu_hat,
                              scenario = scenario.x,
                              stringsAsFactors = FALSE)
      
      # clean up space
      rm(z, Uz, v, Qv, parms, p, XB, lp_hat, mu_hat)
      
      # - drop the 29 plots (28 really) ---------------------------------- ####
      maps_mean <- maps_mean %>% filter(!plot_id %in% drop_plot_ids)
      newdat.x  <- newdat.x %>% filter(!plot_id %in% drop_plot_ids)
      
      # Unless keep plots associated with ES and fiahex_ids, drop those associated
      # with these drop plots
      mean_z <- mean_z %>% filter(fiahex_id %in% newdat.x$fiahex_id)
      mean_v <- mean_v %>% filter(ES %in% newdat.x$ES)
      
      # = return predictions ============================================= ####
      return(list(newdat.x = newdat.x,
                  maps_2k = maps_2k,
                  maps_mean = maps_mean,
                  maps_mean.etc = list(mean_theta = mean_theta,
                                       mean_z = mean_z,
                                       mean_v = mean_v)))
    }) %>%
    # add pps_pa & se vars...
    map(function(preds.x) {
      # n obsv for pps_pa's rbinom call
      n.x <- nrow(preds.x$maps_mean)
      
      # 'mean-pred' map
      # - add pps_pa
      preds.x$maps_mean$pps_pa <- rbinom(n.x, 1, preds.x$maps_mean$mu_hat)
      # - add se_mu based on 2k maps to 'mean-pred' map df
      #   (kk: this is where in previous versions I had incorrectly calc 
      #    se_lp as sd(lp_hat) -- so removed!)
      preds.x$maps_mean <- preds.x$maps_mean %>%
        inner_join(preds.x$maps_2k %>%
                     bind_rows %>%
                     group_by(plot_id) %>%
                     summarize(se_mu = sd(mu_hat)) %>%
                     ungroup,
                   by = "plot_id")
      
      # 2k maps (MCMC)
      # - add pps_pa to each run
      preds.x$maps_2k <- preds.x$maps_2k %>% map(function(df.x) {
        df.x$pps_pa <- rbinom(n.x, 1, df.x$mu_hat)
        return(df.x)
      })
      # - add MCMC_samp_n
      preds.x$maps_2k <- list(preds.x$maps_2k, c(1:length(preds.x$maps_2k))) %>%
        pmap(function(df.x, x) {
          df.x <- df.x %>% mutate(MCMC_samp_n = x)
        })
      
      return(preds.x)
    })
})
run_time <- Sys.time() - start_time
run_time

# save predictions!
lapply(future_scenarios, function(x) {
  saveRDS(list(preds_CNA = preds_CNA[[x]],
               run_name = run_name,
               run_time = run_time),
          file = paste0(path_output, "ch2-7-preds-CNA_", x, ".rds"))
})

# save future prediction names (for ease in loading later on)
saveRDS(future_scenarios,
        file = paste0(path_output, "ch2-7-future_scenarios_chr.rds"))



