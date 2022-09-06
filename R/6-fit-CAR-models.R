###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 6-fit-CAR-models.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Climate change induced shifts in suitable
## habitat projected for PNW tree species with spatial-Bayesian models' by 
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
##
## About - this script:
## - tune CAR SDM sampler for the spp
## - Note on backward use:
##   - to examine the organization of the spatial random effects as tuning
##     progresses, this script loads some mapping objects created in script 
##     8-voronoi-polygons.R; these objects do not rely on script 6 & 7 to be
##     generated, but made more sense to include in that script than here).
## - Appendix B, Table B1: the very end of this script has the code for the
##   numbers in Table B1.
##
## About - output:
## - ch2-6-final-MCMC-samples.rds
##   - location: path_output
##   - final MCMC samples 
## - ch2-6-final-MCMC-burnin.rds
##   - location: path_output
##   - burn-in for the finalMCMC samples 
## - mod_t*.rds
##   - location: path_mod
##   - sampler progress where * numbers are sequential run-segments 
##  
###############################################################################
# data manipulation
library(magrittr)    # for set_rownames()
library(dplyr)       # for ... E V E R Y T H I N G ...
library(purrr)       # working with lists
library(reshape2)    # melt() and dcast()
library(stringr)     # for string help (e.g. str_detect in functions)

# plotting
library(ggplot2)

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_output <- paste0(path_top, "output/")
path_modRDS <- paste0(path_output, "tune fits - rds/")

# path to functions script
path_fn <- paste0(path_top, "scripts - real data/")

# == Load data & fns ===================================================== ####
# do names manually...
spp <- c("abpr", "psmem", "qudo", "quga4", "quke")
dat_ls <- lapply(spp, function(sp.x) {
  x <- readRDS(paste0(path_output, "ch2-5-sampler-ready-data_", sp.x, ".rds"))
  y <- readRDS(paste0(path_output, "ch2-5-rholookup-table_", sp.x, ".rds"))
  x <- x %>% list_modify("rho_lookup" = y$rho_lookup)
  return(x)
}) %>% setNames(spp)

source(paste0(path_fn, "0-functions.R"))

# == Define sampler-helper fns =========================================== ####
# grab new starting values from previous model run 
# (e.g. start at the end of where that chain left off...)
fn.parms_restart <- function(mod.x, dat.x, nMCMC) {
  parms_names <- dat.x$parms_names
  z_index.x   <- dat.x$W_ls$z_index
  n_theta     <- sum(!parms_names %in% c("z_ob", "z_un", "v"))
  
  parms_restart <- vector("list", length(parms_names)) %>% setNames(parms_names)
  parms_restart[1:n_theta] <- mod.x$MCMC_theta[nMCMC, ] %>% select(-n_sim) %>% as.list
  parms_restart[["z_ob"]]  <- mod.x$MCMC_z[nMCMC, as.logical(z_index.x)] %>% unlist
  parms_restart[["z_un"]]  <- mod.x$MCMC_z[nMCMC, !as.logical(z_index.x)] %>% unlist
  parms_restart[["v"]] <- mod.x$MCMC_v[nMCMC, ] %>% unlist
  
  dat.x$parms_start <- parms_restart
  
  return(dat.x)
}

# parallel wrapper for fn.CAR_fit (from 0-...R script)
fn_par_wrapper.CAR_fit <- function(dat_ls, nMCMC, nMCMC_z, nThin, 
                                   n_cores = NULL) {
  # - set up (parallel) -------------------------------------------------- ####
  if(is.null(n_cores)) {
    n_cores <- min(max(1, (detectCores() - 1)), length(dat_ls))
  }
  cl <- makeCluster(n_cores)               # initiate cluster
  clusterEvalQ(cl, {  
    library(Matrix) 
    library(dplyr)
    library(purrr)
    NULL # don't print out packages loaded on each cluster...
  })
  clusterExport(cl, # load global enviro vars to nodes
                varlist = c("fn.CAR_fit", "fn.llike", "fn.draw", 
                            "fn.inv_logit","nMCMC", "nMCMC_z", "nThin"))
  
  # - run sampler -------------------------------------------------------- ####
  mod_ls <- parLapply(cl, dat_ls, function(x) {
    fn.CAR_fit(x$moddat, x$W_ls, x$U, x$Q, x$rho_lookup, 
               x$parms_names, x$parms_start, x$tune_start, 
               nMCMC, nMCMC_z, nThin)
  })
  
  # - shut down (parallel) ----------------------------------------------- ####
  stopCluster(cl) # close cluster
  
  # - export the run ----------------------------------------------------- ####
  return(mod_ls)
}

# == Start sampler ======================================================= ####
# grab number of betas in each data set for ease of tune-updates
n_betas <- dat_ls %>% map(~length(grep("beta_", .x$parms_names)))

# - previous runs ... minimize here to save on space... ------------------ ####
# - MCMC: no thin  + 1k samps + 20:1 z's {1:6} --------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 1

# - tune 0 (e.g. mod_t0) ------------------------------------------------- ####mm
# adjust tuning 
# - (initially setting all to the same values)
dat_ls[["abpr"]]$tune_start  <- c(rep(.1, n_betas$abpr),  10, .01, .1, .005, .05)
dat_ls[["psmem"]]$tune_start <- c(rep(.1, n_betas$psmem), 10, .01, .1, .005, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.1, n_betas$qudo),  10, .01, .1, .005, .005, .05)
dat_ls[["quga4"]]$tune_start <- c(rep(.1, n_betas$quga4), 10, .01, .1, .005, .005, .05)
dat_ls[["quke"]]$tune_start  <- c(rep(.1, n_betas$quke),  10, .01, .1, .005, .05)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# run sampler (parallel: 5 cores ~ 10 mins)
start_time <- Sys.time()
mod_t0 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t0, file = paste0(path_modRDS, "mod_t0.rds"))

# - tune 1 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.10, n_betas$abpr),  12, .01, .2, .004, .06)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 12, .01, .2, .004, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.15, n_betas$qudo),  12, .01, .2, .005, .005, .06)
dat_ls[["quga4"]]$tune_start <- c(rep(.15, n_betas$quga4), 12, .01, .2, .004, .004, .06)
dat_ls[["quke"]]$tune_start  <- c(rep(.10, n_betas$quke),  12, .01, .2, .005, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t0, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t1 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t1, file = paste0(path_modRDS, "mod_t1.rds"))

# - tune 2 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.15, n_betas$abpr),  14, .01, .4, .002, .08)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 14, .01, .4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.20, n_betas$qudo),  14, .01, .4, .005, .005, .08)
dat_ls[["quga4"]]$tune_start <- c(rep(.20, n_betas$quga4), 14, .01, .4, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  14, .01, .4, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t1, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t2 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t2, file = paste0(path_modRDS, "mod_t2.rds"))

# - tune 3 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.14, n_betas$abpr),  14, .01, .8, .003, .09)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 16, .01, .6, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.30, n_betas$qudo),  15, .01, .6, .005, .005, .10)
dat_ls[["quga4"]]$tune_start <- c(rep(.30, n_betas$quga4), 16, .01, .6, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  16, .01, .6, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t2, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t3 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t3, file = paste0(path_modRDS, "mod_t3.rds"))

# - tune 4 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.14, n_betas$abpr),  16, .01, .8, .003, .09)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 18, .01, .8, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.35, n_betas$qudo),  17, .01, .8, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 18, .01, .8, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  16, .01, .6, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t3, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t4 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t4, file = paste0(path_modRDS, "mod_t4.rds"))

# - tune 5 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.14, n_betas$abpr),  18, .01, 0.8, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 20, .01, 1.0, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.40, n_betas$qudo),  19, .01, 1.0, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 20, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  15, .01, 0.6, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t4, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t5 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t5, file = paste0(path_modRDS, "mod_t5.rds"))

# - tune 6 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.14, n_betas$abpr),  20, .01, 0.8, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 22, .01, 1.2, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  21, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 20, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  15, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t5, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
# mod_t6 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t6, file = paste0(path_modRDS, "mod_t6.rds"))

# - MCMC: thin 50% + ... {7:9} ------------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 2 # i.e. keep 1 out of every two...

# - tune 7 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.13, n_betas$abpr),  22, .01, 1.0, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 24, .01, 1.2, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  23, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 22, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  13, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t6, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t7 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t7, file = paste0(path_modRDS, "mod_t7.rds"))

# - tune 8 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.13, n_betas$abpr),  22, .01, 1.0, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 24, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  25, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 24, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  10, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t7, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t8 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t8, file = paste0(path_modRDS, "mod_t8.rds"))

# - tune 9 --------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.13, n_betas$abpr),  18, .01, 1.0, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 25, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  26, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  10, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t8, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t9 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t9, file = paste0(path_modRDS, "mod_t9.rds"))

# - MCMC: thin 75% + ... {10:13} ----------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 4 # i.e. keep 1 out of every 4...

# - tune 10 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.13, n_betas$abpr),  15, .01, 1.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 26, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  26, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 26, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  11, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t9, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t10 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t10, file = paste0(path_modRDS, "mod_t10.rds"))

# - tune 11 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.13, n_betas$abpr),  14, .01, 1.4, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 27, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  26, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 26, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  13, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t10, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t11 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t11, file = paste0(path_modRDS, "mod_t11.rds"))

# - tune 12 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  14, .01, 1.6, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 26, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  26, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 26, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  13, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t11, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t12 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t12, file = paste0(path_modRDS, "mod_t12.rds"))

# - tune 13 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  15, .01, 1.8, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 25, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  26, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 26, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  14, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t12, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t13 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t13, file = paste0(path_modRDS, "mod_t13.rds"))

# - MCMC: thin 50% + ... {14:17} ----------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 2 # i.e. keep 1 out of every two...

# - tune 14 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  16, .01, 2.0, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 23, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  23, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 26, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  16, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t13, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t14 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t14, file = paste0(path_modRDS, "mod_t14.rds"))

# - tune 15 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  17, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 21, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  21, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  18, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t14, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t15 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t15, file = paste0(path_modRDS, "mod_t15.rds"))

# - tune 16 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  17, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 18, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  18, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  17, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t15, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t16 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t16, file = paste0(path_modRDS, "mod_t16.rds"))

# - tune 17 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  18, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 16, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  16, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  18, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t16, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t17 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t17, file = paste0(path_modRDS, "mod_t17.rds"))

# - MCMC: no thin  + ... {18:21} ----------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 1

# - tune 18 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  20, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 14, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  14, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  20, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t17, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t18 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t18, file = paste0(path_modRDS, "mod_t18.rds"))

# - tune 19 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  22, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 10, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.45, n_betas$qudo),  10, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  22, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t18, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t19 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t19, file = paste0(path_modRDS, "mod_t19.rds"))

# - tune 20 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  22, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 5, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  8, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  24, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t19, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t20 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t20, file = paste0(path_modRDS, "mod_t20.rds"))

# - tune 21 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  24, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 3, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  8, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  26, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t20, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t21 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t21, file = paste0(path_modRDS, "mod_t21.rds"))

# - MCMC: thin 50% + ... {22:25} ----------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 2 # i.e. keep 1 out of every two...

# - tune 22 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  24, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 3, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  8, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  27, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t21, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t22 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t22, file = paste0(path_modRDS, "mod_t22.rds"))

# - tune 23 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  24, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 3, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  9, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  27, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t22, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t23 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t23, file = paste0(path_modRDS, "mod_t23.rds"))

# - tune 24 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 3, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  10, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  26, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t23, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t24 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t24, file = paste0(path_modRDS, "mod_t24.rds"))

# - tune 25 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 3, .01, 1.3, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  11, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  25, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t24, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t25 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t25, file = paste0(path_modRDS, "mod_t25.rds"))

# - MCMC: thin 75% + ... {26} -------------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 4 # i.e. keep 1 out of every 4...

# - tune 26 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  24, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 3,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  12, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  23, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t25, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t26 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
saveRDS(mod_t26, file = paste0(path_modRDS, "mod_t26.rds"))

# - MCMC: thin 99% + ... {27} -------------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 100 # i.e. keep 1 out of every 4...

# - tune 27 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  22, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 5,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  13, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  24, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t26, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t27 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t27, file = paste0(path_modRDS, "mod_t27.rds"))

# - MCMC: thin 50% + ... {28:30} ----------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 2 # i.e. keep 1 out of every two...

# - tune 28 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  22, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  14, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .004, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  24, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t27, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t28 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
saveRDS(mod_t28, file = paste0(path_modRDS, "mod_t28.rds"))

# - tune 29 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  23, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  15, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  25, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t28, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t29 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
saveRDS(mod_t29, file = paste0(path_modRDS, "mod_t29.rds"))

# - tune 30 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  17, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  27, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t29, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t30 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
saveRDS(mod_t30, file = paste0(path_modRDS, "mod_t30.rds"))

# - MCMC: thin 99% + ... {31:34} ----------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 100 # i.e. keep 1 out of every 4...

# - tune 31 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  20, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  28, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t30, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
mod_t31 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)

# save rds
# saveRDS(mod_t31, file = paste0(path_modRDS, "mod_t31.rds"))

# - tune 32 (no change needed from tune t31!) ---------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .01, 1.4, .002, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  20, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  28, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t31, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t32 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t32, file = paste0(path_modRDS, "mod_t32.rds"))

# - tune 33 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .01, 1.4, .003, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  20, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  28, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t32, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t33 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t33, file = paste0(path_modRDS, "mod_t33.rds"))

# - tune 34 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01,  2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .015, 1.5, .004, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  20, .01,  1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 27, .01,  1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  28, .01,  0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t33, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t34 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t34, file = paste0(path_modRDS, "mod_t34.rds"))

# - MCMC: thin 50% + ... {35} -------------------------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 20
nThin   <- 2 # i.e. keep 1 out of every two...

# - tune 35 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  25, .01, 2.2, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .03, 1.7, .006, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  20, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  26, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t34, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t35 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t35, file = paste0(path_modRDS, "mod_t35.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {36:39} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two...

# - tune 36 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  26, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .03, 1.9, .008, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  21, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  24, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t35, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t36 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t36, file = paste0(path_modRDS, "mod_t36.rds"))

# - tune 37 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  26, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .03, 1.9, .008, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  23, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  20, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t36, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t37 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t37, file = paste0(path_modRDS, "mod_t37.rds"))

# - tune 38 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  26, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .03, 1.9, .008, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  23, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  10, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t37, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t38 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t38, file = paste0(path_modRDS, "mod_t38.rds"))

# - tune 39 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  26, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .03, 1.9, .008, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  25, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  5, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t38, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t39 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t39, file = paste0(path_modRDS, "mod_t39.rds"))

# - MCMC: thin 90% + 1k samps + 50:1 z's {40} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 10 # i.e. keep 1 out of every 10

# - tune 40 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  26, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.05, n_betas$psmem), 4,  .03, 1.9, .008, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  25, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4, .01, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t39, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t40 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t40, file = paste0(path_modRDS, "mod_t40.rds"))


# - MCMC: thin 99% + 1k samps + 50:1 z's {41} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 41 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  27, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .03, 2.0, .008, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.44, n_betas$qudo),  25, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t40, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t41 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t41, file = paste0(path_modRDS, "mod_t41.rds"))

# - MCMC: no thin  + 1k samps + 50:1 z's {42:44} ------------------------- ####m
# why?
# - dropping things down here just to tune some parameters before heading
#   back up to the 99+1k+50:1...
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 1

# - tune 42 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  27, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .04, 2.1, .009, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.43, n_betas$qudo),  21, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t41, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t42 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t42, file = paste0(path_modRDS, "mod_t42.rds"))

# - tune 43 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  10, .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .04, 2.4, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  10, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 10, .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t42, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t43 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t43, file = paste0(path_modRDS, "mod_t43.rds"))

# - tune 44 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .04, 2.6, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  10, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 5,  .01, 1.0, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t43, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t44 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t44, file = paste0(path_modRDS, "mod_t44.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {45} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 45 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .04, 2.6, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  15, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 4,  .01, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t44, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t45 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t45, file = paste0(path_modRDS, "mod_t45.rds"))

# - MCMC: thin 90% + 1k samps + 50:1 z's {46} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 10 # i.e. keep 1 out of every 10

# - tune 46 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .04, 2.6, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  17, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 4,  .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  6,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t45, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t46 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t46, file = paste0(path_modRDS, "mod_t46.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {47} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 47 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .04, 2.6, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  19, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 6,  .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  6,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t46, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t47 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t47, file = paste0(path_modRDS, "mod_t47.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {48:49} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 48 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .01, 2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05, 2.7, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  19, .01, 1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 5,  .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  8,  .02, 0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t47, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t48 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t48, file = paste0(path_modRDS, "mod_t48.rds"))

# - tune 49 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .02,  2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  2.7, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  19, .01,  1.2, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 4,  .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  20, .02,  0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t48, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t49 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t49, file = paste0(path_modRDS, "mod_t49.rds"))

# - MCMC: thin 99% + 2k samps + 50:1 z's {50} ---------------------------- ####m
nMCMC   <- 2000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 50 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .02,  2.1, .003, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  2.7, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  20, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 4,  .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  20, .02,  0.8, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
# - had to hard code this one... we want to grab the last iter of the previous
#   run... but there weren't 2k in the previous run... barrrrely caught this
#   in time hahaha!
dat_ls <- list(mod_t49, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t50 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t50, file = paste0(path_modRDS, "mod_t50.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {51:54} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 51 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  6,  .02,  2.1, .004, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  2.8, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  20, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 4,  .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  10, .015, 0.8, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t50, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 2000))

# run sampler
start_time <- Sys.time()
mod_t51 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t51, file = paste0(path_modRDS, "mod_t51.rds"))

# - tune 52 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  6,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  3.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  20, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 10, .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  5,  .015, 0.8, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t51, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t52 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t52, file = paste0(path_modRDS, "mod_t52.rds"))

# - tune 53 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  7,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  3.4, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  22, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 15, .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  4,  .015, 0.8, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t52, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t53 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t53, file = paste0(path_modRDS, "mod_t53.rds"))

# - tune 54 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  15,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  3.8, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  24, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .009, 1.0, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3,  .015, 0.8, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t53, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, nMCMC))

# run sampler
start_time <- Sys.time()
mod_t54 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
# saveRDS(mod_t54, file = paste0(path_modRDS, "mod_t54.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {55} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 55 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  15,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  3.8, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  24, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 25, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3,  .015, 0.8, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t54, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t55 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t55, file = paste0(path_modRDS, "mod_t55.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {56:58} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 56 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  5,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  4.4, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  24, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 5,  .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3,  .015, 0.8, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t55, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t56 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t56, file = paste0(path_modRDS, "mod_t56.rds"))

# - tune 57 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  4,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  20, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 5,  .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3,  .015, 1.1, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t56, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t57 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t57, file = paste0(path_modRDS, "mod_t57.rds"))

# - tune 58 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3,  .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4,  .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  10, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 5,  .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3,  .015, 1.1, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t57, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t58 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t58, file = paste0(path_modRDS, "mod_t58.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {59} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 59 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  8, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 5, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t58, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t59 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t59, file = paste0(path_modRDS, "mod_t59.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {60:61} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 60 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  6, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 4, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t59, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t60 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t60, file = paste0(path_modRDS, "mod_t60.rds"))

# - tune 61 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t60, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t61 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t61, file = paste0(path_modRDS, "mod_t61.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {62} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 62 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .01,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .003, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t61, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t62 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t62, file = paste0(path_modRDS, "mod_t62.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {63:64} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 63 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .02,  1.3, .005, .005, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t62, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t63 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t63, file = paste0(path_modRDS, "mod_t63.rds"))

# - tune 64 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .02,  1.3, .006, .006, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t63, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t64 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t64, file = paste0(path_modRDS, "mod_t64.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {65} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 65 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .02,  1.3, .006, .006, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .015, 1.1, .004, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t64, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t65 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t65, file = paste0(path_modRDS, "mod_t65.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {66:67} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 66 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .02,  1.3, .006, .006, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .02,  1.1, .005, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t65, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t66 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t66, file = paste0(path_modRDS, "mod_t66.rds"))

# - tune 67 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .02,  1.3, .006, .006, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .02,  1.1, .006, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t66, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t67 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t67, file = paste0(path_modRDS, "mod_t67.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {68} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 68 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .02,  2.1, .006, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .02,  1.3, .006, .006, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .02,  1.1, .006, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t67, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t68 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t68, file = paste0(path_modRDS, "mod_t68.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {69:70} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 69 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  2.1, .007, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .007, .007, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .007, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t68, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t69 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t69, file = paste0(path_modRDS, "mod_t69.rds"))

# - tune 70 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  2.1, .01, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01, .01, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .007, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t69, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t70 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t70, file = paste0(path_modRDS, "mod_t70.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {71} ---------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 71 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  2.1, .01, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01, .01, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .007, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t70, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t71 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t71, file = paste0(path_modRDS, "mod_t71.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {72:75} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 72 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  2.5, .01, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01, .01, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .008, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t71, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t72 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t72, file = paste0(path_modRDS, "mod_t72.rds"))

# - tune 73 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  3.5, .01, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01, .01, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .008, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t72, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t73 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t73, file = paste0(path_modRDS, "mod_t73.rds"))

# - tune 74 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  4.5, .01, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01, .01, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .008, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t73, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t74 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t74, file = paste0(path_modRDS, "mod_t74.rds"))

# - tune 75 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  5.5, .01, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01, .01, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .008, .06)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t74, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t75 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t75, file = paste0(path_modRDS, "mod_t75.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {76:77} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 76 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .03,  5.5, .01,  .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01,  .01,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .002, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .009, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t75, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t76 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t76, file = paste0(path_modRDS, "mod_t76.rds"))

# - tune 77 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  5.5, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .03,  1.3, .01,  .01,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.1, .01,  .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t76, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t77 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t77, file = paste0(path_modRDS, "mod_t77.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {78:79} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 78 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  5.8, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  1.5, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .03,  1.3, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t77, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t78 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t78, file = paste0(path_modRDS, "mod_t78.rds"))

# - tune 79 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  6.5, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.5, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .04,  1.3, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t78, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t79 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t79, file = paste0(path_modRDS, "mod_t79.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {80:83} ------------------------- ####m
# Note -- checked after t80, adjusted and then ran 3 times over weekend (81:83)
# - to avoid loss of data saving along the way as 81, 82, 83, each 1k long
#   instead of having 1 massive 3k run...
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 80 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  7.0, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.5, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .04,  1.3, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t79, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t80 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t80, file = paste0(path_modRDS, "mod_t80.rds"))

# - tune 81 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  7.5, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.5, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .04,  1.3, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t80, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# - clean-up space
rm(mod_t80)

# run sampler
start_time <- Sys.time()
mod_t81 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t81, file = paste0(path_modRDS, "mod_t81.rds"))

# - tune 82 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  7.5, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.5, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .04,  1.3, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t81, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# - clean-up space
rm(mod_t81)

# run sampler
start_time <- Sys.time()
mod_t82 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t82, file = paste0(path_modRDS, "mod_t82.rds"))

# - tune 83 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  7.5, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.5, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .002, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .04,  1.3, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t82, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# - clean-up space
rm(mod_t82)

# run sampler
start_time <- Sys.time()
mod_t83 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t83, file = paste0(path_modRDS, "mod_t83.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {84:87} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 84 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  7.5, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.4, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .05,  1.4, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t83, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t84 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t84, file = paste0(path_modRDS, "mod_t84.rds"))

# - tune 85 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .04,  7.0, .015, .10)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.4, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .003, .08)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .05,  2.0, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t84, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t85 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t85, file = paste0(path_modRDS, "mod_t85.rds"))

# - tune 86 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05,  7.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.3, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .003, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .06,  2.0, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t85, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t86 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t86, file = paste0(path_modRDS, "mod_t86.rds"))

# - tune 87 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05,  7.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.0, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .004, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .06,  2.0, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t86, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t87 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t87, file = paste0(path_modRDS, "mod_t87.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {88:90} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 88 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05,  7.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.0, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .004, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .06,  2.0, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t87, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t88 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t88, file = paste0(path_modRDS, "mod_t88.rds"))

# - tune 89 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05,  7.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.0, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .004, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .06,  2.0, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t88, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t89 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t89, file = paste0(path_modRDS, "mod_t89.rds"))

# - tune 90 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05,  7.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.0, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .003, .004, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .06,  2.0, .015, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t89, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t90 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t90, file = paste0(path_modRDS, "mod_t90.rds"))

# - MCMC: thin 50% + 1k samps + 50:1 z's {91:93} ------------------------- ####m
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 2 # i.e. keep 1 out of every two

# - tune 91 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05,  8.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05,  5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  3, .04,  2.0, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .009, 1.1, .004, .005, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .06,  2.0, .016, .07)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t90, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t91 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t91, file = paste0(path_modRDS, "mod_t91.rds"))

# - tune 92 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 2.0, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .004, .005, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t91, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t92 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t92, file = paste0(path_modRDS, "mod_t92.rds"))

# - tune 93 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .004, .005, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t92, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t93 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t93, file = paste0(path_modRDS, "mod_t93.rds"))

# - MCMC: thin 99% + 1k samps + 50:1 z's {94:113} ------------------------ ####
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 94 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .004, .005, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t93, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t94 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t94, file = paste0(path_modRDS, "mod_t94.rds"))

# - tune 95 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .005, .006, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t94, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t95 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t95, file = paste0(path_modRDS, "mod_t95.rds"))

# - tune 96 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .005, .006, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t95, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t96 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t96, file = paste0(path_modRDS, "mod_t96.rds"))

# - tune 97 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .006, .007, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t96, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t97 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t97, file = paste0(path_modRDS, "mod_t97.rds"))

# - tune 98 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .006, .007, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t97, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t98 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t98, file = paste0(path_modRDS, "mod_t98.rds"))

# - tune 99 -------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .006, .007, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t98, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t99 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t99, file = paste0(path_modRDS, "mod_t99.rds"))

# - tune 100 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.5, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.5, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015,  .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .007, .008, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .018, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t99, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t100 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t100, file = paste0(path_modRDS, "mod_t100.rds"))

# - tune 101 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.5, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.2, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .02, 1.1, .008, .01,  .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t100, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t101 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t101, file = paste0(path_modRDS, "mod_t101.rds"))

# - tune 102 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.2, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .010, .012, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t101, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t102 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t102, file = paste0(path_modRDS, "mod_t102.rds"))

# - tune 103 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.2, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .010, .012, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t102, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t103 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t103, file = paste0(path_modRDS, "mod_t103.rds"))

# - tune 104 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.2, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .010, .012, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t103, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t104 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t104, file = paste0(path_modRDS, "mod_t104.rds"))

# - tune 105 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 9.2, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .010, .012, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t104, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t105 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t105, file = paste0(path_modRDS, "mod_t105.rds"))

# - tune 106 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.7, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .012, .014, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t105, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t106 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t106, file = paste0(path_modRDS, "mod_t106.rds"))


# - tune 107 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.7, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .012, .014, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t106, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t107 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t107, file = paste0(path_modRDS, "mod_t107.rds"))


# - tune 108 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.7, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.1, .012, .014, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t107, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t108 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t108, file = paste0(path_modRDS, "mod_t108.rds"))

# - tune 109 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 7.5, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.6, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .03, 1.4, .012, .014, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t108, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t109 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t109, file = paste0(path_modRDS, "mod_t109.rds"))

# - tune 110 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.7, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.1, .013, .015, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t109, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t110 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t110, file = paste0(path_modRDS, "mod_t110.rds"))

# - tune 111 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.7, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.1, .013, .015, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t110, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t111 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t111, file = paste0(path_modRDS, "mod_t111.rds"))

# - tune 112 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.0, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.40, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.1, .013, .015, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t111, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t112 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t112, file = paste0(path_modRDS, "mod_t112.rds"))

# - tune 113 ------------------------------------------------------------- ####mm
# adjust tuning
dat_ls[["abpr"]]$tune_start  <- c(rep(.12, n_betas$abpr),  3, .05, 8.7, .015, .12)
dat_ls[["psmem"]]$tune_start <- c(rep(.06, n_betas$psmem), 4, .05, 5.0, .015, .05)
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.8, .015, .015, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.1, .013, .015, .1)
dat_ls[["quke"]]$tune_start  <- c(rep(.15, n_betas$quke),  3, .07, 2.0, .019, .08)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t112, dat_ls) %>% pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t113 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, n_cores = 5)
Sys.time() - start_time

# save rds
saveRDS(mod_t113, file = paste0(path_modRDS, "mod_t113.rds"))

# == Keep running {qudo, quga4}; done with {abpr, psmem, quke} =========== ####
run_spp <- c("qudo", "quga4")

# - MCMC: thin 99% + 1k samps + 50:1 z's {114:117} ----------------------- ####
nMCMC   <- 1000
nMCMC_z <- 50
nThin   <- 100 # i.e. keep 1 out of every 100

# - tune 114 ------------------------------------------------------------- ####
# adjust tuning
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.9, .016, .016, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.2, .014, .016, .1)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t113[run_spp], dat_ls[run_spp]) %>% 
  pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t114 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, 
                                   n_cores = length(dat_ls))
Sys.time() - start_time

# save rds
saveRDS(mod_t114, file = paste0(path_modRDS, "mod_t114.rds"))

# - tune 115 ------------------------------------------------------------- ####
# adjust tuning
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.9, .016, .016, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.2, .014, .016, .1)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t114[run_spp], dat_ls[run_spp]) %>% 
  pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t115 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, 
                                   n_cores = length(dat_ls))
Sys.time() - start_time

# save rds
saveRDS(mod_t115, file = paste0(path_modRDS, "mod_t115.rds"))

# - tune 116 ------------------------------------------------------------- ####
# adjust tuning
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.9, .016, .016, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.2, .014, .016, .1)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t115[run_spp], dat_ls[run_spp]) %>% 
  pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t116 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, 
                                   n_cores = length(dat_ls))
Sys.time() - start_time

# save rds
saveRDS(mod_t116, file = paste0(path_modRDS, "mod_t116.rds"))

# - tune 117 ------------------------------------------------------------- ####
# adjust tuning
dat_ls[["qudo"]]$tune_start  <- c(rep(.42, n_betas$qudo),  4, .04, 1.9, .016, .016, .12)
dat_ls[["quga4"]]$tune_start <- c(rep(.35, n_betas$quga4), 3, .04, 1.2, .014, .016, .1)
dat_ls <- dat_ls %>% map(function(x) {
  x$tune_start <- x$tune_start %>% setNames(x$parms_names)
  return(x)
})

# grab new starts
dat_ls <- list(mod_t116[run_spp], dat_ls[run_spp]) %>% 
  pmap(~fn.parms_restart(.x, .y, 1000))

# run sampler
start_time <- Sys.time()
mod_t117 <- fn_par_wrapper.CAR_fit(dat_ls, nMCMC, nMCMC_z, nThin, 
                                   n_cores = length(dat_ls))
Sys.time() - start_time

# save rds
saveRDS(mod_t117, file = paste0(path_modRDS, "mod_t117.rds"))

# == Peak at: current results ============================================ ####
m_temp <- mod_t117   # so we don't have to change the code below...

# - general summary (theta) ---------------------------------------------- ####
# trace (theta) & summaries
m_temp %>% map(~fn.mod_sum(.x)) # tune, accept rate, theta estimates
list(m_temp, dat_ls) %>% pmap(function(x, y) {
  fn.plot_trace_theta(x, NA, plot_actual = F, 
                      m_name = y$sp.x)
})  # theta trace
list(m_temp, dat_ls) %>% pmap(function(x, y) {
  fn.plot_density_theta(x, NA, plot_actual = F, m_name = y$sp.x)
}) # theta density

# mean acceptance rates by kept iter
m_temp %>% imap(~fn.plot_ar_theta(.x, .y))

# - are z's and v's mixing? ---------------------------------------------- ####
# identify z's and v's of interest (this only needs be done once...)
v_ids <- dat_ls %>% map(~fn.trace_v_ids(.x$moddat))
z_ids <- dat_ls %>% map(~fn.trace_z_ids(.x$moddat))

# trace (a few v_i's) & v summaries (all v's)
# list(m_temp, names(m_temp), v_ids) %>% pmap(~fn.plot_trace_zv(..1, ..2, ..3, "v"))
m_temp %>% map(~fn.est_zv(.x, "v") %>% select(-ES) %>% apply(2, summary))
# m_temp %>% map(~fn.est_zv(.x, "v", by_iter_T = TRUE))

# trace (a few z_i's) & v summaries (all z's)
# list(m_temp, names(m_temp), z_ids) %>% pmap(~fn.plot_trace_zv(..1, ..2, ..3, "z"))
m_temp %>% map(~fn.est_zv(.x, "z") %>% select(-fiahex_id) %>% apply(2, summary))
# m_temp %>% map(~fn.est_zv(.x, "z", by_iter_T = TRUE))

# # distribution of mean z or v
# m_temp %>% imap(~fn.plot_density_zv_mean(.x, .y, "v"))
# m_temp %>% imap(~fn.plot_density_zv_mean(.x, .y, "z"))

# # - z & v details ------------------------------------------------------ ####
# # extract deets: z
# list(dat_ls, z_ids) %>% pmap(function(dat.x, z.x) {
#   vars.x <- names(dat.x$moddat)
#   vars.x <- vars.x[grepl("bio", vars.x) |grepl("gdd5", vars.x) ]
#   vars.x <- vars.x[!grepl("_", vars.x)]
#   
#   dat.x$moddat %>%
#     select(fiahex_id, p_a, y_index, fiahex_obsv, all_of(vars.x)) %>%
#     filter(fiahex_id %in% z.x) %>%
#     group_by(fiahex_id) %>%
#     mutate(prop_p = sum(p_a, na.rm = TRUE)/n(),
#            prop_a = sum(!p_a, na.rm = TRUE)/n(), 
#            n_plots = n(),
#            prop_ob = sum(y_index)/n()) %>%
#     group_by(fiahex_id, prop_p, prop_a, prop_ob, n_plots) %>%
#     summarise_at(vars.x, list(~min(.)))
#   
# })
# 
# # extract deets: v 
# list(dat_ls, v_ids) %>% pmap(function(dat.x, v.x) {
#   vars.x <- names(dat.x$moddat)
#   vars.x <- vars.x[grepl("bio", vars.x) |grepl("gdd5", vars.x) ]
#   vars.x <- vars.x[!grepl("_", vars.x)]
#   
#   dat.x$moddat %>%
#     select(ES, p_a, y_index, fiahex_obsv, all_of(vars.x), fiahex_id) %>%
#     filter(ES %in% v.x) %>%
#     group_by(ES) %>%
#     mutate(prop_p = (sum(p_a, na.rm = TRUE)/n()) %>% round(2),
#            prop_a = (sum(!p_a, na.rm = TRUE)/n()) %>% round(2), 
#            n_plots = n(),
#            n_hex = n_distinct(fiahex_id),
#            prop_ob = (sum(y_index)/n()) %>% round(2)) %>%
#     group_by(ES, prop_p, prop_a, prop_ob, n_plots, n_hex) %>%
#     summarise_at(vars.x, list(~min(.) %>% round(2)))
#   
# })
# 
# # - map a z sample ----------------------------------------------------- ####
# load spatial plotting data if needed...
# - note that the commented files below are produced in script 8... moved 
#   to that script later on...
# library(sf)
# hex_polygons_public.sf <- readRDS(paste0(path_output,
#                                  "ch2-8-hex_polygons_public_sf.rds"))
# study_area.sf  <-  readRDS(paste0(path_output,
#                                   "ch2-8-study_area_sf.rds"))
# 
# # side by side - z sample at iter 1000 & expected z maps
# list(m_temp, hex_polygons_public.sf, names(m_temp)) %>% 
#   pmap(~cowplot::plot_grid(
#     # map based on iteration 1000
#     fn.map_z_sample_sf(..1, ..2, ..3, 
#                        study_area.sf,
#                        MCMC_iter = nMCMC,
#                        expect_T = FALSE,
#                        color_mu_T = TRUE),
#     # map based on expected z
#     fn.map_z_sample_sf(..1, ..2, ..3, 
#                        study_area.sf,
#                        MCMC_iter = nMCMC,
#                        expect_T = TRUE,
#                        color_mu_T = TRUE),
#     nrow = 1))
# == Peak at: long results =============================================== ####
# - figure out how many samples to keep and their ids -------------------- ####
# runs to sample from & their lengths & how many runs will we merge
run_ids   <- 108:117
run_nMCMC <- 1000
run_n     <- length(run_ids)

# number of samples we want at the end
# - i.e. how many samples after merging & thinning these runs do we want?
n_samps <- 2000

# iteration ids to keep from each run to get n_samps
n_subthin <- (run_n*run_nMCMC)/n_samps # this is like our new nThin.. keep one out of this many
keep_iter_id <- n_subthin*c(1:(n_samps/run_n))

# clean-up space
rm(run_nMCMC, run_n, n_subthin)

# - sequentially subsample and merge model results ----------------------- ####
# read in and subset model-runs to keep_iter_id
m_temp <- lapply(run_ids, function(i) {
  readRDS(paste0(path_modRDS, "mod_t", i, ".rds")) %>%
    map(function(sp.x) {
      sp.x[c("MCMC_theta", "MCMC_z", "MCMC_v")] %>%
        modify_at("MCMC_theta", ~.x[keep_iter_id, ]) %>%
        modify_at("MCMC_z", ~.x[keep_iter_id, ]) %>%
        modify_at("MCMC_v", ~.x[keep_iter_id, ])
    })
}) %>% 
  setNames(paste0("t", run_ids))

# merge mod runs by spp, add new 'n_sim' for order to MCMC_theta element
m_temp <- m_temp %>% 
  transpose %>%
  map(function(sp.x) {
    sp.x %>% 
      transpose %>%
      modify_at("MCMC_theta", function(mod.x) {
        mod.x %>%
          bind_rows(.id = "n_iter") %>%
          mutate(n_iter = gsub("t", "", n_iter) %>% as.numeric) %>%
          arrange(n_iter, n_sim) %>%
          mutate(n_sim = 1:n()) %>%
          select(-n_iter)
      }) %>%
      modify_at(c("MCMC_z", "MCMC_v"), ~.x %>% bind_rows)
  })

# - general summary (theta) ---------------------------------------------- ####
# trace (theta) & summaries
m_temp %>% map(~fn.est_theta(.x)) # tune, accept rate, theta estimates
list(m_temp, dat_ls) %>% pmap(function(x, y) {
  fn.plot_trace_theta(x, NA, plot_actual = F, 
                      m_name = y$sp.x)
})  # theta trace
list(m_temp, dat_ls) %>% pmap(function(x, y) {
  fn.plot_density_theta(x, NA, plot_actual = F, m_name = y$sp.x)
}) # theta density

# - are z's and v's mixing? ---------------------------------------------- ####
# identify z's and v's of interest (this only needs be done once...)
v_ids <- dat_ls %>% map(~fn.trace_v_ids(.x$moddat))
z_ids <- dat_ls %>% map(~fn.trace_z_ids(.x$moddat))

# trace (a few v_i's) & v summaries (all v's)
list(m_temp, names(m_temp), v_ids) %>% pmap(~fn.plot_trace_zv(..1, ..2, ..3, "v"))
list(m_temp, names(m_temp), v_ids) %>% pmap(~fn.plot_trace_zv(..1, ..2, ..3, "v", density_T = TRUE))
m_temp %>% map(~fn.est_zv(.x, "v") %>% select(-ES) %>% apply(2, summary))
m_temp %>% map(~fn.est_zv(.x, "v", by_iter_T = TRUE))
               
# trace (a few z_i's) & v summaries (all z's)
list(m_temp, names(m_temp), z_ids) %>% pmap(~fn.plot_trace_zv(..1, ..2, ..3, "z"))
list(m_temp, names(m_temp), z_ids) %>% pmap(~fn.plot_trace_zv(..1, ..2, ..3, "z", density_T = TRUE))
m_temp %>% map(~fn.est_zv(.x, "z") %>% select(-fiahex_id) %>% apply(2, summary))
m_temp %>% map(~fn.est_zv(.x, "z", by_iter_T = TRUE))

# - map a z sample ------------------------------------------------------- ####
# load spatial plotting data if needed...
library(sf)
hex_polygons_public.sf <- readRDS(paste0(path_output,
                                 "ch2-8-hex_polygons_public_sf.rds"))
study_area.sf  <-  readRDS(paste0(path_output,
                                  "ch2-8-study_area_sf.rds"))

# side by side - z sample at iter 1000 & expected z maps
list(m_temp, hex_polygons_public.sf, names(m_temp)) %>% 
  pmap(~cowplot::plot_grid(
    # map based on iteration 1000
    fn.map_z_sample_sf(..1, ..2, ..3, 
                       study_area.sf,
                       MCMC_iter = n_samps,
                       expect_T = FALSE,
                       color_mu_T = TRUE),
    # map based on expected z
    fn.map_z_sample_sf(..1, ..2, ..3, 
                       study_area.sf,
                       MCMC_iter = n_samps,
                       expect_T = TRUE,
                       color_mu_T = TRUE),
    nrow = 1))

# - z & v details -------------------------------------------------------- ####
# extract deets: z
list(dat_ls, z_ids) %>% pmap(function(dat.x, z.x) {
  vars.x <- names(dat.x$moddat)
  vars.x <- vars.x[grepl("bio", vars.x) |grepl("gdd5", vars.x) ]
  vars.x <- vars.x[!grepl("_", vars.x)]

  dat.x$moddat %>%
    select(fiahex_id, p_a, y_index, fiahex_obsv, all_of(vars.x)) %>%
    filter(fiahex_id %in% z.x) %>%
    group_by(fiahex_id) %>%
    mutate(prop_p = sum(p_a, na.rm = TRUE)/n(),
           prop_a = sum(!p_a, na.rm = TRUE)/n(),
           n_plots = n(),
           prop_ob = sum(y_index)/n()) %>%
    group_by(fiahex_id, prop_p, prop_a, prop_ob, n_plots) %>%
    summarise_at(vars.x, list(~min(.)))

})

# extract deets: v
list(dat_ls, v_ids) %>% pmap(function(dat.x, v.x) {
  vars.x <- names(dat.x$moddat)
  vars.x <- vars.x[grepl("bio", vars.x) |grepl("gdd5", vars.x) ]
  vars.x <- vars.x[!grepl("_", vars.x)]

  dat.x$moddat %>%
    select(ES, p_a, y_index, fiahex_obsv, all_of(vars.x), fiahex_id) %>%
    filter(ES %in% v.x) %>%
    group_by(ES) %>%
    mutate(prop_p = (sum(p_a, na.rm = TRUE)/n()) %>% round(2),
           prop_a = (sum(!p_a, na.rm = TRUE)/n()) %>% round(2),
           n_plots = n(),
           n_hex = n_distinct(fiahex_id),
           prop_ob = (sum(y_index)/n()) %>% round(2)) %>%
    group_by(ES, prop_p, prop_a, prop_ob, n_plots, n_hex) %>%
    summarise_at(vars.x, list(~min(.) %>% round(2)))

})

# == Assemble final MCMC samples ========================================= ####
# final MCMC samples & burn-in come from these model-runs for these species:
# - {abpr, psmem, quke}: 104:113 (mod_spp) & 102:103 (burn-in)
# - {qudo, quga4}:       108:117 (mod_spp) & 106:107 (burn-in)

mod_spp <- vector("list", length(spp)) %>% setNames(spp)
burnin_spp <- vector("list", length(spp)) %>% setNames(spp)

# - figure out how many samples to keep and their ids -------------------- ####
# runs to sample from & their lengths & how many runs will we merge
run_n     <- 10
run_nMCMC <- 1000

# number of samples we want at the end
# - i.e. how many samples after merging & thinning these runs do we want?
n_samps <- 2000

# iteration ids to keep from each run to get n_samps
n_subthin <- (run_n*run_nMCMC)/n_samps # this is like our new nThin.. keep one out of this many
keep_iter_id <- n_subthin*c(1:(n_samps/run_n))

# clean-up space
rm(run_nMCMC, run_n, n_subthin)

# - load final runs: {abpr, psmem, quke} --------------------------------- ####
spp.x <- c("abpr", "psmem", "quke")
run_ids.x <- 104:113
burnin_ids.x <- 102:103

# read in and subset model-runs to keep_iter_id, merge mod runs by spp, add 
# new 'n_sim' for order to MCMC_theta element
mod_spp[spp.x] <- lapply(run_ids.x, function(i) {
  (readRDS(paste0(path_modRDS, "mod_t", i, ".rds")))[spp.x] %>%
    map(function(sp.x) {
      sp.x[c("MCMC_theta", "MCMC_z", "MCMC_v")] %>%
        modify_at("MCMC_theta", ~.x[keep_iter_id, ]) %>%
        modify_at("MCMC_z", ~.x[keep_iter_id, ]) %>%
        modify_at("MCMC_v", ~.x[keep_iter_id, ])
    })
}) %>% 
  setNames(paste0("t", run_ids.x)) %>% 
  transpose %>%
  map(function(sp.x) {
    sp.x %>% 
      transpose %>%
      modify_at("MCMC_theta", function(mod.x) {
        mod.x %>%
          bind_rows(.id = "n_iter") %>%
          mutate(n_iter = gsub("t", "", n_iter) %>% as.numeric) %>%
          arrange(n_iter, n_sim) %>%
          mutate(n_sim = 1:n()) %>%
          select(-n_iter)
      }) %>%
      modify_at(c("MCMC_z", "MCMC_v"), ~.x %>% bind_rows)
  })

# make a burn-in runs object too
burnin_spp[spp.x] <- lapply(burnin_ids.x, function(i) {
  (readRDS(paste0(path_modRDS, "mod_t", i, ".rds")))[spp.x]
}) %>% 
  setNames(paste0("t", burnin_ids.x)) %>% 
  transpose %>%
  map(function(sp.x) {
    sp.x %>% 
      transpose %>%
      modify_at("MCMC_theta", function(mod.x) {
        mod.x %>%
          bind_rows(.id = "n_iter") %>%
          mutate(n_iter = gsub("t", "", n_iter) %>% as.numeric) %>%
          arrange(n_iter, n_sim) %>%
          mutate(n_sim = 1:n()) %>%
          select(-n_iter)
      }) %>%
      modify_at(c("MCMC_z", "MCMC_v"), ~.x %>% bind_rows)
  })

# - load final runs: {qudo, quga4} --------------------------------------- ####
spp.x <- c("qudo", "quga4")
run_ids.x <- 108:117
burnin_ids.x <- 106:107

# read in and subset model-runs to keep_iter_id, merge mod runs by spp, add 
# new 'n_sim' for order to MCMC_theta element
mod_spp[spp.x] <- lapply(run_ids.x, function(i) {
  (readRDS(paste0(path_modRDS, "mod_t", i, ".rds")))[spp.x] %>%
    map(function(sp.x) {
      sp.x[c("MCMC_theta", "MCMC_z", "MCMC_v")] %>%
        modify_at("MCMC_theta", ~.x[keep_iter_id, ]) %>%
        modify_at("MCMC_z", ~.x[keep_iter_id, ]) %>%
        modify_at("MCMC_v", ~.x[keep_iter_id, ])
    })
}) %>% 
  setNames(paste0("t", run_ids.x)) %>% 
  transpose %>%
  map(function(sp.x) {
    sp.x %>% 
      transpose %>%
      modify_at("MCMC_theta", function(mod.x) {
        mod.x %>%
          bind_rows(.id = "n_iter") %>%
          mutate(n_iter = gsub("t", "", n_iter) %>% as.numeric) %>%
          arrange(n_iter, n_sim) %>%
          mutate(n_sim = 1:n()) %>%
          select(-n_iter)
      }) %>%
      modify_at(c("MCMC_z", "MCMC_v"), ~.x %>% bind_rows)
  })


# make a burn-in runs object too
burnin_spp[spp.x] <- lapply(burnin_ids.x, function(i) {
  (readRDS(paste0(path_modRDS, "mod_t", i, ".rds")))[spp.x]
}) %>% 
  setNames(paste0("t", burnin_ids.x)) %>% 
  transpose %>%
  map(function(sp.x) {
    sp.x %>% 
      transpose %>%
      modify_at("MCMC_theta", function(mod.x) {
        mod.x %>%
          bind_rows(.id = "n_iter") %>%
          mutate(n_iter = gsub("t", "", n_iter) %>% as.numeric) %>%
          arrange(n_iter, n_sim) %>%
          mutate(n_sim = 1:n()) %>%
          select(-n_iter)
      }) %>%
      modify_at(c("MCMC_z", "MCMC_v"), ~.x %>% bind_rows)
  })

# - Save mod_spp (final MCMC samples) ------------------------------------ ####
saveRDS(mod_spp, 
        file = paste0(path_output, "ch2-6-final-MCMC-samples.rds"))

# - Save burnin_spp (burn-in for final MCMC samples) --------------------- ####
saveRDS(burnin_spp, 
        file = paste0(path_output, "ch2-6-final-MCMC-burnin.rds"))

# == FINAL -- Total number of iterations ================================= ####
iter_counter <- matrix(
  c(# when I was doing 20:1 on the z's
    rep(c(1000, 20, 1),  (6+4)),
    rep(c(1000, 20, 2), (3+4+4+3+1)),
    rep(c(1000, 20, 75), (4+1)),
    rep(c(1000, 20, 100), (5)),
    # when I was doing 50:1 on the z's
    rep(c(1000, 50, 1),  (3)),
    rep(c(1000, 50, 2), (4+1+2)),
    rep(c(1000, 50, 10), (2)),
    rep(c(1000, 50, 100), (2)),
    # starting at iter 50...
    rep(c(2000, 50, 100), (1)),
    rep(c(1000, 50, 2), (4+3+2+2+2+2)),
    rep(c(1000, 50, 100), (1+1+1+1+1+1)),
    # starting at iter 72..
    rep(c(1000, 50, 2), (4+2+4+3)),
    rep(c(1000, 50, 100), (2+4+3+20))),
  # end abpr, psmem, quke
  ncol = 3, byrow = TRUE
) %>% 
  as.data.frame
names(iter_counter) <- c("nMCMC", "nMCMC_z", "thin_mult")
iter_counter$run_type <- "all_sp"

# add qudo & quga4 extra runs...
add_q2 <- as.data.frame(matrix(rep(c(1000, 50, 100), (4)), 
                               ncol = 3, byrow = TRUE))
names(add_q2) <- c("nMCMC", "nMCMC_z", "thin_mult")
add_q2$run_type <- "q2"

# put these together...
iter_counter <- rbind.data.frame(iter_counter, add_q2) %>%
  mutate(run_n = 1:n())
rm(add_q2)

# summarize number of runs...
iter_counter %>%
  ungroup() %>%
  # filter(run_n < 102) %>% # iters before convergence for abpr/psmem/quke
  # filter(run_n < 106) %>% # iters before convergence for qudo/quga4
  mutate(iter_non_z = nMCMC * thin_mult,
         iter_z = nMCMC * nMCMC_z * thin_mult) %>%
  group_by(run_type) %>%
  summarize(total_iternonz = sum(iter_non_z)/1000000,
            total_iterz = sum(iter_z)/1000000)
