###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 0-functions.R
## Updated - 01-28-2023
## Author  - Karin Kralicek (karin.kralicek@udsa.gov)
## 
## This R script is in support of 'Spatial-Bayesian models project shifts in 
## suitable habitat for Pacific Northwest tree species under climate change' by
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - functions used in this project (stored here to keep things tidy)
##
###############################################################################
library(magrittr)    # for set_rownames()
library(nabor)       # for `knn()` neighbor finding function
library(Matrix)      # for `sparseMatrix()`
library(dplyr)       # ... for E V E R Y T H I N G ...
library(purrr)       # working with lists
library(reshape2)    # melt() and dcast(); alt. update with tidyr's pivot fns
library(stringr)
library(lme4)
library(parallel)
library(sf)

# == Fit CAR models ====================================================== ####
# - fn.inv_logit --------------------------------------------------------- ####
fn.inv_logit <- function(x) exp(x) / (1 + exp(x))

# - fn.draw: draw candidates from proposal ------------------------------- ####
fn.draw <- list(
  beta    = function(b, v.sd) rnorm(1, b, v.sd),
  rho     = function(b, v.half_width, rho_values) {
    # v.half_width: sets the number of values to pull from the look-up table
    # (eg, v.half_width = 3 for a discrete uniform over 7 values)
    b_index <- which.min(abs(rho_values - b))
    if(b_index < v.half_width) {
      sample(rho_values[1:(v.half_width*2+1)], 1)
    } else if (b_index > (length(rho_values) - v.half_width)) {
      sample(rho_values[(length(rho_values) - v.half_width*2):length(rho_values)], 1)
    } else {
      sample(rho_values[(b_index - v.half_width):(b_index + v.half_width)], 1)
    }
  },
  sigma_z = function(b, v.half_width, ...) {
    if(b < v.half_width) {
      runif(1, 0, v.half_width*2)
    } else {
      runif(1, (b-v.half_width), (b+v.half_width))
    }
  },
  sigma_v = function(b, v.half_width, ...) {
    if(b < v.half_width) {
      runif(1, 0, v.half_width*2)
    } else {
      runif(1, (b-v.half_width), (b+v.half_width))
    }
  },
  z_ob    = function(b, v.sd) {
    # block proposal of the z_ob's
    rnorm(length(b), b, v.sd)
  },
  z_un    = function(b, v.sd) {
    # block proposal of the z_un's
    rnorm(length(b), b, v.sd)
  },
  v       = function(b, v.sd, ...) {
    # block proposal of the v's
    # - constrained like z's, but centered to agree w mean=0 assumption
    new_b <- rnorm(length(b), b, v.sd)
    new_b <- new_b - mean(new_b)
    return(new_b)
  })

# - fn.llike: calculate the (prop) log-likelihood ------------------------ ####
fn.llike <- list(
  beta    = function(parms, moddat, Q, U, W, z_index, y_index) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- parms[["z_ob"]]
    z[!as.logical(z_index)] <- parms[["z_un"]]
    
    # vector with correct z repeated (e.g. line up fiahex_id's with plot_id's)
    Uz <- as.vector(U %*% z)
    
    # vector with correct v repeated (e.g. line up ES's with plot_id's)
    Qv <- as.vector(Q %*% parms[["v"]])
    
    # betas
    p <- sum(grepl("beta_", names(parms)))
    XB <- parms[["beta_0"]] + as.matrix(moddat[,1:(p-1)]) %*% unlist(parms[2:p])
    
    # mu_hat
    mu <- fn.inv_logit(XB + Qv + Uz)
    
    # only observed y
    p_a <- moddat$p_a[as.logical(y_index)]
    mu  <- mu[as.logical(y_index)]
    
    # return the prop log-like
    sum(dbinom(p_a, 1, mu, log = TRUE))
  },
  rho     = function(parms, W, z_index, rho_lookup, rho_values, ...) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- parms[["z_ob"]]
    z[!as.logical(z_index)] <- parms[["z_un"]]
    
    b_index <- which.min(abs(rho_values - parms[["rho"]]))
    R <- rho_lookup[[b_index]]$R
    R.det_log <- rho_lookup[[b_index]]$R.det_log
    
    sig_z <- parms[["sigma_z"]]
    
    ((t(z) %*% R %*% z) / (-2 * sig_z^2) + R.det_log/2)[1]
  },
  sigma_z = function(parms, W, z_index, rho_lookup, rho_values, ...) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- parms[["z_ob"]]
    z[!as.logical(z_index)] <- parms[["z_un"]]
    
    b_index <- which.min(abs(rho_values - parms[["rho"]]))
    R <- rho_lookup[[b_index]]$R
    
    sig_z <- parms[["sigma_z"]]
    
    ((t(z) %*% R %*% z) / (-2 * sig_z^2) - length(z) * log(sig_z))[1]
  },
  sigma_v = function(parms, ...) {
    v <- parms[["v"]]
    sig_v <- parms[["sigma_v"]]
    
    ((t(v) %*% v) / (-2 * sig_v^2) - length(v) * log(sig_v))[1]
  },
  z_ob    = function(parms, moddat, Q, U, W, z_index, y_index, rho_lookup, rho_values) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- parms[["z_ob"]]
    z[!as.logical(z_index)] <- parms[["z_un"]]
    
    # vector with correct z repeated (e.g. line up fiahex_id's with plot_id's)
    Uz <- as.vector(U %*% z)
    
    # vector with correct v repeated (e.g. line up ES's with plot_id's)
    Qv <- as.vector(Q %*% parms[["v"]])
    
    # betas
    p <- sum(grepl("beta_", names(parms)))
    XB <- parms[["beta_0"]] + as.matrix(moddat[,1:(p-1)]) %*% unlist(parms[2:p])
    
    # mu_hat
    mu <- fn.inv_logit(XB + Qv + Uz)
    
    # only observed y
    p_a <- moddat$p_a[as.logical(y_index)]
    mu  <- mu[as.logical(y_index)]
    
    # return the prop log-like
    b_index <- which.min(abs(rho_values - parms[["rho"]]))
    R <- rho_lookup[[b_index]]$R
    sig_z <- parms[["sigma_z"]]
    
    sum(dbinom(p_a, 1, mu, log = TRUE)) + ((t(z) %*% R %*% z) / (-2 * sig_z^2))[1]
  },
  z_un    = function(parms, W, z_index, rho_lookup, rho_values, ...) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- parms[["z_ob"]]
    z[!as.logical(z_index)] <- parms[["z_un"]]
    
    b_index <- which.min(abs(rho_values - parms[["rho"]]))
    R <- rho_lookup[[b_index]]$R
    sig_z <- parms[["sigma_z"]]
    
    ((t(z) %*% R %*% z) / (-2 * sig_z^2))[1]
  },
  v       = function(parms, moddat, Q, U, W, z_index, y_index, ...) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- parms[["z_ob"]]
    z[!as.logical(z_index)] <- parms[["z_un"]]
    
    # vector with correct z repeated (e.g. line up fiahex_id's with plot_id's)
    Uz <- as.vector(U %*% z)
    
    # vector with correct v repeated (e.g. line up ES's with plot_id's)
    v  <- parms[["v"]]
    Qv <- as.vector(Q %*% v)
    
    # betas
    p <- sum(grepl("beta_", names(parms)))
    XB <- parms[["beta_0"]] + as.matrix(moddat[,1:(p-1)]) %*% unlist(parms[2:p])
    
    # mu_hat
    mu <- fn.inv_logit(XB + Qv + Uz)
    
    # only observed y
    p_a <- moddat$p_a[as.logical(y_index)]
    mu  <- mu[as.logical(y_index)]
    
    # return the prop log-like
    sig_v <- parms[["sigma_v"]]
    
    sum(dbinom(p_a, 1, mu, log = TRUE)) + ((t(v) %*% v) / (-2 * sig_v^2))[1]
  }
)

# - fn.CAR_fit: the MH sampler ------------------------------------------- ####
# note: moved run_time immediately around sampler (ie, inside this fn)
fn.CAR_fit <- function(moddat, W_ls, U, Q, rho_lookup,
                       parms_names, parms_start, tune_start,
                       nMCMC, nMCMC_z, nThin,
                       fix_rho_TF = FALSE) {
  ## About: ----
  ##  This is Metropolis-Hastings MCMC sampler to estimate model parameters 
  ##  and spatial random effects (CAR) for a Bernoulli response (0/1) that
  ##  may have some missing response data.
  ## 
  ## Requires: ----
  ##  Raw data info: ----
  ##  - moddat: data.frame of data required for modeling
  ##    - plot_id: unique identifier for each plot observation
  ##    - y_index: indicates whether a plot was observed (1) or unobserved (0)
  ##    - fiahex_id: unique identifier for each FIA hex-polygon; spatial 
  ##      random effects (e.g. z's) are associated at this level
  ##    - covariate data already manipulated into required variable forms
  ##      (e.g. if interaction between 'X1' and 'X2' required, a column exists
  ##      in moddat that equals X1*X2) 
  ##    - ES: ecological section membership (RE in model)
  ##    - p_a: response (potentially with missing data)
  ##  - W_ls: named list containing
  ##    - W: row standardized sparse matrix of neighbor relations
  ##    - z_index: indicates whether all plots in a hexagon were unobserved (0)
  ##      or if they were either all observed or mixed (1)
  ##  - U: design matrix relating plot_id's to fiahex_id's
  ##  - Q: design matrix relating plot_id's to ES's
  ##  - rho_lookup: pre-created list of look-up tables for rho/R/R.det_log
  ##  MCMC integers: ----
  ##  - nMCMC: number of iterations (excluding thins)
  ##  - nMCMC_z: # of iters for z_ob & z_un for every 1 iter of other prams
  ##  - nThin: number of iterations to thin 
  ##    (i.e., # before storing and moving on to next iteration)
  ##  Named parms vectors: ----
  ##  - parms_names: vector of parameter names
  ##  - parms_start: starting values for the parameters, (z & v as NAs)
  ##  - tune_start: tuning parameter values for drawing samples
  ##    - betas ~ N(b, *) where (*) are sd tuning parameters
  ##    - rho ~ discrete Unif(b, *) where * is the interval half_width tuning 
  ##      parameter centered on b
  ##    - sigma_z ~ Unif(b, *) where * is the interval half_width tuning 
  ##      parameter centered on b
  ##    - sigma_v ~ Unif(b, *) where * is the interval half_width tuning 
  ##      parameter centered on b
  ##    - z_ob ~ N(b, *) where * is the sd tuning parameter centered on b
  ##    - z_un ~ N(b, *) where * is the sd tuning parameter centered on b
  ##    - v ~ N(0, *) where * is the sd tuning parameter
  ##  Development testing: ----
  ##  - fix_rho_TF: T/F should rho be fixed (testing rho estimation options)
  ## 
  
  # == Set up ============================================================ ####
  n_parms <- length(parms_names)
  W <- W_ls$W
  
  # get rho options (for use in fn.draw[["rho"]] & fn.llike fn's)
  rho_values <- rho_lookup %>% map(~.x$rho) %>% unlist
  
  # pull out z_index and y_index
  z_index <- W_ls$z_index
  y_index <- moddat$y_index
  
  # create a candidate vector for iters & add z,v starting values if needed
  # - convert to a list, so that z- and v-blocks can be stored...
  # - rho note: sampling from look-up table, so make sure rho starting value
  #   is moved over to nearest similar look-up table value...
  # - if-statement note: b/c these are empty on initial start before tuning
  parms_try <- parms_start
  
  rho_index <- which.min(abs(rho_values - parms_try$rho)) 
  parms_try$rho <- rho_values[rho_index]
  
  if(is.na(sum(parms_try$z_ob))) {
    parms_try$z_ob <- rnorm(sum(z_index), 0, tune_start["z_ob"])
  }
  if(any(y_index == 0) & any(is.na(parms_try$z_un))) {
    parms_try$z_un <- rnorm(sum(!z_index), 0, tune_start["z_un"])
  }
  if(is.na(sum(parms_try$v))) {
    parms_try$v    <- rnorm(ncol(Q), 0, tune_start["v"])
  }
  
  # where to store things...
  # - accepted samples
  parms_current <- vector("list", nMCMC)
  parms_current <- lapply(parms_current, function(x) {
    x <- vector("list", length = n_parms)
    names(x) <- parms_names
    return(x)
  })
  parms_current[[1]] <- parms_try
  
  # - accept/reject rates
  parms_accept_rate <- lapply(vector("list", nMCMC), function(x) {
    lapply(vector("list", nThin), function(x) {
      ar_base <- vector("list", length(parms_names)) %>% setNames(parms_names)
      ar_base[["z_ob"]] <- vector("list", nMCMC_z)
      if("z_un" %in% parms_names) {
        ar_base[["z_un"]] <- vector("list", nMCMC_z)
      }
      
      return(ar_base)
      
    })
  })
  
  # == The sampler ======================================================= ####
  start_time <- Sys.time()
  
  for(i in 1:nMCMC) {
    # - Print every 10th iter -------------------------------------------- ####
    if(i %% 10==0) cat(paste0("MCMC iteration: ", i, "\n"))
    
    # - Thinning iters --------------------------------------------------- ####
    for (i_thin in 1:nThin) {
      # Basically, here we sit at index "i" nThin times before exiting the loop
      # and updating (e.g. before moving on to iter "i + 1")
      for (x in 1:length(parms_names)) {
        # - MH sampler --------------------------------------------------- ####
        # Because the z's are such large block proposals, we need to sample
        # longer for z_ob & z_un... we don't need to sample as much for the 
        # others, so break these parameters apart with an if-statement to save
        # time...
        
        if (grepl("beta", parms_names[x])) {
          # MH sampler for the betas
          # 1. Draw a candidate from the proposal ------------------------ ####
          # Draw candidates for a beta - function name hard coded
          parms_try[[x]] <- fn.draw[["beta"]](parms_current[[i]][[x]], 
                                              tune_start[x])
          
          # 2. Calc dif in the log-likelihood ---------------------------- ####
          # LLdif for a beta
          LLdif <- 
            (fn.llike[["beta"]](parms = parms_try, 
                                moddat = moddat,
                                Q = Q,
                                U = U,
                                W = W,
                                z_index = z_index,
                                y_index = y_index) -
               fn.llike[["beta"]](parms = parms_current[[i]], 
                                  moddat = moddat,
                                  Q = Q,
                                  U = U,
                                  W = W,
                                  z_index = z_index,
                                  y_index = y_index))
          
          # 3. Accept or reject update ----------------------------------- ####
          rU    <- log(runif(1))
          
          # Store sample value for this iteration
          if(rU < LLdif) {parms_current[[i]][[x]] <- parms_try[[x]]}
          parms_try[[x]] <- parms_current[[i]][[x]]
          
          # Store accept/reject decision
          parms_accept_rate[[i]][[i_thin]][[x]] <- rU < LLdif
          
        } else if (parms_names[x] == "rho") {
          # MH sampler for rho
          if (fix_rho_TF) {
            # Fixing rho (i.e. no accept/reject step; e.g. fix_rho_TF == TRUE)
            parms_current[[i]][[x]] <- parms_try[[x]]
            parms_accept_rate[[i]][x] <- NA
          } else {
            # Sampling for rho (e.g. fix_rho_TF == FALSE)
            # 1. Draw a candidate from the proposal ------------------------ ####
            parms_try[[x]] <- fn.draw[[parms_names[x]]](parms_current[[i]][[x]], 
                                                        tune_start[x],
                                                        rho_values)
            
            # 2. Calc dif in the log-likelihood ---------------------------- ####
            LLdif <- 
              (fn.llike[[parms_names[x]]](parms = parms_try, 
                                          moddat = moddat,
                                          Q = Q,
                                          U = U,
                                          W = W,
                                          z_index = z_index,
                                          y_index = y_index,
                                          rho_lookup = rho_lookup, 
                                          rho_values = rho_values) -
                 fn.llike[[parms_names[x]]](parms = parms_current[[i]],
                                            moddat = moddat,
                                            Q = Q,
                                            U = U,
                                            W = W,
                                            z_index = z_index,
                                            y_index = y_index,
                                            rho_lookup = rho_lookup, 
                                            rho_values = rho_values))
            
            # 3. Accept or reject update ----------------------------------- ####
            rU    <- log(runif(1))
            
            # Store sample value for this iteration
            if(rU < LLdif) {parms_current[[i]][[x]] <- parms_try[[x]]}
            parms_try[[x]] <- parms_current[[i]][[x]]
            
            # Store accept/reject decision
            parms_accept_rate[[i]][[i_thin]][[x]] <- rU < LLdif
            
          }
        } else if (parms_names[x] %in% c("sigma_z", "sigma_v", "v")) {
          # MH sampler for {sigma_z, sigma_v, v}
          # 1. Draw a candidate from the proposal ------------------------ ####
          parms_try[[x]] <- fn.draw[[parms_names[x]]](parms_current[[i]][[x]], 
                                                      tune_start[x],
                                                      rho_values)
          
          # 2. Calc dif in the log-likelihood ---------------------------- ####
          LLdif <- 
            (fn.llike[[parms_names[x]]](parms = parms_try, 
                                        moddat = moddat,
                                        Q = Q,
                                        U = U,
                                        W = W,
                                        z_index = z_index,
                                        y_index = y_index,
                                        rho_lookup = rho_lookup, 
                                        rho_values = rho_values) -
               fn.llike[[parms_names[x]]](parms = parms_current[[i]],
                                          moddat = moddat,
                                          Q = Q,
                                          U = U,
                                          W = W,
                                          z_index = z_index,
                                          y_index = y_index,
                                          rho_lookup = rho_lookup, 
                                          rho_values = rho_values))
          
          # 3. Accept or reject update ----------------------------------- ####
          rU    <- log(runif(1))
          
          # Store sample value for this iteration
          if(rU < LLdif) {parms_current[[i]][[x]] <- parms_try[[x]]}
          parms_try[[x]] <- parms_current[[i]][[x]]
          
          # Store accept/reject decision
          parms_accept_rate[[i]][[i_thin]][[x]] <- rU < LLdif
          
        } else {
          # This is the MH sampler for {z_ob, z_un}
          for (j in 1:nMCMC_z) {
            # 1. Draw a candidate from the proposal ---------------------- ####
            # updating a z_ob or z_un block...
            parms_try[[x]] <- fn.draw[[parms_names[x]]](parms_current[[i]][[x]], 
                                                        tune_start[x])
            
            # 2. Calc dif in the log-likelihood -------------------------- ####
            LLdif <- 
              (fn.llike[[parms_names[x]]](parms = parms_try, 
                                          moddat = moddat,
                                          Q = Q,
                                          U = U,
                                          W = W,
                                          z_index = z_index,
                                          y_index = y_index,
                                          rho_lookup = rho_lookup,
                                          rho_values = rho_values) -
                 fn.llike[[parms_names[x]]](parms = parms_current[[i]],
                                            moddat = moddat,
                                            Q = Q,
                                            U = U,
                                            W = W,
                                            z_index = z_index,
                                            y_index = y_index,
                                            rho_lookup = rho_lookup,
                                            rho_values = rho_values))
            
            # 3. Accept or reject update --------------------------------- ####
            rU    <- log(runif(1))
            
            # Store sample value for this iteration
            # - because we're referencing index "i" here, we'll continue 
            #   re-writing over this same spot, except for the last iter, which
            #   will stick and we'll move out of this for-loop.
            if(rU < LLdif) {parms_current[[i]][[x]] <- parms_try[[x]]}
            parms_try[[x]] <- parms_current[[i]][[x]]
            
            # Store accept/reject decision
            parms_accept_rate[[i]][[i_thin]][[x]][[j]] <- rU < LLdif
          }
        }
        }
    }
    # - Get ready to move to iter "i + 1" -------------------------------- ####
    # update where to draw from for the next iter...
    parms_current[[i + 1]] <- parms_current[[i]]
  }
  
  end_time <- Sys.time()
  
  # == Clean up output =================================================== ####
  # Create df with MCMC iters (rows) by parameters (cols)...
  # - note: by default our code fills in one extra iter of NA that we don't 
  #   need for parms_current - so drop that with this first line of code...
  
  # (acceptance rate: T/F whether proposed sample accepted by iter (n_sim))
  # - the 'acceptance rate' associated with a sim number here will be the 
  #   (mean) acceptance rate for the iterations leading up to that kept
  #   sim (e.g. if 75% thin, drop 3 and keep 4th sample, so acceptance rate 
  #   listed here (for say sim_n = 1), would be the average of those 4 T/F 
  #   values.)
  accept_rate <- parms_accept_rate %>% 
    # nMCMC level
    map(~.x %>% 
        map(~.x %>% map(~.x %>% unlist)) %>%
        transpose %>%
        map(~.x %>% unlist %>% mean)) %>%
    transpose %>%
    map(~.x %>% unlist) %>%
    bind_rows(.id = "n_sim")
  
  # (sample value by iter (n_sim))
  # - identify non-block update parms (e.g. not z_ob, z_un, or v)
  theta_names <- c(grep("beta_", parms_names, value = TRUE), "rho",
                   grep("sigma_", parms_names, value = TRUE))
  
  MCMC_theta <- parms_current[1:nMCMC] %>% map(~.x[theta_names])
  MCMC_theta <- do.call(rbind, MCMC_theta) %>% apply(2, unlist) %>% as.data.frame()
  MCMC_theta$n_sim <- 1:nMCMC
  
  # (matrix of z: MCMC iters (rows) by fiahex_id (cols))
  # - reassemble z vecotr from blocks
  MCMC_z <- parms_current[1:nMCMC] %>% map(function(x) {
    # re-assemble z vector from z_ob- and z_un-blocks
    z <- rep(NA, nrow(W))
    z[as.logical(z_index)]  <- x[["z_ob"]]
    z[!as.logical(z_index)] <- x[["z_un"]]
    
    return(z)
  })
  MCMC_z <- do.call(rbind, MCMC_z) %>% as.data.frame
  names(MCMC_z) <- colnames(U)
  
  # (matrix of v: MCMC iters (rows) by ES (cols))
  MCMC_v <- parms_current[1:nMCMC] %>% map(~.x$v)
  MCMC_v <- do.call(rbind, MCMC_v) %>% as.data.frame
  names(MCMC_v) <- colnames(Q)
  
  # == Return items of interest ========================================== ####
  return(list(run_time = end_time - start_time,
              accept_rate = accept_rate,
              MCMC_theta = MCMC_theta,
              MCMC_z = MCMC_z,
              MCMC_v = MCMC_v,
              parms_start = parms_start,
              tune_start  = tune_start,
              settings = c(nMCMC = nMCMC,
                           nThin = nThin,
                           nMCMC_z = nMCMC_z)))
}

# == Summary & plotting ================================================== ####
# All of the fn's below work with SDM output from fn.CAR_fit...
# - fn: summaries, plots, & maps used (primarily) during tuning ---------- ####
# summary stats (works with mod.x)
# - acceptance rates
fn.ar <- function(mod) {
  # acceptance rates
  # Note:
  # - the 'acceptance rate' associated with a sim_n here will be the 
  #   (mean) acceptance rate for the iterations leading up to that kept
  #   sim (e.g. if 75% thin, drop 3 and keep 4th sample, so acceptance rate 
  #   listed here (for say sim_n = 1), would be the average of those 4 T/F 
  #   values.)
  # - because each of those sim_n mean values are averaging the same number of
  #   T/F values, we can then take the mean of those means as below...
  mod$accept_rate %>% 
    select(-n_sim) %>%
    apply(2, function(x) round(mean(x), 2))
}

# - estimates
# (theta)
fn.est_theta <- function(mod) {
  # parms ests
  mod$MCMC_theta %>% 
    select(-n_sim) %>%
    apply(2, function(x) round(mean(x), 2))
}
# (z & v)
fn.est_zv  <- function(mod, z_or_v.chr, by_iter_T = FALSE) {
  # about:
  # - mod is an output from fn.CAR_fit
  # - z_or_v.chr is a string (either "z" or "v")
  
  if (z_or_v.chr == "z") {
    if (by_iter_T) {
      mod$MCMC_z %>% 
        apply(1, function(x) {
          data.frame(iter_mean = mean(x),
                     iter_median = median(x),
                     iter_IQR = IQR(x))
        }) %>% 
        bind_rows %>% 
        apply(2, summary)
    } else {
      mod$MCMC_z %>% 
        apply(2, function(x) {
          # notes: 
          # - `rle` stands for 'run length encoding'
          # - nonupdate of 1 means updated every time
          # - update of (n_sim - 1) means updated every time
          data.frame(
            mean = mean(x),
            iqr  = IQR(x),
            n_changes = sum(ifelse(x == lag(x), 0, 1), na.rm = TRUE),
            longest_nonupdate = rle(x)$lengths %>% max,
            avg_nonupdate = rle(x)$lengths %>% mean,
            longest_update = with(rle((x != lag(x))[-1]), 
                                  tapply(lengths, values, max))["TRUE"] %>% unname())
        }) %>% 
        bind_rows(.id = "fiahex_id") %>%
        mutate(fiahex_id = gsub("fiahex_id", "", fiahex_id))
    }
  } else if (z_or_v.chr == "v") {
    if (by_iter_T) {
      mod$MCMC_v %>%
        apply(1, function(x) {
          data.frame(iter_mean = mean(x),
                     iter_median = median(x),
                     iter_IQR = IQR(x))
        }) %>% 
        bind_rows %>% 
        apply(2, summary)
    } else {
      mod$MCMC_v %>% 
        apply(2, function(x) {
          # notes: 
          # - `rle` stands for 'run length encoding'
          # - nonupdate of 1 means updated every time
          # - update of (n_sim - 1) means updated every time
          data.frame(
            mean = mean(x),
            iqr  = IQR(x),
            n_changes = sum(ifelse(x == lag(x), 0, 1), na.rm = TRUE),
            longest_nonupdate = rle(x)$lengths %>% max,
            avg_nonupdate = rle(x)$lengths %>% mean,
            longest_update = with(rle((x != lag(x))[-1]), 
                                  tapply(lengths, values, max))["TRUE"] %>% unname())
        }) %>% 
        bind_rows(.id = "ES") %>%
        mutate(ES = gsub("ES", "", ES))
    }
  } else {"Not a valid z_or_v.chr (e.g. z or v chr)"}
}

# - return tibble with tuning params, accept rates, and theta ests
fn.mod_sum <- function(mod) {
  m_name <- deparse(substitute(mod))
  list(accept_rate = fn.ar(mod), 
       param_est = fn.est_theta(mod),
       tune = mod$tune_start) %>% 
    map(bind_rows) %>% 
    bind_rows(.id = "summary_type") %>%
    mutate(m_run = m_name)
}

# - identify z_i & v_i to take a deeper look at
# (this mostly used with fn.plot_trace_zv)
fn.trace_z_ids <- function(moddat) {
  # get the fiahex_id for hex's with these characteristics 
  hex <- moddat %>% 
    group_by(fiahex_id) %>%
    summarize(p_plots = sum(p_a, na.rm = TRUE),
              a_plots = sum(!p_a, na.rm = TRUE),
              n_plots = n(),
              y_unobsv = sum(!y_index)) %>%
    ungroup
  
  # pure hex (all observed)
  x_out <- c(many_p.ob = (hex %>% 
                            filter(y_unobsv == 0) %>% 
                            filter(p_plots == max(p_plots)))$fiahex_id[1])
  x_out <- c(x_out,
             many_a.ob = (hex %>% 
                            filter(y_unobsv == 0 & !fiahex_id %in% x_out) %>% 
                            filter(a_plots == max(a_plots)))$fiahex_id[1])
  x_out <- c(x_out,
             many_mix.ob = (hex %>% 
                              filter(y_unobsv == 0 & !fiahex_id %in% x_out) %>% 
                              filter(a_plots > 1 & p_plots > 1) %>%
                              filter(n_plots == max(n_plots)))$fiahex_id[1])
  if(is.na(x_out["many_mix.ob"])) {
    x_out["many_mix.ob"] <- (hex %>% 
                               filter(y_unobsv == 0 & !fiahex_id %in% x_out) %>% 
                               filter(n_plots == max(n_plots)))$fiahex_id[1]
  }
  x_out <- c(x_out,
             few.ob = (hex %>% 
                         filter(y_unobsv == 0 & !fiahex_id %in% x_out) %>% 
                         filter(n_plots == min(n_plots)))$fiahex_id[1])
  
  # different characteristics if z_un exist...
  if (any(moddat$y_index == 0)) {
    x_out <- c(x_out,
               many.mix = (hex %>% 
                             filter((y_unobsv > 0) & 
                                      (y_unobsv < n_plots) & 
                                      (!fiahex_id %in% x_out)) %>% 
                             filter(n_plots == max(n_plots)))$fiahex_id[1])
    x_out <- c(x_out,
               un = (hex %>% 
                       filter(y_unobsv == n_plots) %>% 
                       filter(n_plots == max(n_plots)))$fiahex_id[1])
    
    return(x_out)
  } else {
    tmp <- (hex %>% filter(y_unobsv == 0 & !fiahex_id %in% x_out) %>% arrange(desc(n_plots)))$fiahex_id
    x_out <- c(x_out,
               many.ob = tmp[1],
               mid.ob = tmp[floor(length(tmp)/2)])
    
    return(x_out)
  }
}
fn.trace_v_ids <- function(moddat) {
  # get ES id's for ES with the most, the least, and a mid-range for 
  # - the total no. plots
  # - number of p_plots
  ES_p <- moddat %>%
    group_by(ES) %>%
    summarize(n_plots = n(),
              p_plots = sum(p_a, na.rm = TRUE)) %>%
    ungroup %>%
    mutate(ES = as.character(ES))
  
  
  # grab ES to represent p-plot position
  ES_by_p_plots <- (ES_p %>% arrange(desc(p_plots)))$ES
  
  x_out <- c(many_p_plots = ES_by_p_plots[1],
             mid_p_plots  = ES_by_p_plots[floor(length(ES_by_p_plots)/2)],
             few_p_plots  = ES_by_p_plots[length(ES_by_p_plots)])
  
  # filter these out so we don't get the same id's 2x
  # grab ES to represent n-plot position (outside of those ES already examined)
  ES_by_plots <- (ES_p %>% 
                    filter(!ES %in% x_out) %>% 
                    arrange(desc(n_plots)))$ES
  
  x_out <- c(x_out,
             many_plots = ES_by_plots[1],
             mid_plots  = ES_by_plots[floor(length(ES_by_plots)/2)],
             few_plots  = ES_by_plots[length(ES_by_plots)])
  
  x_out
}

# - plotting while tuning
fn.plot_trace_theta <- function(mod, parms_actual, plot_actual = TRUE, m_name = NULL) {
  if(is.null(m_name)) {m_name <- deparse(substitute(mod))}
  
  df_plot <- mod$MCMC_theta %>% 
    melt(id.vars = "n_sim", variable.name = "parameter")
  
  df_plot.quant <- df_plot %>%
    group_by(parameter) %>%
    summarize(lower = quantile(value, probs = .05),
              upper = quantile(value, probs = .95))
  
  df_plot.actual <- as.data.frame(t(parms_actual)) %>% 
    melt(variable.name = "parameter")
  
  p <-     ggplot(df_plot, aes(y = value, x = n_sim)) +
    facet_wrap(~parameter, scales = "free", ncol = 4) +
    geom_line() +
    geom_hline(data = df_plot.quant, aes(yintercept = lower), color = "blue") +
    geom_hline(data = df_plot.quant, aes(yintercept = upper), color = "blue") +
    ggtitle(paste0(m_name, ": MCMC trace & equal-tailed 90% credible interval")) +
    xlab("MCMC iteration")
  
  # - chain convergence
  if(plot_actual) {
    p + geom_hline(data = df_plot.actual, aes(yintercept = value), color = "red")
  } else {
    p
  }
}
fn.plot_density_theta <- function(mod, parms_actual, plot_actual = TRUE, m_name = NULL) {
  if(is.null(m_name)) {m_name <- deparse(substitute(mod))}
  
  df_plot <- mod$MCMC_theta %>% 
    melt(id.vars = "n_sim", variable.name = "parameter")
  
  df_plot.quant <- df_plot %>%
    group_by(parameter) %>%
    summarize(lower = quantile(value, probs = .05),
              upper = quantile(value, probs = .95))
  
  df_plot.actual <- as.data.frame(t(parms_actual)) %>% 
    melt(variable.name = "parameter")
  
  # - chain convergence
  p <-ggplot(df_plot, aes(x = value)) +
    facet_wrap(~parameter, scales = "free", ncol = 4) +
    # geom_density(aes(fill = parameter, color = parameter), alpha = .5) +
    geom_density(fill = "black", color = "black", alpha = .5) +
    geom_vline(data = df_plot.quant, aes(xintercept = lower), color = "blue") +
    geom_vline(data = df_plot.quant, aes(xintercept = upper), color = "blue") +
    # geom_rug() +
    ggtitle(paste0(m_name, ": posterior densities & equal-tailed 90% credible interval")) +
    # scale_color_viridis_d() +
    # scale_fill_viridis_d() +
    theme(legend.position = "none") 
  
  # add red-line with true value if simulation (plot_actual == TRUE)
  if(plot_actual) {
    p + geom_vline(data = df_plot.actual, aes(xintercept = value), color = "red")
  } else {
    p
  }
}
fn.plot_trace_zv <- function(mod, mod_name, zv_for_trace, z_or_v.chr, density_T = FALSE) {
  # grab the z_i (or v_i) corresponding to this hex (or ES)
  if (z_or_v.chr == "z") {
    df_plot<- mod$MCMC_z[paste0("fiahex_id", zv_for_trace)]
    df_plot$n_sim <- 1:nrow(df_plot)
    df_plot <- df_plot %>% 
      melt(id.vars = "n_sim") %>%
      mutate(variable = gsub("fiahex_id", "", variable)) %>%
      inner_join(data.frame(type = factor(names(zv_for_trace), levels = names(zv_for_trace)),
                            variable = zv_for_trace),
                 by = "variable") %>%
      mutate(variable = paste0(type, ": ", variable))
    
  } else if (z_or_v.chr == "v") {
    df_plot<- mod$MCMC_v[paste0("ES", zv_for_trace)]
    df_plot$n_sim <- 1:nrow(df_plot)
    df_plot <- df_plot %>% 
      melt(id.vars = "n_sim") %>%
      mutate(variable = gsub("ES", "", variable)) %>%
      inner_join(data.frame(type = factor(names(zv_for_trace), levels = names(zv_for_trace)),
                            variable = zv_for_trace),
                 by = "variable") %>%
      mutate(variable = paste0(type, ": ", variable))
  } else {
    "Not a valid z_or_v.chr (e.g. z or v chr)"
  }
  
  # for quantile lines..
  df_plot.quant <- df_plot %>%
    group_by(variable) %>%
    summarize(lower = quantile(value, probs = .05),
              upper = quantile(value, probs = .95))
  
  if (!density_T) {
    # trace plots for these hex
    ggplot(df_plot, aes(y = value, x = n_sim)) +
      facet_wrap(~variable, scales = "free", nrow = 2) +
      geom_line() +
      geom_hline(data = df_plot.quant, aes(yintercept = lower), color = "blue") +
      geom_hline(data = df_plot.quant, aes(yintercept = upper), color = "blue") +
      ggtitle(paste0(mod_name, ": ", z_or_v.chr,
                     " MCMC trace & equal-tailed 90% credible interval"))
  } else {
    # density plots for these hex
    ggplot(df_plot, aes(x = value)) +
      facet_wrap(~variable, scales = "free", nrow = 2) +
      geom_density(fill = "black", color = "black", alpha = .5) +
      geom_vline(data = df_plot.quant, aes(xintercept = lower), color = "blue") +
      geom_vline(data = df_plot.quant, aes(xintercept = upper), color = "blue") +
      ggtitle(paste0(mod_name, ": ", z_or_v.chr,
                     " posterior densities & equal-tailed 90% credible interval"))
  }
  
}
fn.plot_density_zv_mean <- function(mod, m_name, z_or_v.chr) {
  # grab the z_i (or v_i) corresponding to this hex (or ES)
  if (z_or_v.chr == "z") {
    # get means of samples for each z or v
    df_plot <- mod$MCMC_z %>%
      apply(2, function(x) mean(x)) %>%
      bind_rows %>%
      melt(variable.name = "fiahex_id",
           value.name = "z_mean") %>%
      mutate(fiahex_id = gsub("fiahex_id", "", fiahex_id))
    
    # get deciles for geom_rug (strip plot)
    df_deciles <- data.frame(decile = df_plot$z_mean %>% quantile(c(0:10)*.1))
    
    ggplot(data = df_plot, aes(x = z_mean)) +
      # geom_histogram(aes(y=..density..),
      #                bins = 200,
      #                color = "black",
      #                fill = "darkgrey",
      #                alpha = .5) +
      geom_density(aes(y=..density..), alpha = .2, fill = "#2D708EFF") +
      geom_rug(data = df_deciles, aes(decile)) + 
      # geom_point(aes(y = 0), alpha = .05, shape = 73, size = 6) +
      ggtitle(paste0(m_name, ": distribution of ", z_or_v.chr, "_i means"),
              subtitle = "split plot of deciles, max, and min values")
    
  } else if (z_or_v.chr == "v") {
    # get means of samples for each z or v
    df_plot <- mod$MCMC_v %>%
      apply(2, function(x) mean(x)) %>%
      bind_rows %>%
      melt(variable.name = "ES",
           value.name = "v_mean") %>%
      mutate(ES = gsub("ES", "", ES))
    
    # get deciles for geom_rug (strip plot)
    df_deciles <- data.frame(decile = df_plot$v_mean %>% quantile(c(0:10)*.1))
      
    ggplot(data = df_plot, aes(x = v_mean)) +
      # geom_histogram(aes(y=..density..),
      #                binwidth = .25,
      #                color = "black",
      #                fill = "darkgrey",
      #                alpha = .5) +
      geom_density(aes(y=..density..), alpha = .2, fill = "#2D708EFF") +
      geom_rug(data = df_deciles, aes(decile)) + 
      # geom_point(aes(y = 0), alpha = .4, shape = 73, size = 6) +
      ggtitle(paste0(m_name, ": distribution of ", z_or_v.chr, "_i means"),
              subtitle = "split plot of deciles, max, and min values")
    
  } else {
    "Not a valid z_or_v.chr (e.g. z or v chr)"
  }
}
fn.plot_ar_theta <- function(mod, m_name = NULL) {
  if(is.null(m_name)) {m_name <- deparse(substitute(mod))}
  
  # fix n_sim & melt for faceting by parameter
  df_plot <- mod$accept_rate
  df_plot$n_sim <- 1:nrow(df_plot)
  df_plot <- df_plot %>% melt(id.vars = "n_sim", variable.name = "parameter")
  
  df_plot.quant <- df_plot %>%
    group_by(parameter) %>%
    summarize(lower = quantile(value, probs = .05),
              upper = quantile(value, probs = .95))

  # plot acceptance rates for each parameter by iteration
    ggplot(df_plot, aes(y = value, x = n_sim)) +
      facet_wrap(~parameter, scales = "free") +
      geom_line() +
      geom_hline(data = df_plot.quant, aes(yintercept = lower), color = "blue") +
      geom_hline(data = df_plot.quant, aes(yintercept = upper), color = "blue") +
      ggtitle(m_name, subtitle = "lines correspond to equal-tailed 90% interval")
}

# - sf: z_maps (true & public) 
fn.map_z_sample_sf <- function(mod.x, hex.x, sp.x, 
                               study_area.sf,
                               MCMC_iter,
                               expect_T = TRUE,
                               color_mu_T = TRUE) {
  # About: this function maps the in-sample hex-polygons (hex_polygons_*.sf)
  #        and colors (fills) by a sample (MCMC_iter) of the z's 
  # - expect_T: should the mean of the z's be mapped (T) or should 
  #   a sample of z's from the last iteration of the model run be mapped (F)
  # - color_mu: should hexagons be colored based on the scale of the mean 
  #   (T; e.g. mu_hat), or on the scale of the linear predictor (F; e.g. lp_hat)
  
  # extract type of z's to be used
  z.df <- switch(expect_T + 1,
                 mod.x$MCMC_z[MCMC_iter, ],
                 mod.x$MCMC_z %>% apply(2, mean) %>% bind_rows)
  
  # associate the last z sample with their centroid coordinates
  z.df <- z.df %>%
    melt %>%
    rename(fiahex_id = variable,
           z_value = value) %>%
    mutate(fiahex_id = gsub("fiahex_id", "", fiahex_id))
  
  # atm z's on lp scale, change to mu-scale if specified
  if (color_mu_T) {z.df <- z.df %>% mutate(z_value = fn.inv_logit(z_value))}
  
  # everything is in UTM11N...
  hex.x <- hex.x %>%
    mutate(fiahex_id = as.character(fiahex_id)) %>%
    left_join(z.df, by = "fiahex_id")
  
  # make custom titles
  # - plot title
  p_title <- switch(expect_T + 1, 
                    paste0(sp.x, ": z sample from MCMC iter", MCMC_iter),
                    paste0(sp.x, ": mean z based on ?k MCMC iters"))
  # - legend title
  l_title <- switch(color_mu_T + 1,
                    "lp-scale",
                    "mu-scale")
  
  # plot it!
  ggplot() +
    geom_sf(data = study_area.sf, fill = "white") + # plot first so points are on top
    geom_sf(data = hex.x,
            aes(fill = z_value),
            color = NA) +
    scale_fill_viridis_c(l_title)+
    theme_bw()  +
    ggtitle(p_title)
}

# - fn: plots & maps used when analyzing results (scripts 9+) ------------ ####
# sp-climate relationships based on conditional probability:
# - returns summaries, figures, and raw optimal data
#   - mean climate conditions for peak conditional probability
#   - the mode of (optimal) climate conditions for peak conditional probability
fn.cond_spclim_rel <- function(MCMC_theta.x, dat_ls.x, dat_std.x, covar_lims.x, 
                               crosswalk.x, 
                               grid_size = 100, long_name_var,
                               long_names = FALSE) {
  # About                                                                  ####
  # - this function returns plots of conditional predicted probabilities,
  #   which are predicted by holding random effect (z, v) and any uninvolved
  #   climate covariates at zero (prior means & std means, respectively)
  # - for plotting, the argument `grid_size` is used to one dim of seq or 
  #   both dims for square grid
  # - `long_names = TRUE` will set plot axis names to Ecosphere styled names
  #   (e.g. 'ppt-wetQ (mm)' for bio16)
  
  # = Set up ============================================================= ####
  # - define univariate mode grabbing fn --------------------------------- ####
  fn.modes_univariate <- function(x) {
    x <- round(x, 1)
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)] 
  }
  
  # - identify univariate or bivariate relationships --------------------- ####
  # For each beta parameter, add a label to crosswalk.x indicating the 
  # highest-level variable-form its covariate is associated with: 
  # - bivariate relationships: interaction
  # - univariate relationships: quadratic (for bio* terms) or linear (gdd5)
  
  # Make temporary df with labels to add to crosswalk.x....
  # - intercept (repeat 3 times so will match up with different beta-types)
  beta_type <- data.frame(orig_names = rep("intercept", 3), 
                          beta_ind = c("inter", "quadlin", "lin"))
  
  # - betas involved with interaction terms
  if (any(grepl("_b", crosswalk.x$orig_names))) {
    beta_type <- crosswalk.x$orig_names %>%
      # names of covariates involved in interaction terms
      grep("_b", ., value = TRUE) %>%
      str_split(., "_") %>% 
      unlist %>%
      # vector of variable forms (L, Q, I) associated with these covariates
      sapply(function(x) grep(x, crosswalk.x$orig_names, value = TRUE)) %>% 
      as.vector %>%
      unique %>%
      # make into a df we can update beta_type with
      data.frame(orig_names = ., 
                 beta_ind = "inter") %>%
      rbind.data.frame(beta_type, .)
  }
  # - all other bio* vars (not involved in interaction terms)
  if (any(grepl("bio", (crosswalk.x %>% 
                        filter(!orig_names %in% beta_type$orig_names) %>% 
                        pull(orig_names))))) {
    beta_type <- beta_type %>% rbind.data.frame(
      data.frame(orig_names = crosswalk.x %>% 
                   filter(!orig_names %in% beta_type$orig_names) %>%
                   pull(orig_names) %>%
                   grep("bio", ., value = TRUE),
                 beta_ind = "quadlin"))
  }
  # - non-bio vars (i.e. gdd5 for us) 
  if (nrow(crosswalk.x %>% 
           filter(!orig_names %in% beta_type$orig_names)) > 0) {
    beta_type <- beta_type %>% rbind.data.frame(
      data.frame(orig_names = crosswalk.x %>% 
                   filter(!orig_names %in% beta_type$orig_names) %>%
                   pull(orig_names),
                 beta_ind = "lin"))
  }
  
  # Add to crosswalk.x and delete
  crosswalk.x <- crosswalk.x %>% inner_join(beta_type, by = "orig_names")
  rm(beta_type)
  
  # drop non-used beta_ind terms for intercept
  # - b/c not all model forms have all types remove leftover 'intercept' terms
  crosswalk.x <- crosswalk.x %>%
    group_by(beta_ind) %>%
    mutate(n_beta = n()) %>%
    ungroup %>%
    filter(n_beta > 1)
  
  # - initiate list to hold plots & data by relationship ----------------- ####
  n_plots <- sum(grepl("_b", crosswalk.x$orig_names)) +
    sum(crosswalk.x[!crosswalk.x$mod_name == "beta_0", ]$beta_ind == "quadlin")/2 +
    sum(crosswalk.x[!crosswalk.x$mod_name == "beta_0", ]$beta_ind == "lin")
  
  # create plot-holding list 
  cond_plots <- vector("list", length = n_plots)
  
  # - set index for list to 1 -------------------------------------------- ####
  plot_ind <- 1 # plot set's index
  
  # = Bivariate relationships ============================================ ####
  if (any(crosswalk.x$beta_ind == "inter")) {
    # identify types of interactions (b/c there can be more than one)
    inter_base <- grep("_b", crosswalk.x$orig_names, value = TRUE)
    
    # for each interaction, develop the corresponding var-plot
    cond_plots[1:length(inter_base)] <- lapply(inter_base, function(x) {
      # - what variables are at the core of this interaction? ------------ ####
      base.x  <- str_split(x, "_") %>% unlist
      
      # - make grid of variable forms for prediction --------------------- ####
      # Note: we want to visualize on the original scale, so seq evenly along
      #   the original scale...
      
      # - transform bounds back to original scale
      range.x <- covar_lims.x %>% 
        # grab the right stats and nonstd values
        filter(stat_type %in% c("max_value", "min_value")) %>%
        filter(std_type == "nonstd") %>% 
        # grab the right base-vars & drop unneeded cols
        filter(clim_var %in% base.x) %>%
        select(-scenario, -std_type) %>%
        # get climate bounds for simulation
        group_by(clim_var, stat_type) %>%
        summarize(value = ifelse(stat_type == "max_value", 
                                 max(value), 
                                 min(value))) %>%
        distinct(value) %>%
        ungroup %>%
        # recast into an easier format
        dcast(clim_var ~ stat_type, value.var = "value")
      range.x$expand_name <- c("Var1", "Var2")
      
      # - if bio* 16:19 exist, move back to monthly scale 
      #   (do this so that things work when re-std to pred below... will
      #   move it all back to total from monthly post prediction)
      adj_vars <- paste0("bio", 16:19)
      range.x <- range.x %>%
        mutate_at(c("max_value", "min_value"), 
                  ~ifelse(clim_var %in% adj_vars, .x/3, .x))
      
      # - expand grid on original scale
      df.x <- expand.grid(seq(range.x[1, 3], range.x[1, 2], length.out = grid_size),
                          seq(range.x[2, 3], range.x[2, 2], length.out = grid_size)) 
      
      # - transform back to std scale
      df.x <- df.x %>% 
        mutate(dummy = 1:n()) %>%
        melt(id.vars = "dummy",
             variable.name = "expand_name") %>%
        inner_join(range.x %>% select(clim_var, expand_name),
                   by = "expand_name") %>%
        select(-expand_name) %>%
        left_join(dat_std.x, by = "clim_var") %>%
        mutate(value_std = (value - std_mean) / std_sd) %>%
        select(dummy, clim_var, value_std) %>%
        dcast(dummy ~ clim_var, value.var = "value_std") %>%
        select(-dummy)
      
      # - get quad/inter forms
      df.x <- df.x %>%
        mutate(!!paste0(base.x[1], "_2") := (!!(sym(base.x[1])))^2,
               !!paste0(base.x[2], "_2") := (!!(sym(base.x[2])))^2,
               !!paste0(base.x[1], "_", base.x[2]) := (!!(sym(base.x[1]))) * (!!(sym(base.x[2]))))
      
      # - add predictions & inflection points by MCMC sample ------------- ####
      # - order df.x in same way as MCMC_theta 
      #   (e.g. beta_* sequentially in R, so bio13 is before bio8, etc.)
      vars.x <- crosswalk.x %>% 
        filter(orig_names %in% c(names(df.x), "intercept") & 
                 beta_ind == "inter")
      p.x <- nrow(vars.x)
      df.x <- df.x %>% 
        select(one_of(crosswalk.x %>% 
                        filter(orig_names %in% names(df.x)) %>%
                        pull(orig_names)))
      # - add preds by MCMC sample
      df.x <- c(1:nrow(MCMC_theta.x)) %>% 
        map(function(i) {
          parms.i <- (MCMC_theta.x %>% select(one_of(vars.x$mod_name)))[i, ] %>% unlist
          df_out <- df.x %>% select(all_of(base.x))
          df_out$mu_cond <- 
            (parms.i[1] + as.matrix(df.x[, 1:(p.x-1)]) %*% parms.i[2:p.x]) %>%
            as.vector %>% fn.inv_logit %>% as.vector
          df_out$MCMC_run <- i
          
          df_out %>%
            mutate(MCMC_run = i,       # add MCMC sample number
                   dummy_id = 1:n())   # 'plot_id', avoids later dcast agg.
        }) %>% 
        bind_rows
      
      # - get inflection points (including optimal conditions) by MCMC sample
      df_opt <- lapply(1:nrow(MCMC_theta.x), function(i) {
        parms.i <- (MCMC_theta.x %>% select(one_of(vars.x$mod_name)))[i, ] %>% unlist
        
        # store the type of inflection point...
        D <- (4*parms.i[4]*parms.i[5] - (parms.i[6])^2)
        fxx <- sign(2*parms.i[4])
        
        opt_type <- case_when(D == 0 ~ "unk",
                              D < 0 ~ "saddle",
                              (D > 0 & fxx < 0) ~ "max", # these are 'optimals'
                              (D > 0 & fxx > 0) ~ "min")
        
        opt_vec <- (-1/D *
                      matrix(c(2*parms.i[5], -parms.i[6], 
                               -parms.i[6], 2*parms.i[4]),
                             nrow = 2,
                             byrow = TRUE) %*% 
                      c(parms.i[2], parms.i[3])) %>% as.vector 
        
        data.frame(value = opt_vec,
                   opt_type = opt_type,
                   clim_var = sort(base.x), # ensure order (e.g. bio13, bio8)
                   MCMC_run = i)
      }) %>%
        bind_rows
      
      # - transform variables back to original scale --------------------- ####
      df.x <- df.x %>%
        melt(id.vars = c("mu_cond", "MCMC_run", "dummy_id"),
             variable.name = "clim_var") %>%
        left_join(dat_std.x, by = "clim_var") %>%
        mutate(value_nonstd = value * std_sd + std_mean) %>%
        dcast(mu_cond + MCMC_run + dummy_id ~ clim_var, 
              value.var = "value_nonstd") %>%
        select(-dummy_id)
      
      df_opt <- df_opt %>%
        left_join(dat_std.x, by = "clim_var") %>%
        mutate(value_nonstd = value * std_sd + std_mean) %>%
        dcast(MCMC_run + opt_type ~ clim_var, 
              value.var = "value_nonstd")
      
      # - make bio16 & bio18 for qtr total (not mean) -------------------- ####
      # ATM bio16 & bio18 are "mean monthly ppt in * quarter", but we want
      # to report "total ppt in * quarter", so these variables need to be
      # multiplied by 3 (3 months in a quarter)... bio17 & bio19 are also
      # quarterly ppt, but they don't appear in our models... however
      # include those below to be thorough.
      adj_vars <- paste0("bio", 16:19)
      
      df.x <- df.x %>% mutate_at(vars(any_of(adj_vars)), ~ .x*3)
      df_opt <- df_opt %>% mutate_at(vars(any_of(adj_vars)), ~ .x*3)
      
      # - get mean-predictions & prediction uncertainty (se) ------------- ####
      # data manipulation for the plotting methods below...
      df_plot.x <- df.x %>%
        group_by(across(all_of(base.x))) %>%
        summarize(se_mu_cond = sd(mu_cond),
                  mean_mu_cond = mean(mu_cond)) %>%
        ungroup
      
      # - mean conditions for peak conditional probability --------------- ####
      mean_cond_peak <- df_plot.x %>% filter(mean_mu_cond == max(mean_mu_cond))
      
      # - mode of optimal conditions (on non-std data) ------------------- ####
      # using ks::kde b/c MASS::kde2d H selector returns 0 for some of our
      # species' bandwidths... which means we can't estimate mode of the
      # optimum...
      
      # for kde, constrain optimums to range of datasets (current & future)
      # - search only within area where 2k optimums are estimated
      opt_bounds <- df_plot.x %>% 
        select(all_of(base.x)) %>% 
        apply(., 2, function(x) c(min(x), max(x))) %>% 
        as.data.frame
      
      df_opt_restricted <- df_opt %>%
        # drop any 'optimums' that are not a maximum
        filter(opt_type == "max") %>%
        # restrict to only values within our data sets...
        filter((!!as.symbol(base.x[1])) > opt_bounds[[base.x[1]]][1]) %>%
        filter((!!as.symbol(base.x[1])) < opt_bounds[[base.x[1]]][2]) %>%
        filter((!!as.symbol(base.x[2])) > opt_bounds[[base.x[2]]][1]) %>%
        filter((!!as.symbol(base.x[2])) < opt_bounds[[base.x[2]]][2])
      
      # constrain kde density to coords within the restricted optimums
      # - this way, if all optimums in a smaller area, we'll have a more 
      #   zoomed-in finer grid.. if they ranged widely this will still be
      #   restricting thins to the range of conditions with in our datasets
      opt_lims <- df_opt_restricted %>% 
        select(all_of(base.x)) %>% 
        apply(., 2, function(x) c(min(x), max(x))) %>% 
        as.data.frame %>% 
        unlist %>% 
        unname
      
      # calculate the density
      # - calls ks::Hpi by default for selection of H
      #   (which is good as these need a full bandwidth matrix)
      df_opt.kde <- ks::kde(x = df_opt_restricted %>% select(base.x),
                            gridsize = 151, # default
                            xmin = opt_lims[c(1,3)],
                            xmax = opt_lims[c(2,4)])
      
      df_opt.kde <- with(df_opt.kde,
                         data.frame(expand.grid(eval.points[[1]],
                                                eval.points[[2]]),
                                    as.vector(estimate))) %>%
        set_colnames(c(base.x, "density"))
      
      # clean-up space 
      rm(opt_lims)
      
      # - final plots: mean-predictions & SE & optimals ------------------ ####
      # create new x & y axis label shortcuts
      x_var <- base.x[1]
      y_var <- base.x[2]
      z_var <- "mean_mu_cond"
      se_z  <- "se_mu_cond"
      
      # set axis titles
      if (long_names) {
        x_lab <- long_name_var %>% 
          filter(orig_names == x_var) %>% 
          pull(both_names)
        y_lab <- long_name_var %>% 
          filter(orig_names == y_var) %>% 
          pull(both_names)
      } else {
        x_lab <- long_name_var %>% 
          filter(orig_names == x_var) %>% 
          mutate(short_name = paste0(orig_names, " (", units, ")")) %>% 
          pull(short_name)
        y_lab <- long_name_var %>% 
          filter(orig_names == y_var) %>% 
          mutate(short_name = paste0(orig_names, " (", units, ")")) %>% 
          pull(short_name)
      }
      
      # mean-pred plot
      p1_meanpred <- ggplot(df_plot.x) +
        geom_tile(aes_string(x = x_var, y = y_var, fill = z_var),
                  alpha = 0.9) +
        scale_fill_gradient("probability",
                            low = "black", high = "white") +
        xlab(x_lab) +
        ylab(y_lab) +
        theme_bw() + 
        theme(aspect.ratio = 1)
      
      # prediction uncertianty
      p1_se <- ggplot(df_plot.x) +
        geom_tile(aes_string(x = x_var, y = y_var, fill = se_z),
                  alpha = 0.9) +
        scale_fill_gradient(paste0("SE of\n", "probability"),
                            low = "black", high = "white") +
        xlab(x_lab) +
        ylab(y_lab) +
        theme_bw() + 
        theme(aspect.ratio = 1)
      
      # optimal conditions
      p1_opt <- ggplot(df_opt.kde) +
        geom_tile(aes_string(x = x_var, y = y_var, fill = "density"),
                  alpha = 0.9) +
        scale_fill_gradient(low = "black", high = "white") +
        xlab(x_lab) +
        ylab(y_lab) +
        theme_bw() + 
        theme(aspect.ratio = 1)
      
      # - make a list of info to return ---------------------------------- ####
      return(list(p1_meanpred = p1_meanpred,
                  p1_se = p1_se,
                  p1_opt = p1_opt,
                  mean_cond_peak = mean_cond_peak,
                  # raw 2k inflection (critical: saddle, min, max, etc.) points
                  opt_cond_mcmc = df_opt, 
                  # optimals within dataset bounds, used for mode calculation
                  opt_cond_restricted = df_opt_restricted, 
                  # bivariate mode
                  mode = df_opt.kde %>% filter(density == max(density))))
      
    }) %>% 
      setNames(inter_base)
    
    # set up the next plot-set's index
    plot_ind <- length(inter_base) + 1
  }
  
  # = Univariate relationships =========================================== ####
  # - covariates with linear & quadratic terms --------------------------- ####
  # any remaining bio* vars at this point will have a linear & quadratic forms
  if (any(crosswalk.x$beta_ind == "quadlin")) {
    # identify types of these vars (i.e. more than one?)
    quadlin_base <- crosswalk.x %>%
      filter(beta_ind == "quadlin") %>%
      pull(orig_names) %>%
      grep("_2", ., value = TRUE) %>%
      gsub("_2", "", .)
    
    # for each interaction, develop the corresponding var-plot
    cond_plots[plot_ind:(plot_ind + length(quadlin_base) - 1)] <- 
      lapply(quadlin_base, function(base.x) {
        # - make grid of variable forms for prediction ------------------- ####
        # Note: we want to visualize on the original scale, so seq evenly along
        #   the original scale...
        
        # - transform bounds back to original scale
        range.x <- covar_lims.x %>% 
          # grab the right stats and nonstd values
          filter(stat_type %in% c("max_value", "min_value")) %>%
          filter(std_type == "nonstd") %>% 
          # grab the right base-vars & drop unneeded cols
          filter(clim_var == base.x) %>%
          select(-scenario, -std_type) %>%
          # get climate bounds for simulation
          group_by(clim_var, stat_type) %>%
          summarize(value = ifelse(stat_type == "max_value", 
                                   max(value), 
                                   min(value))) %>%
          distinct(value) %>%
          ungroup %>%
          # recast into an easier format
          dcast(clim_var ~ stat_type, value.var = "value") 
        range.x$expand_name <- c("Var1")
        
        # - if bio* 16:19 exist, move back to monthly scale 
        #   (do this so that things work when re-std to pred below... will
        #   move it all back to total from monthly post prediction)
        adj_vars <- paste0("bio", 16:19)
        range.x <- range.x %>%
          mutate_at(c("max_value", "min_value"), 
                    ~ifelse(clim_var %in% adj_vars, .x/3, .x))
        
        
        # - expand grid on original scale
        df.x <- data.frame(Var1 = seq(range.x[1, 3], range.x[1, 2], length.out = grid_size))
        
        # - transform back to std scale
        df.x <- df.x %>% 
          mutate(dummy = 1:n()) %>%
          melt(id.vars = "dummy",
               variable.name = "expand_name") %>%
          inner_join(range.x %>% select(clim_var, expand_name),
                     by = "expand_name") %>%
          select(-expand_name) %>%
          left_join(dat_std.x, by = "clim_var") %>%
          mutate(value_std = (value - std_mean) / std_sd) %>%
          select(dummy, clim_var, value_std) %>%
          dcast(dummy ~ clim_var, value.var = "value_std") %>%
          select(-dummy)
        
        # - get quad forms
        df.x <- df.x %>%
          mutate(!!paste0(base.x, "_2") := (!!(sym(base.x)))^2)
        
        # - add predictions & inflection points by MCMC sample ----------- ####
        # - order df.x in same way as MCMC_theta (e.g. beta_* sequentially)
        vars.x <- crosswalk.x %>% 
          filter(orig_names %in% c(names(df.x), "intercept") & 
                   beta_ind == "quadlin")
        p.x <- nrow(vars.x)
        df.x <- df.x %>% 
          select(one_of(crosswalk.x %>% 
                          filter(orig_names %in% names(df.x)) %>%
                          pull(orig_names)))
        # - add preds by MCMC sample
        df.x <- c(1:nrow(MCMC_theta.x)) %>% 
          map(function(i) {
            parms.i <- (MCMC_theta.x %>% select(one_of(vars.x$mod_name)))[i, ] %>% unlist
            df_out <- df.x %>% select(all_of(base.x))
            df_out$mu_cond <- (parms.i[1] + as.matrix(df.x[, 1:(p.x-1)]) %*% parms.i[2:p.x]) %>% 
              as.vector %>% fn.inv_logit %>% as.vector
            
            df_out %>%
              mutate(MCMC_run = i,       # add MCMC sample number
                     dummy_id = 1:n())   # 'plot_id', avoids later dcast agg.
          }) %>% 
          bind_rows
        
        # - get inflection points (including optimal conditions) by MCMC sample
        df_opt <- lapply(1:nrow(MCMC_theta.x), function(i) {
          parms.i <- (MCMC_theta.x %>% select(one_of(vars.x$mod_name)))[i, ] %>% unlist
          fxx <- sign(2*parms.i[3])
          
          data.frame(value = -parms.i[2]/(2*parms.i[3]),
                     opt_type = case_when(fxx < 0 ~ "max", # these are 'optimals'
                                          fxx > 0 ~ "min",
                                          fxx == 0 ~ "unk"))
        }) %>% 
          bind_rows
        df_opt$clim_var <- base.x
        df_opt$MCMC_run <- 1:nrow(df_opt)
        
        # - transform variables back to original scale ------------------- ####
        df.x <- df.x %>%
          melt(id.vars = c("mu_cond", "MCMC_run", "dummy_id"),
               variable.name = "clim_var") %>%
          left_join(dat_std.x, by = "clim_var") %>%
          mutate(value_nonstd = value * std_sd + std_mean) %>%
          dcast(mu_cond + MCMC_run + dummy_id ~ clim_var, 
                value.var = "value_nonstd") %>%
          select(-dummy_id)
        
        df_opt <- df_opt %>%
          left_join(dat_std.x, by = "clim_var") %>%
          mutate(value_nonstd = value * std_sd + std_mean) %>%
          dcast(MCMC_run + opt_type ~ clim_var, 
                value.var = "value_nonstd")
        
        # - make bio16 & bio18 for qtr total (not mean) ------------------ ####
        # ATM bio16 & bio18 are "mean monthly ppt in * quarter", but we want
        # to report "total ppt in * quarter", so these variables need to be
        # multiplied by 3 (3 months in a quarter)... bio17 & bio19 are also
        # quarterly ppt, but they don't appear in our models... however
        # include those below to be thorough.
        adj_vars <- paste0("bio", 16:19)
        
        df.x <- df.x %>% mutate_at(vars(any_of(adj_vars)), ~ .x*3)
        df_opt <- df_opt %>% mutate_at(vars(any_of(adj_vars)), ~ .x*3)
        
        # - mean peak and optimal mode ----------------------------------- ####
        # get mean conditions for peak conditional probability
        df.x_meanpred <- df.x %>%
          group_by(across(all_of(base.x))) %>%
          summarize(se_mu_cond = sd(mu_cond),
                    mean_mu_cond = mean(mu_cond)) %>%
          ungroup
        
        # get mode of optimal conditions
        # - constrain optimums to range of datasets (current & future) for
        #   realistic mode & posterior optimum plot
        opt_bounds <- c(min(df.x_meanpred[[base.x]]), 
                        max(df.x_meanpred[[base.x]]))
        
        df_opt_restricted <- df_opt %>%
          # drop any 'optimums' that are not a maximum
          filter(opt_type == "max") %>%
          # restrict to only values within our data sets...
          filter((!!as.symbol(base.x[1])) > opt_bounds[1]) %>%
          filter((!!as.symbol(base.x[1])) < opt_bounds[2]) 
        
        # constrain kde density to coords within the restricted optimums
        # - this way, if all optimums in a smaller area, we'll have a more 
        #   zoomed-in finer grid.. if they ranged widely this will still be
        #   restricting thins to the range of conditions with in our datasets
        opt_lims <- c(min(df_opt_restricted[[base.x]] %>% floor), 
                      max(df_opt_restricted[[base.x]]) %>% ceiling)
        
        # calculate the density
        # - calls ks::Hpi by default for selection of H
        #   (which is good as these need a full bandwidth matrix)
        df_opt.kde <- ks::kde(x = df_opt_restricted[[base.x]],
                              gridsize = 151, # default
                              xmin = opt_lims[1],
                              xmax = opt_lims[2])
        
        df_opt.kde <- with(df_opt.kde,
                           data.frame(eval.points,
                                      as.vector(estimate))) %>%
          set_colnames(c(base.x, "density"))
        
        # clean-up space 
        rm(opt_lims)
        
        # - final plot and summaries ------------------------------------- ####
        # set axis titles
        x_lab <- ifelse(long_names,
                        long_name_var %>% 
                          filter(orig_names == base.x) %>% 
                          pull(both_names),
                        long_name_var %>% 
                          filter(orig_names == base.x) %>% 
                          mutate(short_name = paste0(orig_names, " (", units, ")")) %>% 
                          pull(short_name))
        y_lab <- "probability"
        
        p1_meanpred <- ggplot() +
          geom_line(data = df.x,
                    aes_string(x = base.x,
                               y = "mu_cond",
                               group = "MCMC_run"),
                    alpha = .05) + 
          geom_line(data = df.x_meanpred,
                    aes_string(x = base.x,
                               y = "mean_mu_cond"),
                    color = "red",
                    size = 1) + 
          xlab(x_lab) +
          ylab(y_lab) +
          theme_bw() + 
          theme(aspect.ratio = 1)
        
        p1_opt <- ggplot() +
          geom_density(data = df_opt_restricted, 
                       aes_string(base.x),
                       fill = "black", alpha = .5) +
          xlab(x_lab) +
          theme_bw() + 
          theme(aspect.ratio = 1)
        
        # - make a list of info to return ---------------------------------- ####
        return(list(p1_meanpred = p1_meanpred,
                    p1_opt = p1_opt,
                    mean_cond_peak = df.x_meanpred %>% filter(mean_mu_cond == max(mean_mu_cond)),
                    # raw 2k inflection (critical: saddle, min, max, etc.) points
                    opt_cond_mcmc = df_opt, 
                    # optimals within dataset bounds, used for mode calculation
                    opt_cond_restricted = df_opt_restricted, 
                    # mode
                    mode_old = df_opt_restricted[[base.x]] %>% fn.modes_univariate,
                    mode = (df_opt.kde %>% filter(density == max(density)))[[base.x]]))
        
      }) %>%
      setNames(quadlin_base)
    
    # set up next plot's index
    plot_ind <- length(quadlin_base) + plot_ind
  }
  
  # - covariates with only linear terms ---------------------------------- ####
  # any non-bio* vars (should be it!)
  if (any(crosswalk.x$beta_ind == "lin")) {
    # identify types of these vars (i.e. more than one?)
    lin_base <- crosswalk.x %>%
      filter(beta_ind == "lin" & mod_name != "beta_0") %>%
      pull(orig_names)
    
    # for each interaction, develop the corresponding var-plot
    cond_plots[plot_ind:(plot_ind + length(lin_base) - 1)] <- 
      lapply(lin_base, function(base.x) {
        # - make grid of variable forms for prediction ------------------- ####
        # Note: we want to visualize on the original scale, so seq evenly along
        #   the original scale...
        
        # - transform bounds back to original scale
        range.x <- covar_lims.x %>% 
          # grab the right stats and nonstd values
          filter(stat_type %in% c("max_value", "min_value")) %>%
          filter(std_type == "nonstd") %>% 
          # grab the right base-vars & drop unneeded cols
          filter(clim_var == base.x) %>%
          select(-scenario, -std_type) %>%
          # get climate bounds for simulation
          group_by(clim_var, stat_type) %>%
          summarize(value = ifelse(stat_type == "max_value", 
                                   max(value), 
                                   min(value))) %>%
          distinct(value) %>%
          ungroup %>%
          # recast into an easier format
          dcast(clim_var ~ stat_type, value.var = "value") 
        range.x$expand_name <- c("Var1")
        
        # - expand grid on original scale
        df.x <- data.frame(Var1 = seq(range.x[1, 3], range.x[1, 2], length.out = grid_size))
        
        # - transform back to std scale
        df.x <- df.x %>% 
          mutate(dummy = 1:n()) %>%
          melt(id.vars = "dummy",
               variable.name = "expand_name") %>%
          inner_join(range.x %>% select(clim_var, expand_name),
                     by = "expand_name") %>%
          select(-expand_name) %>%
          left_join(dat_std.x, by = "clim_var") %>%
          mutate(value_std = (value - std_mean) / std_sd) %>%
          select(dummy, clim_var, value_std) %>%
          dcast(dummy ~ clim_var, value.var = "value_std") %>%
          select(-dummy)
        
        # - add predictions by MCMC run ---------------------------------- ####
        # - order df.x in same way as MCMC_theta (e.g. beta_* sequentially)
        vars.x <- crosswalk.x %>% 
          filter(orig_names %in% c(names(df.x), "intercept") & 
                   beta_ind == "lin")
        p.x <- nrow(vars.x)
        df.x <- df.x %>% 
          select(one_of(crosswalk.x %>% 
                          filter(orig_names %in% names(df.x)) %>%
                          pull(orig_names)))
        # - add preds by MCMC sample
        df.x <- c(1:nrow(MCMC_theta.x)) %>% 
          map(function(i) {
            parms.i <- (MCMC_theta.x %>% select(one_of(vars.x$mod_name)))[i, ] %>% unlist
            df_out <- df.x %>% select(all_of(base.x))
            df_out$mu_cond <- (parms.i[1] + as.matrix(df.x[, 1:(p.x-1)]) %*% parms.i[2:p.x]) %>% 
              as.vector %>% fn.inv_logit %>% as.vector
            
            df_out %>%
              mutate(MCMC_run = i,       # add MCMC sample number
                     dummy_id = 1:n())   # 'plot_id', avoids later dcast agg.
          }) %>% 
          bind_rows
        
        # - note: no true optimal conditions b/c only linear variable form
        
        # - transform variables back to original scale ------------------- ####
        df.x <- df.x %>%
          melt(id.vars = c("mu_cond", "MCMC_run", "dummy_id"),
               variable.name = "clim_var") %>%
          left_join(dat_std.x, by = "clim_var") %>%
          mutate(value_nonstd = value * std_sd + std_mean) %>%
          dcast(mu_cond + MCMC_run + dummy_id ~ clim_var, 
                value.var = "value_nonstd") %>%
          select(-dummy_id)
        
        # final plot and summaries --------------------------------------- ####
        # set axis titles
        x_lab <- ifelse(long_names,
                        long_name_var %>% 
                          filter(orig_names == base.x) %>% 
                          pull(both_names),
                        long_name_var %>% 
                          filter(orig_names == base.x) %>% 
                          mutate(short_name = paste0(orig_names, " (", units, ")")) %>% 
                          pull(short_name))
        y_lab <- "probability"
        
        df.x_meanpred <- df.x %>%
          group_by(across(all_of(base.x))) %>%
          summarize(se_mu_cond = sd(mu_cond),
                    mean_mu_cond = mean(mu_cond)) %>%
          ungroup
        
        p1_meanpred <- ggplot() +
          geom_line(data = df.x,
                    aes_string(x = base.x,
                               y = "mu_cond",
                               group = "MCMC_run"),
                    alpha = .05) + 
          geom_line(data = df.x_meanpred,
                    aes_string(x = base.x,
                               y = "mean_mu_cond"),
                    color = "red",
                    size = 1) + 
          xlab(x_lab) +
          ylab(y_lab) +
          theme_bw() + 
          theme(aspect.ratio = 1)
        
        return(list(p1_meanpred = p1_meanpred,
                    mean_cond_peak = df.x_meanpred %>% 
                      filter(mean_mu_cond == max(mean_mu_cond))))
      }) %>%
      setNames(lin_base)
  }
  
  # = send back the results ============================================== ####
  # switch instead of ifelse to be able to return 'NULL'...
  cond_plot_names <- c(switch(exists("inter_base") + 1, NULL, inter_base),
                       switch(exists("quadlin_base") + 1, NULL, quadlin_base),
                       switch(exists("lin_base") + 1, NULL, lin_base))
  # returned named vector of plots
  cond_plots %>% setNames(cond_plot_names)
}

# mean-prediction maps & uncertainty maps in 2x5 grid for a given species
# - note, this uses cowplot fns, but calling like `cowplot::` to avoid issues
#   with patchwork, which is also commonly used in these scripts.
fn.map_mp_uncert <- function(preds.x, sp_name, plot_sf.x,
                             study_area.sf, 
                             T_fnodiff_used = FALSE,
                             EJ_options_check = TRUE) {
  # About:
  # Generate a 2 (mean-predictions & uncertainty) x 5 (scenario) grid of maps
  # for each species. Generated with facet_grid (1x5) for mean-predictions & 
  # uncertainty and then these two figures joined via cowplot::plot_grid().
  # 
  # Note width/ht/size opts were originally hard-coded, updated this fn but 
  # still should check that objects like EJ_dpi have been created (should be
  # created earlier in the scripts in which this function is called, e.g., 
  # 10...R (these are Fig requirements for the journal).
  
  # - Check for 'EJ_' objects -------------------------------------------- ####
  EJ_options_chr <- c("EJ_axis_text_size",
                      "EJ_strip_text_size",
                      "EJ_axis_title_size",
                      "EJ_label_size",
                      "EJ_dim_units",
                      "EJ_dim_1col_width_max",
                      "EJ_dim_2col_width_max",
                      "EJ_dim_2col_ht_max",
                      "EJ_app_dim_units",
                      "EJ_app_dim_width_max",
                      "EJ_app_dim_ht_max",
                      "EJ_dpi")
  
  if (!all(EJ_options_chr %in% ls(envir = .GlobalEnv))) {
    stop(paste0("Missing object:\n",
                paste0(EJ_options_chr[!EJ_options_chr %in% ls(envir = .GlobalEnv)],
                       collapse = "\n")))
  }
  
  # - make the base figures ---------------------------------------------- ####
  tmp_margin <- margin(1, 0, 0, 0)
  leg_margin <- margin(0, -7, 0, -9)
  future_lab_levels <- c("RCP 4.5\nCCSM4",
                         "RCP 4.5\nHadGEM2-ES",
                         "RCP 8.5\nCCSM4",
                         "RCP 8.5\nHadGEM2-ES")
  # Note: using pe = mean-predictions & pu = uncertainty plots
  preds.x <- plot_sf.x %>%
    select(-p_a, -NAME) %>%
    right_join(preds.x, by = "plot_id") %>%
    mutate(scenario = case_when(
      scenario == "current" ~ "current\n",
      scenario == "CCSM4_rcp45_m2085" ~ future_lab_levels[1],
      scenario == "HadGEM2_ES_rcp45_m2085" ~ future_lab_levels[2],
      scenario == "CCSM4_rcp85_m2085" ~ future_lab_levels[3],
      scenario == "HadGEM2_ES_rcp85_m2085" ~ future_lab_levels[4]))
  
  future.df <- preds.x %>%
    filter(scenario != "current\n") %>%
    mutate(scenario = factor(scenario, levels = future_lab_levels))
  
  # current part of the figure ----
  pe_current <- (preds.x %>% filter(scenario == "current\n")) %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(data = study_area.sf, fill = "grey50") +
    geom_sf(aes(fill = mu_hat), color = NA)  +
    scale_fill_viridis_c("prob.", option = "viridis",
                         limits = c(0, 1),
                         breaks = c(0.0, 0.5, 1.0)) +
    facet_grid(. ~ scenario) +
    scale_x_continuous(breaks = c(-122, -118)) +
    theme_bw() +
    theme(plot.margin = tmp_margin,
          axis.title = element_blank(),
          strip.text.x = element_text(size = EJ_strip_text_size),
          axis.text  = element_text(size = EJ_axis_text_size,
                                    color = "black"),
          legend.title = element_text(size = EJ_axis_title_size),
          legend.text  = element_text(size = EJ_strip_text_size),
          legend.key.height = unit(4, "mm"),
          legend.key.width = unit(3, "mm"),
          legend.box.margin = leg_margin)
  pu_current <- (preds.x %>% filter(scenario == "current\n")) %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(data = study_area.sf, fill = "grey50") +
    geom_sf(aes(fill = se_mu), color = NA)  +
    scale_fill_viridis_c("SE", option = "magma") +
    facet_grid(. ~ scenario) +
    scale_x_continuous(breaks = c(-122, -118)) +
    theme_bw()  +
    theme(plot.margin = tmp_margin,
          axis.title = element_blank(),
          strip.text.x = element_text(size = EJ_strip_text_size),
          axis.text  = element_text(size = EJ_axis_text_size,
                                    color = "black"),
          legend.title = element_text(size = EJ_axis_title_size),
          legend.text  = element_text(size = EJ_strip_text_size),
          legend.key.height = unit(4, "mm"),
          legend.key.width = unit(3, "mm"),
          legend.box.margin = leg_margin)
  
  # future part of figure ----
  if (T_fnodiff_used) {
    pe_future <- future.df %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, fill = "grey50") +
      geom_sf(aes(fill = mu_hat), color = NA)  +
      scale_fill_viridis_c("prob.", option = "viridis",
                           limits = c(0, 1),
                           breaks = c(0.0, 0.5, 1.0)) +
      facet_grid(. ~ scenario) +
      scale_x_continuous(breaks = c(-122, -118,)) +
      theme_bw() +
      theme(plot.margin = tmp_margin, 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            strip.text.x = element_text(size = EJ_strip_text_size),
            axis.text.x  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
            legend.title = element_text(size = EJ_axis_title_size),
            legend.text  = element_text(size = EJ_strip_text_size),
            legend.key.height = unit(4, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.box.margin = leg_margin) 
    pu_future <- future.df %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, fill = "grey50") +
      geom_sf(aes(fill = se_mu), color = NA)  +
      scale_fill_viridis_c("SE", option = "magma") +
      facet_grid(. ~ scenario) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() +
      theme(plot.margin = tmp_margin, 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            strip.text.x = element_text(size = EJ_strip_text_size),
            axis.text.x  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
            legend.title = element_text(size = EJ_axis_title_size),
            legend.text  = element_text(size = EJ_strip_text_size),
            legend.key.height = unit(4, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.box.margin = leg_margin) 
  } else {
    pe_future <- future.df %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, fill = "grey50") +
      geom_sf(aes(fill = mu_hat), color = NA)  +
      scale_fill_gradient2("diff prob.",
                           low = "#39568cff",
                           mid = "white",
                           high = "#FDE725FF",
                           midpoint = 0,
                           limits = c(-1, 1),
                           breaks = c(-1.0, 0.0, 1.0)) +
      facet_grid(. ~ scenario) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() +
      theme(plot.margin = tmp_margin, 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            strip.text.x = element_text(size = EJ_strip_text_size),
            axis.text.x  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
            legend.title = element_text(size = EJ_axis_title_size),
            legend.text  = element_text(size = EJ_strip_text_size),
            legend.key.height = unit(4, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.box.margin = leg_margin) 
    pu_future <- future.df %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, fill = "grey50") +
      geom_sf(aes(fill = se_mu), color = NA)  +
      scale_fill_gradient2("diff SE",
                           low = "#451077FF",
                           mid = "white",
                           high = "#EB5760FF",
                           midpoint = 0,
                           limits = (c(-1, 1)*max(abs(future.df$se_mu)))) +
      facet_grid(. ~ scenario) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() +
      theme(plot.margin = tmp_margin, 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            strip.text.x = element_text(size = EJ_strip_text_size),
            axis.text.x  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
            legend.title = element_text(size = EJ_axis_title_size),
            legend.text  = element_text(size = EJ_strip_text_size),
            legend.key.height = unit(4, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.box.margin = leg_margin) 
  }
  
  # - join the two figures into one -------------------------------------- ####
  pe_row <- cowplot::plot_grid(pe_current, pe_future,
                               nrow = 1, axis = "bt", rel_widths = c(1.4, 4))
  pu_row <- cowplot::plot_grid(pu_current, pu_future,
                               nrow = 1, axis = "bt", rel_widths = c(1.4, 4))
  # p_body <- 
  cowplot::plot_grid(pe_row, pu_row,
                     nrow = 2, axis = "l", rel_heights = c(1, 1)) + 
    theme(plot.background=element_rect(fill="white", color=NA))
  
  # # make joint title & subtitle
  # p_title <- ifelse(T_fnodiff_used,
  #                   paste0(sp_name, 
  #                          ": mean-predictions & uncertianty",
  #                          "; future values on original probability scale"),
  #                   paste0(sp_name, 
  #                          ": mean-predictions & uncertianty",
  #                          "; future values differenced from current"))
  # 
  # p_title <- cowplot::ggdraw() +
  #   cowplot::draw_label(label = p_title, x = 0, hjust = 0) +
  #   theme(plot.margin = margin(0, 0, 0, 7))
  # 
  # # add title to the joint plot & return
  # cowplot::plot_grid(p_title, p_body, ncol = 1, rel_heights = c(0.025, 1))
  
}

# quantile maps (range-size est with VP method) w. filled polygons by mu_hat
# - note, this uses cowplot fns, but calling like `cowplot::` to avoid issues
#   with patchwork, which is also commonly used in these scripts.
fn.map_rangesize <- function(preds.x, sp_name.x, plot_sf.x, study_area.sf,
                             diff_map_T = FALSE,
                             EJ_options_check = TRUE) {
  # About:                                                                -----
  # Makes a 2x5 paneled figure, where each panel shows the map cooresponding
  # to the upper or lower quantile (row) of range-size for a given scenario
  # (column); insert on each map is the range-size estimate (thousand km^2)
  # for that MCMC sample.
  # - preds.x: element from `map_base` which contains only the MCMC samples
  #   corresponding to the upper and lower quantiles (atm, set throughout
  #   to the 5th and 95th percentiles)
  # - sp_name.x: species code (lower case)
  # - plot_sf.x: public polygons (sf) used for plotting plot-level info
  # - study_area.sf: outline of OR/WA/CA used for plotting maps
  # - diff_map_T: logical; if 'FALSE' for all maps to be on the prob-scale 
  #   or 'TRUE' for future-maps to be vizualized as the difference from the 
  #   current-map
  # 
  # Note width/ht/size opts were originally hard-coded, updated this fn but 
  # still should check that objects like EJ_dpi have been created (should be
  # created earlier in the scripts in which this function is called, e.g., 
  # 10...R (these are Fig requirements for the journal).
  
  
  # - Check for 'EJ_' objects -------------------------------------------- ####
  EJ_options_chr <- c("EJ_axis_text_size",
                      "EJ_strip_text_size",
                      "EJ_axis_title_size",
                      "EJ_label_size",
                      "EJ_dim_units",
                      "EJ_dim_1col_width_max",
                      "EJ_dim_2col_width_max",
                      "EJ_dim_2col_ht_max",
                      "EJ_app_dim_units",
                      "EJ_app_dim_width_max",
                      "EJ_app_dim_ht_max",
                      "EJ_dpi")
  
  if (!all(EJ_options_chr %in% ls(envir = .GlobalEnv))) {
    stop(paste0("Missing object:\n",
                paste0(EJ_options_chr[!EJ_options_chr %in% ls(envir = .GlobalEnv)],
                       collapse = "\n")))
  }
  
  # - set legend & percentile labels ------------------------------------- ####
  # legend labels
  new_scen_names <- c("current\n",
                      "RCP 4.5\n CCSM4",
                      "RCP 4.5\n HadGEM2-ES",
                      "RCP 8.5\n CCSM4",
                      "RCP 8.5\n HadGEM2-ES")
  old_scen_names <- c("current",
                      "CCSM4_rcp45_m2085",
                      "HadGEM2_ES_rcp45_m2085",
                      "CCSM4_rcp85_m2085",
                      "HadGEM2_ES_rcp85_m2085")
  
  # percentile labels
  q_type_levels <- preds.x[[1]] %>% 
    map(~.x$q_type[1]) %>% 
    unlist %>% 
    sort(decreasing = TRUE) %>% # so upper quantile on top
    gsub("q_", "", .) %>% 
    paste0(.,"th percentile")
  
  # - make base dfs ------------------------------------------------------ ####
  preds.x <- preds.x %>% map_depth(2, function(df.x) {
    # now working within a scenario and then within the MCMC sample coor. to 
    # one of the quantiles
    df.x %>%
      # convert the range-size estimate to thousand km^2
      mutate(range_size = (range_size_ha/100/1000) %>% round(1)) %>%
      # add right label/levels for percentiles
      mutate(q_type = paste0(gsub("q_", "", q_type), "th percentile")) %>%
      mutate(q_type = factor(q_type, levels = q_type_levels)) %>%
      # add new label/levels for scenarios
      mutate(scenario = case_when(
        scenario == "current" ~ "current\n",
        scenario == old_scen_names[2] ~ new_scen_names[2],
        scenario == old_scen_names[3] ~ new_scen_names[3],
        scenario == old_scen_names[4] ~ new_scen_names[4],
        scenario == old_scen_names[5] ~ new_scen_names[5])) %>%
      select(plot_id, mu_hat, scenario, range_size, q_type)
  })
  plot_sf.x <- plot_sf.x %>% select(plot_id, geometry)
  
  # - Make figure -------------------------------------------------------- ####
  if (diff_map_T) {
    # make current part of the figure
    p_current <- (preds.x$current %>%
                    bind_rows %>%
                    inner_join(plot_sf.x, by = "plot_id")) %>%
      ggplot(aes(geometry = geometry))+
      geom_sf(data = study_area.sf, fill = "grey50") + 
      geom_sf(aes(fill = mu_hat), color = NA) +
      # place label with range size over northern NV-ish...
      geom_sf_label(data = . %>% 
                      distinct(scenario, q_type, .keep_all = TRUE) %>%
                      mutate(x = 491178, y = 4555448) %>% # these are UTM11N coords
                      st_as_sf(., 
                               coords = c("x", "y"), 
                               crs = st_crs(study_area.sf)),
                    aes(label = range_size),
                    size = (5/14) * 7) +
      facet_grid(q_type ~ scenario,
                 labeller = labeller(q_type = NULL)) +
      scale_fill_viridis_c("prob.", option = "viridis",
                           limits = c(0, 1),
                           breaks = c(0.0, 0.5, 1.0)) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.margin = margin(5.5,0,0,0),
            axis.title = element_blank(),
            strip.text = element_text(size = EJ_strip_text_size),
            axis.text  = element_text(size = EJ_axis_text_size,
                                      color = "black"),
            legend.title = element_text(size = EJ_axis_title_size),
            legend.text  = element_text(size = EJ_strip_text_size),
            legend.key.height = unit(3, "mm"),
            legend.key.width = unit(4, "mm"),
            legend.box.spacing = unit(0.1, "mm"))
    
    # make future part of the figure
    p_future <- (list(preds.x$current,
                      preds.x[names(preds.x) != "current"] %>% transpose) %>%
                   # work within the maps for either the upper or lower quantile
                   pmap(function(df.c, ls.f) {
                     ls.f %>% map(function(df.f) {
                       df.f %>%
                         inner_join(df.c %>% select(plot_id, mu_hat), by = "plot_id") %>%
                         mutate(pred_diff = mu_hat.x - mu_hat.y) %>%
                         select(plot_id, pred_diff, range_size, q_type, scenario)
                     }) %>%
                       bind_rows
                   }) %>%
                   bind_rows %>%
                   inner_join(plot_sf.x, by = "plot_id")) %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = study_area.sf, fill = "grey50") + 
      geom_sf(aes(fill = pred_diff), color = NA) +
      # place label with range size over northern NV-ish...
      geom_sf_label(data = . %>% 
                      distinct(scenario, q_type, .keep_all = TRUE) %>%
                      mutate(range_size = if_else(range_size == 0, "< 0", as.character(range_size))) %>%
                      mutate(x = 491178, y = 4555448) %>% # these are UTM11N coords
                      st_as_sf(., 
                               coords = c("x", "y"), 
                               crs = st_crs(study_area.sf)),
                    aes(label = range_size),
                    color = "#666666",
                    size = (5/14) * 7) +
      facet_grid(q_type ~ scenario) +
      scale_fill_gradient2("diff prob.",
                           low = "#39568cff",
                           mid = "white",
                           high = "#FDE725FF",
                           midpoint = 0,
                           limits = c(-1, 1),
                           breaks = c(-1.0, 0.0, 1.0)) +
      scale_x_continuous(breaks = c(-122, -118)) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.margin = margin(5.5,0,0,0), 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            strip.text.x = element_text(size = EJ_strip_text_size),
            strip.text.y = element_text(size = EJ_strip_text_size),
            axis.text.x  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
            legend.title = element_text(size = EJ_axis_title_size),
            legend.text  = element_text(size = EJ_strip_text_size),
            legend.key.height = unit(3, "mm"),
            legend.key.width = unit(4, "mm"),
            legend.box.spacing = unit(0.1, "mm"))
    
    # make build & return joint figure
    # - tweaked with these settings until they'd align top & bottom... 
    #   but it took a while & rel_widths()...
    # p_body <- 
    cowplot::plot_grid(p_current, p_future,
                       nrow = 1, axis = "bt", rel_widths = c(1.4, 4)) + 
      theme(plot.background=element_rect(fill="white", color=NA))
    
    # make joint title & subtitle
    # p_title <- cowplot::ggdraw() +
    #   cowplot::draw_label(
    #     label = paste0(toupper(sp_name.x), 
    #                    ": upper & lower range-size estimates (thousand km^2; VP method)",
    #                    "; future values differenced from current"),
    #     x = 0,
    #     hjust = 0) +
    #   theme(plot.margin = margin(0, 0, 0, 7))
    # 
    # # add title to the joint plot & return
    # cowplot::plot_grid(p_title, p_body, ncol = 1, rel_heights = c(0.025, 1))
    
  } else {
    (preds.x %>%
       map(~.x %>% bind_rows) %>%
       bind_rows %>%
       inner_join(plot_sf.x, by = "plot_id")) %>%
      ggplot(aes(geometry = geometry))+
      geom_sf(data = study_area.sf) + 
      geom_sf(aes(fill = mu_hat), color = NA) +
      # place label with range size over northern NV-ish...
      geom_sf_label(data = . %>% 
                      distinct(scenario, q_type, .keep_all = TRUE) %>%
                      mutate(range_size = if_else(range_size == 0, "< 0", as.character(range_size))) %>%
                      mutate(x = 491178, y = 4555448) %>% # these are UTM11N coords
                      st_as_sf(., 
                               coords = c("x", "y"), 
                               crs = st_crs(study_area.sf)),
                    aes(label = range_size)) +
      facet_grid(q_type ~ scenario) +
      scale_fill_viridis_c("probability", option = "viridis") +
      scale_x_continuous(breaks = c(-122, -118, -114)) +
      theme_bw() +
      ggtitle(paste0(sp_name.x, ": maps corresponding to upper & lower range-size estimates (thousand km^2)"),
              subtitle = "based on 2k MCMC samples and VP method")
  }
}

# (old) boxplot or density plot of range-size estimates for VP method
fn.plot_rangesize <- function(range_size_ls, quant_range, 
                              plot_type = "density",
                              scen_colPal = NULL) {
  # - make new legend labels --------------------------------------------- ####
  legend_values <- with(scen_colPal, hex_cb_gp %>% setNames(scenario))
  legend_labels <- legend_values %>%
    names %>%
    gsub("_m2085", "", .) %>% 
    gsub("_rcp", " rcp", .) %>% 
    str_split(., " ") %>% 
    lapply(function(x) paste0(x[2], " ", x[1])) %>% 
    unlist %>% gsub("NA ", "", .)
  names(legend_values) <- legend_labels
  
  # - get base data set -------------------------------------------------- ####
  df <- range_size_ls %>%
    imap(function(sp.x, sp_name) {
      sp.x %>%
        filter(map_group == "maps_2k") %>%
        select(-map_group) %>%
        rowwise %>%
        mutate(scenario = legend_labels[grep(scenario, scen_colPal$scenario)]) %>% 
        mutate(scenario = factor(scenario, levels = legend_labels)) %>%
        mutate(range_ksqkm = (range_size_ha/100)/1000) %>% # convert from ha to 1k km^2
        select(-range_size_ha) %>%
        ungroup %>%
        mutate(spp = sp_name)
    }) 
  
  q_bars <- df %>%
    imap(function(sp.x, sp_name) {
      sp.x %>%
        split(.$scenario) %>%
        map(function(scen.x) {
          df.x <- density(scen.x$range_ksqkm)
          df.x <- with(df.x, data.frame(range_ksqkm = x, dens = y))
          q_bars_id <- (nrow(df.x) * quant_range) %>% round(0)
          # return the top (x,y) coord of quantile's line segment...
          # (other end of segment will be at (x,0))
          df.x <- df.x %>% arrange(desc(range_ksqkm))
          df.x <- df.x[q_bars_id, ] %>%
            filter(range_ksqkm >= 0)
        }) %>%
        bind_rows(.id = "scenario") %>%
        mutate(spp = sp_name)
    })
  
  if (plot_type == "boxplot") {
    df <- df %>% map(function(sp.x) {
      sp.x %>%
        group_by(scenario, spp) %>%
        summarize(q_lower  = quantile(range_ksqkm, min(quant_range)) %>% as.numeric,
                  q_upper  = quantile(range_ksqkm, max(quant_range)) %>% as.numeric,
                  avg_mean = mean(range_ksqkm),
                  y_min    = min(range_ksqkm),
                  y_max    = max(range_ksqkm)) %>%
        ungroup %>%
        mutate(scenario = factor(gsub(" ", "\n", scenario), 
                                 levels = gsub(" ", "\n", legend_labels)))
    })
    legend_values <- legend_values %>% setNames(gsub(" ", "\n", legend_labels))
  }
  
  # - plot it! ----------------------------------------------------------- ####
  std_margin <- margin(0, 0, 5.5, 0)
  
  switch(grep(plot_type, c("boxplot", "density"))[1],
         # for custom boxplot (list of plots by spp)
         {
           df %>%
             imap(function(sp.x, sp_name) {
               if (sp_name == "abpr") {
                 # keeps x label, no guide
                 ggplot(sp.x) +
                   geom_boxplot(aes(x = scenario, fill = scenario,
                                    ymin = log1p(y_min), 
                                    ymax = log1p(y_max),
                                    lower = log1p(q_lower), 
                                    upper = log1p(q_upper), 
                                    middle = log1p(avg_mean)),
                                position=position_dodge(width = 0),
                                stat = "identity",
                                width = .5,
                                alpha = .9) +
                   scale_fill_manual(values = legend_values) +
                   facet_grid(spp ~ ., scales = "free_y") +
                   # scale_y_continuous(position = "right") +
                   ylab("log1p(1k km^2)") +
                   theme_bw() +
                   theme(plot.margin = std_margin,
                         axis.title.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.text.x = element_blank()) + 
                   # patchwork workaround with guides instead of theme to remove legend
                   guides(color = "none", fill = "none") 
               } else if (sp_name == "quke") {
                 # keeps x label and guide
                 ggplot(sp.x) +
                   geom_boxplot(aes(x = scenario, fill = scenario,
                                    ymin = y_min, 
                                    ymax = y_max,
                                    lower = q_lower, 
                                    upper = q_upper, 
                                    middle = avg_mean),
                                position=position_dodge(width = 0),
                                stat = "identity",
                                width = .5,
                                alpha = .9) +
                   facet_grid(spp ~ ., scales = "free_y") +
                   scale_fill_manual(values = legend_values) +
                   # scale_y_continuous(position = "right") +
                   ylab("1k km^2") +
                   theme_bw() +
                   guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                   theme(plot.margin = std_margin) + 
                   guides(color = "none", fill = "none")
               } else {
                 # no x label and no guide
                 ggplot(sp.x) +
                   geom_boxplot(aes(x = scenario, fill = scenario,
                                    ymin = y_min, 
                                    ymax = y_max,
                                    lower = q_lower, 
                                    upper = q_upper, 
                                    middle = avg_mean),
                                position=position_dodge(width = 0),
                                stat = "identity",
                                width = .5,
                                alpha = .9) +
                   facet_grid(spp ~ ., scales = "free_y") +
                   scale_fill_manual(values = legend_values) +
                   # scale_y_continuous(position = "right") +
                   ylab("1k km^2") +
                   theme_bw() +
                   theme(plot.margin = std_margin,
                         axis.title.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.text.x = element_blank()) + 
                   guides(color = "none", fill = "none")
               }
             })
         },
         # for density (list of plots by spp)
         {
           df %>%
             imap(function(sp.x, sp_name) {
               if (sp_name == "abpr") {
                 ggplot(sp.x %>% filter(spp == "abpr")) +
                   geom_density(aes(x = log1p(range_ksqkm), 
                                    y = log1p(..density..), 
                                    fill = scenario),
                                alpha = .5,
                                trim = TRUE) +
                   # geom_segment(data = q_bars[[sp_name]],
                   #              aes(yend = log1p(0),
                   #                  xend = log1p(range_ksqkm),
                   #                  y = log1p(dens),
                   #                  x = log1p(range_ksqkm),
                   #                  color = scenario),
                   #              alpha = .8) +
                   facet_grid(spp ~ ., scales = "free_y") +
                   scale_color_manual(values = legend_values, guide = FALSE) +
                   scale_fill_manual(values = legend_values,
                                     labels = gsub(" ", "\n", legend_labels))  +
                   xlab("log1p(total range size)") +
                   ylab("log1p(density)") +
                   theme_bw() +
                   theme(plot.margin = std_margin) + 
                   # patchwork workaround with guides instead of theme to remove legend
                   guides(color = "none", fill = "none") 
               } else if (sp_name == "quke") {
                 ggplot(sp.x %>% filter(spp != "abpr")) +
                   geom_density(aes(x = range_ksqkm, 
                                    y = ..density.., 
                                    fill = scenario),
                                alpha = .5,
                                trim = TRUE) +
                   # geom_segment(data = q_bars[[sp_name]],
                   #              aes(yend = 0,
                   #                  xend = range_ksqkm,
                   #                  y = dens,
                   #                  x = range_ksqkm,
                   #                  color = scenario),
                   #              alpha = .8) +
                   facet_grid(spp ~ ., scales = "free_y") +
                   scale_color_manual(values = legend_values, guide = FALSE) +
                   scale_fill_manual(values = legend_values,
                                     labels = gsub(" ", "\n", legend_labels))  +
                   xlab("total range size") +
                   ylab("density") +
                   theme_bw() +
                   guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                   theme(plot.margin = std_margin,
                         # legend.position = "bottom",
                         # legend.box.margin = margin(3, 0, 0, 0
                         legend.title = element_blank())
               } else {
                 ggplot(sp.x %>% filter(spp != "abpr")) +
                   geom_density(aes(x = range_ksqkm, 
                                    y = ..density.., 
                                    fill = scenario),
                                alpha = .5,
                                trim = TRUE) +
                   # geom_segment(data = q_bars[[sp_name]],
                   #              aes(yend = 0,
                   #                  xend = range_ksqkm,
                   #                  y = dens,
                   #                  x = range_ksqkm,
                   #                  color = scenario),
                   #              alpha = .8) +
                   facet_grid(spp ~ ., scales = "free_y") +
                   scale_color_manual(values = legend_values, guide = FALSE) +
                   scale_fill_manual(values = legend_values,
                                     labels = gsub(" ", "\n", legend_labels))  +
                   xlab("total range size") +
                   ylab("density") +
                   theme_bw() +
                   theme(plot.margin = std_margin,
                         axis.title.x = element_blank()) + 
                   guides(color = "none", fill = "none")
               }
             })
         })
  
}

# == Data simulation function ============================================ ####
# fn.sim_data: simulate spatially autocorrelated data 
# - simulates data sets with similar structure (oddities) to the real data 
#   used in this project
# - can be used to check sampler code while in development
fn.sim_data <- function(grid_x, grid_y,
                        parms_actual, n_RE_groups, RE_group_probs,
                        perc_y_unobsv, plot_to_hex_rate,
                        perc_remove) {
  # About: ----
  # This code simulates data sets with similar structure (oddities) to the 
  # real data used in this project. The function arguments include:
  # - grid_x & grid_y: 
  #   Number of rows (x) and columns (y) the data can occur within. These square
  #   grid cells are equivalent to the hexagon FIA polygons and represent the 
  #   level at which I am specifying the spatial random effects (aka the z's)
  # - parms_actual: 
  #   Named vector of actual parameter values (names below)
  # - n_RE_groups:
  #   Number of groups for the random effect, equivalent to my Ecological 
  #   Sections (ES) that indicate plot-level membership.
  #   - a note here: observations are recorded at the plot-level, membership
  #     to an ES is also assigned at the plot level (e.g. not the hex-level).
  #     this means that two plots can occur in the same grid cell but have still
  #     have different ES group membership
  # - RE_group_probs:
  #   Probability of membership to each RE group -- included to simulate data
  #   that does not have equal-sized RE groups (e.g. some of mine have thousands
  #   of plots, while others have less than 10...)
  # - perc_y_unobserv:
  #   The percentage of the response that will be unobserved/missing; because
  #   of how this code is set up, the actual percent unobserved will be close
  #   but not exact to this percentage (to be specified as between 0-1)
  # - plot_to_hex_rate: 
  #   Average number of plots that occur in a single hexagon or grid cell; to
  #   be similar to my data this should be quite low (1.2 - 2.0)
  # - perc_remove:
  #   The percentage of plots to remove from the grid before returning the final
  #   simulated data. Added this feature to represent the un-even pattern of 
  #   number of neighbors in my data set (e.g. not all hex's have 6 neighbors)
  
  # == Pull apart parms_actual into individual elements ================== ####
  beta_0   <- parms_actual["beta_0"]
  beta_1   <- parms_actual["beta_1"]
  beta_2   <- parms_actual["beta_2"]
  beta_3   <- parms_actual["beta_3"]
  rho     <- parms_actual["rho"]
  sigma_z <- parms_actual["sigma_z"]
  sigma_v <- parms_actual["sigma_v"]
  
  # == Simulate covariate data and RE-group membership =================== ####
  # Here start by filling out the entire grid, then we'll slowly chop bits away
  
  # max grid size
  n_grid <- grid_x * grid_y
  
  # generate the grid
  moddat <- data.frame(xcoord = rep(1:grid_x, times = grid_y, each = 1),
                       ycoord = rep(1:grid_y, times = 1, each = grid_x))
  
  # thin-out the number of grid cells via perc_remove
  moddat <- moddat[sample(c(1:nrow(moddat)),
                          size = n_grid*(1-perc_remove)), ]
  
  # add hex (grid cell) id
  moddat$fiahex_id <- as.character(1:nrow(moddat)) # char for later fns
  
  # generate covariate data
  # - the `*5` because we're going to later thin out to the plot_to_hex_rate
  moddat <- rbind.data.frame(moddat, moddat, moddat, moddat, moddat)
  
  moddat$bio1 <- rnorm(nrow(moddat))
  moddat$bio2 <- rnorm(nrow(moddat))
  moddat$bio3 <- rnorm(nrow(moddat))
  
  # thin out number of plots in each grid cell
  moddat <- moddat[sample(x = c(1:nrow(moddat)), 
                          size = round(n_distinct(moddat$fiahex_id)*plot_to_hex_rate)), ]
  
  # add RE group
  moddat$ES <- sample(x = paste0("ES_", c(1:n_RE_groups)),
                      size = nrow(moddat),
                      replace = TRUE,
                      prob = RE_group_probs)
  
  # add observed/unobserved classification to each plot
  moddat$y_index <- sample(x = c(0,1), 
                           size = nrow(moddat), 
                           prob = c(perc_y_unobsv, 1 - perc_y_unobsv),
                           replace = TRUE)
  
  # add observed/unobserved classification to each hex
  # - 0 if a grid cell only contains plots that are unobserved 
  # - 1 if all observed or a mix of ob/un
  moddat <- moddat %>% 
    group_by(fiahex_id) %>%
    mutate(z_index = ifelse(sum(y_index) == 0, 0, 1)) %>%
    ungroup
  
  # create plot_id and add rownames to moddat
  moddat$plot_id <- 1:nrow(moddat)
  moddat <- as.data.frame(moddat)
  rownames(moddat) <- moddat$plot_id
  
  # == Get W matrix & remove isolated grid cells ========================= ####
  # Do this in a few steps:
  # - identify neighboring hex for each fiahex_id (not plot_id)
  # - identify isolated hex's (i.e. those without any neighbors)
  # - remove isolated hex's (and their associated plots) from data set
  # - re-identify hex neighbors with this clean data set & create sparse W
  #   matrix of neighbor relations from that...
  # Need to remove isolated hex's so that W (and therefore Sigma of the 
  # spatial RE) is positive definite.
  
  # Part 0: grab relavent data & set n_nb and nn_max_dist ---------------- ####
  # grab coords for fiahex_id's (e.g. not plot_id's) 
  df.x <- moddat %>%
    select(fiahex_id, xcoord, ycoord, z_index) %>%
    distinct()
  
  # specify max dist for a nearest neighbor
  nn_max_dist <- 1
  
  # specify number of nbs (plus one for self)
  n_nb <- 4 + 1
  
  # Part 1: identify isolated hex's to drop ------------------------------ ####
  # - grab coord data with rownames -------------------------------------- ####
  xycoord <- data.frame(xcoord = df.x$xcoord,
                        ycoord = df.x$ycoord)
  rownames(xycoord) <- df.x$fiahex_id  # keep track of rownames
  
  # - get neighbor relations --------------------------------------------- ####
  # hex-grid and add one b/c nabor::knn() counts self
  knnout <- knn(xycoord, k = n_nb)
  
  # only keep the immediately adjacent neighbors sharing a side (and ditch self)
  nb_list <- lapply(1:nrow(df.x), function(i) {
    knnout$nn.idx[i, 2:n_nb][knnout$nn.dist[i, 2:n_nb] <= nn_max_dist]
  })
  
  names(nb_list) <- rownames(xycoord) # keep track of rownames
  
  # - identify isolated hexs & drop from df.x & moddat ------------------- ####
  # keep track isolated hexagon's fiahex_id's and remove their corresponding
  # plots from df.x 
  hex_isolated <- names(nb_list)[nb_list %>% map(~length(.x) == 0) %>% unlist]
  
  # resubset list to keep hexagons (if any isolated exist)
  # - wrapped in if statement b/c sometimes in sim there are no isolated cells
  if(!is.null(hex_isolated)) {
    df.x <- df.x[!df.x$fiahex_id %in% hex_isolated, ]
    moddat <- moddat[!moddat$fiahex_id %in% hex_isolated, ]
  }
  
  
  # Part 2: create W sparse matrix --------------------------------------- ####
  # - grab coord data with rownames -------------------------------------- ####
  xycoord <- data.frame(xcoord = df.x$xcoord,
                        ycoord = df.x$ycoord)
  rownames(xycoord) <- df.x$fiahex_id  # keep track of rownames
  
  # - get neighbor relations --------------------------------------------- ####
  # hex-grid and add one b/c nabor::knn() counts self
  knnout <- knn(xycoord, k = n_nb)
  
  # only keep the immediately adjacent neighbors sharing a side (and ditch self)
  nb_list <- lapply(1:nrow(df.x), function(i) {
    knnout$nn.idx[i, 2:n_nb][knnout$nn.dist[i, 2:n_nb] <= nn_max_dist]
  })
  
  names(nb_list) <- rownames(xycoord) # keep track of rownames
  
  # - create W (a sparse matrix of neighbor relations) ------------------- ####
  sparsei <- NULL # - sparsei: index of self
  sparsej <- NULL # - sparsej: index of neighbors
  
  for (i in 1:length(nb_list)) {
    sparsei <- c(sparsei, rep(i, times = length(nb_list[[i]])))
    sparsej <- c(sparsej, nb_list[[i]])
  }
  
  # this is our 'W' that we want out::
  W <- sparseMatrix(x = 1, i = sparsei, j = sparsej)
  
  # - create W_ls (to work with CAR fn) ---------------------------------- ####
  W_ls <- list(W = W,
               z_index = df.x$z_index %>% setNames(df.x$fiahex_id))
  
  # - clean up space ----------------------------------------------------- ####
  rm(sparsei, sparsej, knnout, nb_list, n_nb, nn_max_dist, 
     df.x, xycoord, hex_isolated, W)
  
  # == Get U & Q design matrices ========================================= ####
  # design matrix relating plot-level observations to spatial RE
  U <- sparse.model.matrix(
    object = plot_id ~ fiahex_id - 1, 
    data = moddat %>%
      mutate(fiahex_id = factor(fiahex_id, levels = names(W_ls$z_index))) %>%
      set_rownames(.$plot_id)) # b/c mutate drops them...
  
  # design matrix relating plot-level observations to non-spatial RE
  Q <- sparse.model.matrix(
    object = plot_id ~ ES - 1, 
    data = moddat %>% set_rownames(.$plot_id)) # b/c mutate drops them...
  
  # == Generate response data ============================================ ####
  # Create Sigma and get the chol decomp --------------------------------- ####
  Sigma <- (sigma_z^2) * 
    solve(Diagonal(x = as.vector(W_ls$W %*% rep(1, nrow(W_ls$W)))) - rho*W_ls$W)
  L <- t(chol(Sigma))  # transpose to get lower triangular matrix
  
  # machine precision checks
  # - checks
  Q_posdef <- sum(eigen(Sigma)$values < 0)  # check if positive definite
  Q_mperc  <- max(abs(L %*% t(L) - Sigma))  # check chol decomp
  # - outcomes
  if (Q_posdef != 0) stop("simulated Sigma not positive definite")
  print(paste("machine precision for chol decomp:", Q_mperc))
  
  # Simulate spatial & non-spatial RE; create linear model --------------- ####
  z_actual <- as.vector(L %*% rnorm(ncol(U)))   # for trouble shooting
  v_actual <- rnorm(n_RE_groups, 0, sigma_v)    # ...
  Uz <- as.vector(U %*% z_actual)
  Qv <- as.vector(Q %*% v_actual)
  
  inv_mu <- beta_0 + 
    beta_1*moddat$bio1 + beta_2*moddat$bio2 + beta_3*moddat$bio3 + Qv + Uz
  
  # Sample from the linear model ----------------------------------------- ####
  mu <- exp(inv_mu)/(1 + exp(inv_mu)) # logit link
  moddat$p_a <- rbinom(nrow(moddat), 1, mu)
  
  # Mask out unobserved plot_id's ---------------------------------------- ####
  moddat <- moddat %>% mutate(p_a = ifelse(y_index == 0, NA, p_a))
  
  # == Return needed objects ============================================= ####
  # Order correctly and re-add rownames (2xcheck) ------------------------ ####
  moddat <- moddat %>% select(starts_with("bio"), everything())
  rownames(moddat) <- moddat$plot_id # this should already be good, but why not
  
  # return our things
  return(list(moddat = moddat,
              W_ls = W_ls,
              U = U,
              Q = Q,
              z_actual = z_actual,
              v_actual = v_actual))
}

