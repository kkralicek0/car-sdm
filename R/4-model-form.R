###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 4-model-forms.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Spatial-Bayesian models project shifts in 
## suitable habitat for Pacific Northwest tree species under climate change' by
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
##
## About - this script:
## - identify species-specific model forms
##
## About - output:
## - ch2-4-modeling-data.Rdata
##   - create list with species-specific modeling data sets & mean/sd used for
##     standardizing variables
## - ch2-4-glmm_all_data_fits.rds; ch2-4-glmm_fits_*.rds; 
##   - location: path_mod
##   - intermediate step, saving results from fitting these model forms along
##     the way, as these take a while to fit (and don't want to have to refit)
##   - * is 2 through 10 (e.g., ch2-4-glmm_fits_10.rds)
##
###############################################################################
library(dplyr)        # for ...E V E R Y T H I N G...
library(stringr)      # for string help (e.g. str_detect in functions)
library(reshape2)     # melt() and dcast()
library(purrr)        # working with lists
library(lme4)         # glmer()
library(broom)        # for tidy-ing glm/glmer results...

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_mod    <- paste0(path_data, "det model forms/")
path_images <- paste0(path_top, "images/2-explore-data/")

# == Load data =========================================================== ####
load(paste0(path_data, "ch2-2-data.Rdata"))
spp <- names(FIA_spp)

# == Define formula-fitting fn =========================================== ####
# fn: specify lin/quad terms in formula ---------------------------------- ####
# idea: 
# - for ppt and temp we'd expect a quadratic species response curve, hmi too
# - we would not expect a quadratic form for gdd5
fn.quad_form <- function(biohmi_vars) {
  # silo gdd5 from other vars if it is there...
  bio2_vars <- biohmi_vars[!biohmi_vars == "gdd5"]
  gdd5_vars <- biohmi_vars[biohmi_vars == "gdd5"]
  
  # square the terms that need squaring
  lp  <- paste0(bio2_vars, collapse = " + ") # all the linear terms
  lp2 <- paste0("I(", bio2_vars, "^2)", collapse = " + ")
  
  # if gdd5 there add it in
  paste0(c(lp, lp2, gdd5_vars), collapse = " + ")
}

# fn: add interactions terms if they exist ------------------------------- ####
# we would expect the impact of ppt on a species to depend on temp
# (however we wouldn't expect say min temp to impact the effect of max temp)
fn.inter_form <- function(bio_vars) {
  inter1 <- ("bio10" %in% bio_vars) & ("bio18" %in% bio_vars)
  inter1 <- switch(inter1 + 1, NULL, "bio10:bio18")
  
  inter2 <- ("bio8" %in% bio_vars) & ("bio16" %in% bio_vars)
  inter2 <- switch(inter2 + 1, NULL, "bio8:bio16")
  
  paste0(c(inter1, inter2), collapse = " + ")
}
# fn: (quke-10) add 13:8 interaction terms ------------------------------- ####
fn.inter_form138 <- function(bio_vars) {
  inter1 <- ("bio10" %in% bio_vars) & ("bio18" %in% bio_vars)
  inter1 <- switch(inter1 + 1, NULL, "bio10:bio18")
  
  inter2 <- ("bio8" %in% bio_vars) & ("bio16" %in% bio_vars)
  inter2 <- switch(inter2 + 1, NULL, "bio8:bio16")
  
  inter3 <- ("bio8" %in% bio_vars) & ("bio13" %in% bio_vars)
  inter3 <- switch(inter3 + 1, NULL, "bio8:bio13")
  
  paste0(c(inter1, inter2, inter3), collapse = " + ")
}

# == Subset FIA_spp to vars based on physiology/theory =================== ####
# About:
# - specify vars of interest for each species based on physiology
# - var choice from "Silvics of NA" & "Devine et al. (2012)" tree profiles
# note:
# - approached 2 ways... ran for 'all vars' approach with two sets of vars
#   and compared ... first set omitted some extremes and pptann (also bio13
#   for quke)... second set added these back in... compared results and 
#   decided to go with second set of 'all vars' as our 'all vars', but kept
#   the first set here so I could recall what my process was...

# 'all vars' set 1 (not used in the end... but kept for reference)
# # - from ggpairs plots produced in script 3, decided to drop bio12 b/c of  
# #   direct relationship with hmi and high corr with bio18
# # - note: intuition is to also drop bio6 b/c of high corr with bio8/bio10/gdd5 
# #   for most spp... but will check this through model fitting b/c seems like
# #   something that could change with climate change...
# clim_vars  <- list(
#   "abpr"   = c("hmi", "gdd5", "bio5", "bio8", "bio10", "bio16", "bio18"),
#   "psmem"  = c("hmi", "gdd5", "bio5", "bio6", "bio8",  "bio10", "bio16",  "bio18"),
#   "qudo"   = c("hmi", "gdd5", "bio6", "bio8", "bio10", "bio13", "bio16", "bio18"),
#   "quga4"  = c("hmi", "gdd5", "bio6", "bio8", "bio10", "bio16", "bio18"),
#   "quke"   = c("hmi", "gdd5", "bio6", "bio8", "bio10", "bio16", "bio18"))

# 'all vars' set 2 
# - extreme events (max/min temp snaps; bio5, bio6)
# - annual ppt (bio12)
# - (for quke) bio13 to get at "wet feet" note in reference materials
clim_vars  <- list(
  "abpr"   = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio16", "bio18"),
  "psmem"  = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio16", "bio18"),
  "qudo"   = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio13", "bio16", "bio18"),
  "quga4"  = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio16", "bio18"),
  "quke"   = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio13", "bio16", "bio18"))

data_spp <- list(FIA_spp, clim_vars) %>% pmap(function(sp.x, clim.x) {
  sp.x %>% select("p_a", "ES", clim.x)
})

# == Standardize covariates ============================================== ####
# standardize these variables & keep list of mean/sd used to standardize
# - keeping mean/sd for 'standardizing' CNA future data for pred later on...
data_spp.std <- list(data_spp, clim_vars) %>% pmap(function(sp.x, clim.x) {
  dat <- sp.x %>% 
    mutate_at(clim.x, function(x) (x - mean(x)) / sd(x))
  
  mean_sd <- sp.x %>% 
    select(clim.x) %>% 
    melt(variable.name = "clim_var") %>%
    group_by(clim_var) %>%
    summarize(std_mean = mean(value),
              std_sd   = sd(value))
  
  return(list(dat = dat,
              mean_sd = mean_sd))
})

# == (1) Fit to all pot. vars (w & wo & re ES) =========================== ####
# (glm - w & wo ES) ------------------------------------------------------ ####
# fit model forms & tidy results
# (w/o ES as factor variable)
tidy_gt.1_woES <- list(data_spp.std, clim_vars) %>% 
  pmap(function(sp.x, vars.x) {
  # start by looking at things without ES (our intended random effect)
  sp.x <- sp.x$dat %>% select(-ES)
  
  # specify the formula we'll be fitting
  lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a))),
                 fn.inter_form(names(sp.x %>% select(-p_a)))),
               collapse = " + ")
  formula.x <- paste0("p_a ~ ", lp) %>% as.formula
  
  # fit glm & tidy coef summaries with broom
  glm(formula.x,
      data = sp.x,
      family = binomial(link = "logit")) %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>% 
    filter(p.value > 0.1) %>%      # look at the sig ones...
    rename(p.value_woES = p.value)
})

# (w ES as factor variable)
tidy_gt.1_wES <- list(data_spp.std, clim_vars) %>% 
  pmap(function(sp.x, vars.x) {
  # start by looking at things without ES (our intended random effect)
  # - that was rough! add it back in!!
  sp.x <- sp.x$dat
  
  # grab unique ES names for later...
  ES.x  <- paste0("ES", unique(sp.x$ES))
  
  # specify the formula we'll be fitting
  lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                 fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                 "ES"),
               collapse = " + ")
  formula.x <- paste0("p_a ~ ", lp) %>% as.formula
  
  # fit glm & tidy coef summaries with broom
  glm(formula.x,
      data = sp.x,
      family = binomial(link = "logit")) %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>% 
    filter(!term %in% ES.x) %>% # don't care if these are sig...
    filter(p.value > 0.1) %>% # look at the sig ones...
    rename(p.value_wES = p.value)
})

# (glmm - re ES) --------------------------------------------------------- ####
# these take a long time to fit, so save & load if revisiting this script:
if (file.exists(paste0(path_mod, "ch2-4-glmm_all_data_fits.rds"))) {
  tidy_gt.1_reES <- readRDS(paste0(path_mod, "ch2-4-glmm_all_data_fits.rds"))
} else {
  # fit model forms & tidy results
  tidy_gt.1_reES <- list(data_spp.std, clim_vars) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat

      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula

      # fit glmm
      # - setting nAGQ = 10 is based on work from a related project for these
      #   species (dissertation ch1 work...)
      glmm <- glmer(formula = formula.x,
                    data = sp.x,
                    family = binomial(link = "logit"),
                    nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
                    control = glmerControl(
                      optCtrl = list(maxfun = 1e8), # max iter for optimization
                      check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))

      # fit glm & tidy coef summaries with broom
      glmm %>%
      tidy %>%
      mutate_if(is.numeric, round, 3) %>%
      filter(p.value > 0.1) %>% # look at the sig ones...
      rename(p.value_reES = p.value)
    })

  # quickly save as rds as these take a loooong time to fit...
  saveRDS(tidy_gt.1_reES,
          file = paste0(path_mod, "ch2-4-glmm_all_data_fits.rds"))
  
}

# inspect potentially important covars ----------------------------------- ####
# look at them !!!
tidy_pvalues <- list(tidy_gt.1_woES, tidy_gt.1_wES, tidy_gt.1_reES) %>% 
  pmap(function(x, y, z) {
  full_join(x %>% select(term, p.value_woES),
            y %>% select(term, p.value_wES),
            by = "term") %>%
      full_join(z %>% select(term, p.value_reES),
                by = "term")
})

# == (2) re-fit based on new set of clim_vars ============================ ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio8", "bio10", "bio16"),
  "psmem"  = c("gdd5", "bio8", "bio10", "bio16", "bio18"),
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),
  "quke"   = c("hmi", "bio5", "bio8", "bio16"))

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
  sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
  sp.x
})

# (glmm - re ES) --------------------------------------------------------- ####
# Again... load if already fit as this takes a while...
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_2.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_2.rds"))
} else {
  # fit model forms & tidy results
  glm0 <- list(data_spp.std_tmp, clim_vars) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_2.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0 %>% map(function(x) {
  x %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(p.value > 0) %>% # look at the sig ones...
    rename(p.value_reES = p.value)
})

# == (3) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio8", "bio10", "bio16", "bio18"),          # updated
  "psmem"  = c("gdd5", "bio8", "bio16", "bio18"),           # updated
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),           # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),  # final
  "quke"   = c("bio13", "bio5", "bio8", "bio16"))             # updated

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_3.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_3.rds"))
} else {
  # fit model forms & tidy results
  glm0[c(1,2,5)] <- list(data_spp.std_tmp[c(1,2,5)], clim_vars[c(1,2,5)]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_3.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[c(1,2,5)] %>% map(function(x) {
  x %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(p.value > 0.1) %>% # look at the sig ones...
    rename(p.value_reES = p.value)
})
glm0[c(1,2,5)] %>% map(~.x %>% tidy)

# == (4) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                                # updated
  "psmem"  = c("gdd5", "bio8", "bio12", "bio18"),                # updated           
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),                # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),       # final
  "quke"   = c("bio5", "bio8", "bio16", "gdd5", "hmi", "bio13")) # updated

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_4.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_4.rds"))
} else {
  # fit model forms & tidy results
  glm0[c(1,2,5)] <- list(data_spp.std_tmp[c(1,2,5)], clim_vars[c(1,2,5)]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_4.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[c(1,2,5)] %>% map(function(x) {
  x %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(p.value > 0.1) %>% # look at the sig ones...
    rename(p.value_reES = p.value)
})
glm0[c(1,2,5)] %>% map(~.x %>% tidy)

# == (5) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                                # final
  "psmem"  = c("gdd5", "bio8", "bio18"),                         # updated            
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),                # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),       # final
  "quke"   = c("bio5", "bio8", "bio16", "gdd5", "bio12", "bio13")) # updated

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_5.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_5.rds"))
} else {
  # fit model forms & tidy results
  glm0[c(2,5)] <- list(data_spp.std_tmp[c(2,5)], clim_vars[c(2,5)]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_5.rds"))
  
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[c(2,5)] %>% map(function(x) {
  x %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(p.value > 0.1) %>% # look at the sig ones...
    rename(p.value_reES = p.value)
})
glm0[c(2,5)] %>% map(~.x %>% tidy)

# == (6) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                                # final
  "psmem"  = c("gdd5", "bio8", "bio10", "bio18"),                # updated             
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),                # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),       # final
  "quke"   = c("gdd5", "bio5", "bio8", "bio16", "bio13"))        # updated 

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_6.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_6.rds"))
  
} else {
  # fit model forms & tidy results
  glm0[c(2,5)] <- list(data_spp.std_tmp[c(2,5)], clim_vars[c(2,5)]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_6.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[c(2,5)] %>% map(function(x) {
  x %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(p.value > 0.1) %>% # look at the sig ones...
    rename(p.value_reES = p.value)
})
glm0[c(2,5)] %>% map(~.x %>% tidy)

# == (7) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                           # final
  "psmem"  = c("gdd5", "bio8", "bio10", "bio18"),           # final             
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),           # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),  # final
  "quke"   = c("gdd5", "bio8", "bio16", "bio13"))           # updated 

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_7.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_7.rds"))
} else {
  # fit model forms & tidy results
  glm0[5] <- list(data_spp.std_tmp[5], clim_vars[5]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_7.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[[5]] %>%
    tidy %>%
    mutate_if(is.numeric, round, 3) %>%
    filter(p.value > 0.1) %>% # look at the sig ones...
    rename(p.value_reES = p.value)

glm0[[5]] %>% tidy

# == (8) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                           # final
  "psmem"  = c("gdd5", "bio8", "bio10", "bio18"),           # final             
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),           # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),  # final
  "quke"   = c("gdd5", "bio5", "bio16", "bio13"))           # updated 

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_8.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_8.rds"))
} else {
  # fit model forms & tidy results
  glm0[5] <- list(data_spp.std_tmp[5], clim_vars[5]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_8.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[[5]] %>%
  tidy %>%
  mutate_if(is.numeric, round, 3) %>%
  filter(p.value > 0) %>% # look at the sig ones...
  rename(p.value_reES = p.value)

glm0[[5]] %>% tidy

# == (9) re-fit {abpr, psmem, quke} w. new set =========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                           # final
  "psmem"  = c("gdd5", "bio8", "bio10", "bio18"),           # final             
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),           # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),  # final
  "quke"   = c("gdd5", "bio13", "bio8", "bio10", "bio16", "bio18")) # updated 

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_9.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_9.rds"))
} else {
  # fit model forms & tidy results
  glm0[5] <- list(data_spp.std_tmp[5], clim_vars[5]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_9.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[[5]] %>%
  tidy %>%
  mutate_if(is.numeric, round, 3) %>%
  filter(p.value > 0) %>% # look at the sig ones...
  rename(p.value_reES = p.value)

glm0[[5]] %>% tidy

# == (10) re-fit {abpr, psmem, quke} w. new set ========================== ####
# subset to new vars ----------------------------------------------------- ####
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                           # final
  "psmem"  = c("gdd5", "bio8", "bio10", "bio18"),           # final             
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),           # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),  # final
  "quke"   = c("gdd5", "bio13", "bio8", "bio5")) # updated 

data_spp.std_tmp <- list(data_spp.std, clim_vars) %>%
  pmap(function(sp.x, clim.x) {
    sp.x$dat <- sp.x$dat %>% select("p_a", "ES", clim.x)
    sp.x
  })

# (glmm - re ES) --------------------------------------------------------- ####
if (file.exists(paste0(path_mod, "ch2-4-glmm_fits_10.rds"))) {
  # read back in...
  glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_10.rds"))
} else {
  # fit model forms & tidy results
  glm0[5] <- list(data_spp.std_tmp[5], clim_vars[5]) %>%
    pmap(function(sp.x, vars.x) {
      # start by looking at things without ES (our intended random effect)
      # - that was rough! add it back in!!
      sp.x <- sp.x$dat
      
      # specify the formula we'll be fitting
      lp <- paste0(c(fn.quad_form(names(sp.x %>% select(-p_a, -ES))),
                     fn.inter_form138(names(sp.x %>% select(-p_a, -ES))),
                     "(1 | ES)"),
                   collapse = " + ")
      formula.x <- paste0("p_a ~ ", lp) %>% as.formula
      
      # fit glmm (nAGQ = 10 based on ch1 work...)
      glmer(formula = formula.x,
            data = sp.x,
            family = binomial(link = "logit"),
            nAGQ = 10,  # no. interp pts for adaptive GHQ (approx)
            control = glmerControl(
              optCtrl = list(maxfun = 1e8), # max iter for optimization
              check.conv.grad = .makeCC("warning", tol = 2e-3, relTol = NULL)))
    })
  
  # quickly save as rds as these take a loooong time to fit...
  saveRDS(glm0,
          file = paste0(path_mod, "ch2-4-glmm_fits_10.rds"))
}

# inspect covars --------------------------------------------------------- ####
# fit glm & tidy coef summaries with broom
glm0[[5]] %>%
  tidy %>%
  mutate_if(is.numeric, round, 3) %>%
  filter(p.value > 0) %>% # look at the sig ones...
  rename(p.value_reES = p.value)

glm0[[5]] %>% tidy

# == Subset data based on final model forms ============================== ####
# final model vars 
clim_vars  <- list(
  "abpr"   = c("bio10", "bio18"),                           # final
  "psmem"  = c("gdd5", "bio8", "bio10", "bio18"),           # final             
  "qudo"   = c("bio6", "bio8", "bio16", "bio18"),           # final
  "quga4"  = c("bio5", "bio8", "bio10", "bio16", "bio18"),  # final
  "quke"   = c("gdd5", "bio13", "bio8", "bio5"))            # final

other_vars <- c(grep("bio", names(FIA), value = TRUE), "hmi", "gdd5")
other_vars <- names(FIA %>% select(-all_of(c(other_vars, spp))))

data_spp <- list(FIA_spp, clim_vars) %>% pmap(function(sp.x, clim.x) {
  sp.x %>% select(all_of(c(other_vars, "p_a", clim.x)))
})
data_spp.std <- list(data_spp, clim_vars) %>% pmap(function(sp.x, clim.x) {
  dat <- sp.x %>% 
    mutate_at(clim.x, function(x) (x - mean(x)) / sd(x))
  
  mean_sd <- sp.x %>% 
    select(all_of(clim.x)) %>% 
    melt(variable.name = "clim_var") %>%
    group_by(clim_var) %>%
    summarize(std_mean = mean(value),
              std_sd   = sd(value))
  
  return(list(dat = dat,
              mean_sd = mean_sd))
})

# grab our 'starting values' from this last naive fit
glm0 <- readRDS(paste0(path_mod, "ch2-4-glmm_fits_10.rds"))
parms_start_spp <- glm0 %>% map(function(glm.x) {
  beta_est <- fixef(glm.x)
  sigma_v_est <- sqrt(VarCorr(glm.x)[[1]][1])
  c(beta_est, sigma_v_est) %>% setNames(c(names(beta_est), "sigma_v"))
})

# == save subsetted and std data ========================================= ####
save(clim_vars, data_spp, data_spp.std, parms_start_spp,
     file = paste0(path_data, "ch2-4-modeling-data.Rdata"))
