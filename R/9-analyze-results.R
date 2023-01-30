###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 9-analyze-results.R
## Updated - 01-28-2023
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Spatial-Bayesian models project shifts in 
## suitable habitat for Pacific Northwest tree species under climate change' by
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - Calculate summary statistics based on prediction maps. Figures 
##   correspond to:
##   - Manuscript: Figures 1-3 and Tables 1, 3*, & 4
##     (*the Moran's I model-check summaries for Table 3 are calculated in
##     script 12-...R)
##   - Appendix S1: Figures S7-S13
##
## About - output:
## - ch2-9-pred_covar_lims.rds
##   - location: path_output
##   - saved for use in this script; e.g., max/min ('limits') of covariates  
##     under the current & four future scenarios for each spp, where covariates  
##     are those used in modeling
###############################################################################
# data manipulation
library(magrittr)    # for set_rownames()
library(dplyr)       # for ... E V E R Y T H I N G ...
library(purrr)       # working with lists
library(reshape2)    # melt() and dcast()
library(stringr)     # for string help (e.g. str_detect in functions)

# plotting
library(ggplot2)
library(cowplot)     # for plot_grid, used in sp-clim rel. plots
library(patchwork)   # used for more complicated plot arrangements
library(sf)          # used for map-figures

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_output <- paste0(path_top, "output/")

# path to functions script
path_fn <- paste0(path_top, "scripts - real data/")

# where to ggsave-ed figures go
path_images <- paste0(path_top, "images/")

# == Load data & fns ===================================================== ####
# - fns & drop plot ids -------------------------------------------------- ####
# load fns
source(paste0(path_fn, "0-functions.R"))

# plots dropped for prediction
# - combination of plots outside of ClimateNA's extent (see script 7 comments)
#   and those too near the coast line (see script 8)
# - while we aren't predicting at these locations, they were used in 'fitting'
#   the SDMs
# - all drop plots are non-forested (e.g. unobserved response), so they only 
#   occurred in the qudo and quga4 data sets, but obj is a list of 5 by spp
drop_plot_ids <- readRDS(file = paste0(path_output, "ch2-8-drop_plot_ids.rds"))

# - variable std & training data ----------------------------------------- ####
# variable std (for beta plots)
load(paste0(path_data, "ch2-4-modeling-data.Rdata")) # subset modeling data
rm(data_spp, parms_start_spp)

# current - training data (formatted & standardized)
spp    <- names(data_spp.std)  # extract spp names
dat_ls <- lapply(spp, function(x) {
  readRDS(paste0(path_output, "ch2-5-sampler-ready-data_", x, ".rds"))
}) %>% setNames(spp)

# - final model run ------------------------------------------------------ ####
# model: final MCMC samples/run
mod_spp <- readRDS(paste0(path_output, "ch2-6-final-MCMC-samples.rds")) 

# - create|load lims of covars for preds (current & future) -------------- ####
# About: 
# - for this script we only need the limits of covars under current & future 
#   scenarios
# - all 'newdat.x' from current & future are std according to same scale

# load or make object
if (file.exists(paste0(path_output, "ch2-9-pred_covar_lims.rds"))) {
  pred_covar_lims <- readRDS(paste0(path_output, "ch2-9-pred_covar_lims.rds"))
} else {
  # newdat.x from current climate predictions
  PRISM_dat <- readRDS(paste0(path_output, "ch2-7-preds-PRISM.rds"))$preds_PRISM %>%
    map(function(sp.x) {
      keep_vars <- names(sp.x$newdat.x)
      keep_vars <- c(grep("bio", keep_vars, value = TRUE),
                     grep("gdd5", keep_vars, value = TRUE),
                     grep("hmi", keep_vars, value = TRUE))
      keep_vars <- keep_vars[!grepl("_", keep_vars)]
      
      sp.x$newdat.x %>%
        select(one_of(keep_vars)) %>% # really all_of, but tour can't handle...
        melt %>%
        group_by(variable) %>%
        summarize(max_value = max(value),
                  min_value = min(value),
                  q_10 = quantile(value, .10) %>% as.numeric,
                  q_25 = quantile(value, .25) %>% as.numeric,
                  q_75 = quantile(value, .75) %>% as.numeric,
                  q_90 = quantile(value, .90) %>% as.numeric)
    })
  
  # newdat.x from future climate predictions
  # - names of future scenarios
  future_scenarios <- readRDS(paste0(path_output, "ch2-7-future_scenarios_chr.rds"))
  # - read-in data
  CNA_dat <- lapply(future_scenarios, function(x) {
    readRDS(paste0(path_output, "ch2-7-preds-CNA_", x, ".rds")) %>%
      map(function(sp.x) {
        sp.x$newdat.x %>% 
          select(-plot_id) %>%
          melt %>% 
          group_by(variable) %>%
          summarize(max_value = max(value),
                    min_value = min(value),
                    q_10 = quantile(value, .10) %>% as.numeric,
                    q_25 = quantile(value, .25) %>% as.numeric,
                    q_75 = quantile(value, .75) %>% as.numeric,
                    q_90 = quantile(value, .90) %>% as.numeric)
      }) 
  }) %>% setNames(future_scenarios)
  
  # merge the two lists of newdat.x summaries together by spp
  pred_covar_lims <- append(CNA_dat, list(PRISM_dat)) %>%
    setNames(c(future_scenarios, "current")) %>% 
    transpose %>%
    map(function(scenario.x) scenario.x %>% bind_rows(.id = "scenario"))
  
  # clean-up space
  rm(PRISM_dat, CNA_dat)
  
  # only keep base-vars (e.g. drop quad & inter terms...)
  pred_covar_lims <- pred_covar_lims %>%
    map(function(sp.x) {
      keep_vars <- sp.x$variable %>% unique
      keep_vars <- keep_vars[!grepl("_2", keep_vars)] # drop quad terms
      keep_vars <- keep_vars[!grepl("_b", keep_vars)] # drop inter terms
      
      sp.x %>% filter(variable %in% keep_vars)
    })
  
  # add column with non-std values
  # - here moving bio16 & bio18 to totals (not monthly)
  pred_covar_lims  <- list(pred_covar_lims,
                           data_spp.std %>% map(~.x$mean_sd)) %>%
    pmap(function(sp.x, std.x) {
      sp.x <- sp.x %>%
        inner_join(std.x %>% rename(variable = clim_var), 
                   by = "variable") %>%
        # add column with un-standardized values
        # - here need sapply b/c looking at whole df...
        mutate_if(sapply(., is.numeric) & !grepl("std", names(.)), 
                  list(nonstd = ~(.x * std_sd) + std_mean)) %>%
        # make bio16 & bio18 vars into totals (if they exist)
        mutate_at(vars(contains("nonstd")), 
                  ~ifelse(variable %in% paste0("bio", 16:19), .x * 3, .x)) %>%
        select(-std_mean, -std_sd) 
    })
  
  # collapse down to neater format
  pred_covar_lims <- pred_covar_lims %>%
    map(function(sp.x) {
      sp.x %>%
        melt(id.vars = c("scenario", "variable"),
             variable.name = "stat_type") %>%
        mutate(std_type = ifelse(grepl("nonstd", stat_type), "nonstd", "std")) %>%
        mutate(stat_type = gsub("_nonstd", "", stat_type)) %>%
        rename(clim_var = variable) %>%
        select(scenario, clim_var, stat_type, std_type, value)
    })
  
  # save this file -- process takes a while, so do this to speed things up...
  saveRDS(pred_covar_lims,
          file = paste0(path_output, "ch2-9-pred_covar_lims.rds"))
}

# - create crosswalks between names for covars --------------------------- ####
# covariate long names and units
# - make crosswalk with longer variable names & units for figures & tables 
long_name_var <- data.frame(
  orig_names = c("gdd5", "bio8", "bio6", "bio5", 
                 "bio18", "bio16", "bio13", "bio10",
                 "intercept"),
  long_names = c("growing degree days",
                 "mean temp of wettest qtr",
                 "min temp of coldest month",
                 "max temp of warmest month",
                 "ppt in warmest qtr",
                 "ppt in wettest qtr",
                 "ppt in wettest month",
                 "mean temp of warmest qtr",
                 "FE intercept"),
  # add new short_names for Ecosphere ('EJ') figures
  EJ_names = c("gdd5",
               "tmean-wetQ",
               "tmin-coldM",
               "tmax-warmM",
               "ppt-warmQ",
               "ppt-wetQ",
               "ppt-wetM",
               "tmean-warmQ",
               "intercept"),
  units = c("days",
            rep("C", 3),
            rep("mm", 3), # all ppt vars (PRISM and ClimateNA) are in mm
            "C",
            ""))  %>%
  mutate(both_names = if_else(EJ_names == "intercept",
                              EJ_names,
                              paste0(EJ_names, " (", units, ")")))

# betas to covariate short names
# - crosswalk original beta names to model's "beta_" names for plotting
beta_crosswalk <- list(mod_spp %>% map(~.x$MCMC_theta), 
                       dat_ls %>% map(~.x$parms_names_orig)) %>%
  pmap(function(MCMC_theta.x, parms_orig.x) {
    p <- MCMC_theta.x %>% names %>% grepl("beta", .) %>% sum
    
    data.frame(mod_name = names(MCMC_theta.x)[1:p],
               orig_names = parms_orig.x[1:p])
  }) %>%
  # add names to use in Ecosphere figures
  map(function(df.x) {
    df.x %>% mutate(EJ_names = str_replace_all(
      orig_names,
      with(long_name_var, setNames(EJ_names, orig_names))))
  })

# - spatial data for mapping (z's and ES's) ------------------------------ ####
# study area boundary mapping data
study_area.sf <- readRDS(paste0(path_output, 
                                "ch2-8-study_area_sf.rds"))
# spatial hex mapping data
hex_polygons_public.sf <- readRDS(paste0(path_output, 
                                         "ch2-8-hex_polygons_public_sf.rds"))
# spatial plot mapping data
plot_polygons_public.sf <- readRDS(paste0(path_output, 
                                          "ch2-8-plot_polygons_public_sf.rds"))

# == Arguments for figures =============================================== ####
# common names for species
spp_common <- c("noble fir",
                "coastal Douglas-fir",
                "blue oak",
                "white oak",
                "black oak") %>%
  setNames(spp)

# set general text-size & figure width arguments (based on Ecosphere's author
# guidelines)
# - text: 6-10 pt
EJ_axis_text_size  <- 7
EJ_strip_text_size <- 7
EJ_axis_title_size <- 8
EJ_label_size      <- 10

# - max dimensions
EJ_dim_units <- "cm"
EJ_dim_1col_width_max <- 8.5 #cm
EJ_dim_2col_width_max <- 18  #cm
EJ_dim_2col_ht_max    <- 22  #cm

# - max dim appendices (b/c supplying a pdf of word doc with 1" margins)
EJ_app_dim_units <- "in"
EJ_app_dim_width_max <- 6.5 #in
EJ_app_dim_ht_max    <- 9   #in

# - dpi: min 300-600
EJ_dpi <- 800

# == Table 1 & hex-summaries in text ===================================== ####
# Need to report the mean & maximum number of plots per hexagon...
# (by spp)
dat_ls %>% 
  map(function(df.x) {
    df.x$moddat %>% 
      distinct(fiahex_id, .keep_all = TRUE) %>% 
      # (T/F were plots spatially intensified in this hexagon?)
      mutate(hex_intense = n_plots_hex > 1) %>% 
      # filter(hex_intense == TRUE) %>% 
      # (average number of plots/hex by hex-type)
      summarize(avg_n = mean(n_plots_hex) %>% round(2),
                med_n = median(n_plots_hex) %>% round(2),
                max_n = max(n_plots_hex) %>% round(2))
  }) %>%
  bind_rows(.id = "spp") 

# (across study area ... this is in the paper)
dat_ls %>%
  map(~.x$moddat %>% select(fiahex_id, n_plots_hex)) %>%
  bind_rows %>%
  group_by(fiahex_id) %>%
  filter(n_plots_hex == max(n_plots_hex)) %>%
  ungroup %>%
  distinct(fiahex_id, .keep_all = TRUE) %>% 
  # (T/F were plots spatially intensified in this hexagon?)
  mutate(hex_intense = n_plots_hex > 1) %>% 
  # filter(hex_intense == TRUE) %>% 
  # (average number of plots/hex by hex-type)
  summarize(avg_n = mean(n_plots_hex) %>% round(2),
            med_n = median(n_plots_hex) %>% round(2),
            max_n = max(n_plots_hex) %>% round(2))

# Table 1 summary stats...
dat_ls %>%
  map(function(df.x) {
    any_missing <- (df.x$moddat %>% filter(y_index == 0) %>% nrow) > 0
    
    prevalence <- ((sum((df.x$moddat %>% filter(y_index == 1))$p_a) / 
                     nrow(df.x$moddat %>% filter(y_index == 1))) * 100) %>%
      round(1)
    
    n_plots <- df.x$moddat %>% 
      group_by(y_index) %>% 
      summarize(n_plots = n())
   
    n_hex   <- df.x$moddat %>% 
      group_by(fiahex_id) %>% 
      summarize(z_type = as.numeric(sum(y_index) > 0)) %>% 
      ungroup %>% 
      group_by(z_type) %>% 
      summarize(n_hex = n()) 
    
    n_ES <- n_distinct(df.x$moddat$ES)
    
    if (any_missing) {
      data.frame(
        prevalence = prevalence,
        y_o = n_plots %>% filter(y_index == 1) %>% pull(n_plots),
        y_u = n_plots %>% filter(y_index == 0) %>% pull(n_plots),
        z_o = n_hex %>% filter(z_type == 1) %>% pull(n_hex),
        z_u = n_hex %>% filter(z_type == 0) %>% pull(n_hex),
        ES  = n_ES)
    } else {
      data.frame(
        prevalence = prevalence,
        y_o = n_plots %>% filter(y_index == 1) %>% pull(n_plots),
        y_u = NA,
        z_o = n_hex %>% filter(z_type == 1) %>% pull(n_hex),
        z_u = NA,
        ES  = n_ES)
    }
    
  }) %>%
  bind_rows(.id = "spp")
  

# == Appendix S1: posterior summaries for parameters (not z) ============= ####
# - Figure S8: ES (v) ---------------------------------------------------- ####
p_body <- (mod_spp %>% 
             map(function(mod.x) {
               mod.x$MCMC_v %>%
                 melt(variable.name = "ES") %>%
                 mutate(ES = gsub("ES", "", ES)) %>%
                 group_by(ES) %>%
                 summarize(y_min = min(value),
                           y_max = max(value),
                           ql = quantile(value, .05) %>% unname,
                           qu = quantile(value, .95) %>% unname,
                           y_mean = mean(value))
             }) %>% 
             bind_rows(.id = "spp") %>%
             mutate(spp = recode(spp, !!!spp_common),
                    spp = factor(spp, levels = unname(spp_common)))) %>%
  ggplot() +
  geom_boxplot(aes(x = ES,
                   ymin = y_min,
                   lower = ql,
                   middle = y_mean,
                   upper = qu,
                   ymax = y_max),
               position = position_dodge(width = 0),
               stat = "identity",
               width = .5,
               alpha = .7,
               fill = "grey") +
  geom_hline(yintercept = 0, color = "blue") +
  ylab("value") +
  theme_bw() +
  theme(strip.text.x = element_text(size = EJ_strip_text_size),
        axis.title = element_text(size = EJ_axis_title_size),
        axis.text  = element_text(size = EJ_axis_text_size,
                                  colour="black"),
        legend.position = "none") +
  facet_wrap(~ spp, nrow = 1, scales = "free_x") +
  coord_flip() 

# view it
p_body +
  ggtitle("Posterior densities for v's (random intercept for ES) by species",
          subtitle = "based on 2k MCMC samples")

# save it (w/o title, and as compressed tiff)
ggsave(filename = paste0(path_images, 
                         "FigS8.tiff"),
       plot = p_body,
       height = EJ_app_dim_width_max, 
       width = EJ_app_dim_width_max,
       units = EJ_app_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")

# - Figure S7: betas, sigma_z, sigma_v, rho ------------------------------ ####
p_body <- list(mod_spp %>% map(~.x$MCMC_theta), beta_crosswalk) %>%
  pmap(function(theta, cross) {
    # get bio-style names (e.g. bio10_bio18, etc)
    names(theta)[match(cross$mod_name, names(theta))] <- cross$EJ_names
    
    # set order for FE parameter factor levels
    new_levels <- c("bio5", "bio6", "bio8", 
                    "bio10", "bio13", "bio16", "bio18")
    new_levels <- c("intercept", "gdd5",
                    new_levels,
                    paste0(new_levels, "_2"),
                    "bio8_bio13", "bio8_bio16", "bio10_bio18")
    # - now switch to Ecosphere-styled names
    new_levels <- new_levels %>%
      str_replace_all(with(long_name_var, setNames(EJ_names, orig_names)))
    
    # summarize into boxplot-mashable format...
    theta <- theta %>%
      select(-n_sim) %>%
      melt(variable.name = "theta") %>%
      group_by(theta) %>%
      summarize(y_min = min(value),
                y_max = max(value),
                ql = quantile(value, .05) %>% unname,
                qu = quantile(value, .95) %>% unname,
                y_mean = mean(value))
    
    # split into FE parameters & others
    list(theta %>% 
           filter(!theta %in% c("rho", "sigma_z", "sigma_v")) %>%
           mutate(theta = factor(theta, levels = new_levels)),
         theta %>% 
           filter(theta %in% c("rho", "sigma_z", "sigma_v")) %>%
           mutate(theta = factor(theta, levels = c("rho", "sigma_z", "sigma_v"))))
  }) %>%
  transpose() %>%
  map(function(theta.x) {
    (theta.x %>% 
       bind_rows(.id = "spp") %>%
       mutate(spp = recode(spp, !!!spp_common),
              spp = factor(spp, levels = rev(unname(spp_common))))) %>%
      ggplot() +
      geom_boxplot(aes(x = spp,
                       ymin = y_min,
                       lower = ql,
                       middle = y_mean,
                       upper = qu,
                       ymax = y_max),
                   stat = "identity",
                   width = .5,
                   alpha = .7,
                   fill = "grey") +
      theme_bw() +
      labs(y = "value") +
      theme(strip.text.x = element_text(size = EJ_strip_text_size),
            axis.title.x = element_text(size = EJ_axis_title_size),
            axis.text  = element_text(size = EJ_axis_text_size,
                                      colour="black"),
            legend.position = "none", 
            axis.title.y = element_blank()) +
      facet_wrap(~ theta, ncol = 4, scales = "free_x") +
      coord_flip() 
  }) 

((p_body[[2]] + theme(axis.title.x = element_blank())) + 
    plot_spacer() + 
    plot_layout(widths = c(3.25, 1))) / 
  (p_body[[1]] + geom_hline(yintercept = 0, color = "blue")) + 
  plot_layout(heights = c(1, 7)) +
  plot_annotation(title = "Posterior summaries for parameters (excluding v or z) by species",
                  subtitle = "Based on 2k MCMC samples")

# save it (w/o title, and as compressed tiff)
ggsave(filename = paste0(path_images, 
                         "FigS7.tiff"),
       plot = ((p_body[[2]] + theme(axis.title.x = element_blank())) + 
                 plot_spacer() + 
                 plot_layout(widths = c(3.25, 1))) / 
         (p_body[[1]] + geom_hline(yintercept = 0, color = "blue")) + 
         plot_layout(heights = c(1, 7)),
       height = EJ_app_dim_width_max, 
       width = EJ_app_dim_width_max,
       units = EJ_app_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")

# what are the posterior means?
list(mod_spp, beta_crosswalk) %>% 
  pmap(function(mod.x, crosswalk.x) {
    post_mean <- mod.x$MCMC_theta %>% apply(2, mean)
    post_mean <- post_mean[1:(length(post_mean) - 4)] # don't need rho/sigma's/n_sim
    
    # get bio* style names...
    update_index <- match(crosswalk.x$mod_name, names(post_mean))
    names(post_mean)[update_index] <- crosswalk.x$EJ_names
    
    return(post_mean)
  }) %>%
  map(~.x %>% round(2))

# == sp-climate relationship ============================================= ####
# - conditional sp-climate relationships --------------------------------- ####
cond_plots <- list(mod_spp %>% map(~.x$MCMC_theta),
                   dat_ls %>% map(~.x[c("sp.x", "parms_names_orig")]),
                   data_spp.std %>% map(~.x$mean_sd),
                   pred_covar_lims,
                   beta_crosswalk) %>%
  pmap(~fn.cond_spclim_rel(..1, ..2, ..3, ..4, ..5,
                           grid_size = 100,
                           long_name_var,
                           long_names = T))

# - Summaries for manuscript text (Table 4) ------------------------------ ####
# mean peak & optimal mode (Table 4 info from here)
# - ignore the mean_peak -- not used
cond_plots %>% map(function(sp.x) {
  # summarize by relationship...
  sp.x %>% imap(function(clim.x, clim_name.x) {
    # handle different types of relationships differently
    if(grepl("_b", clim_name.x)) {
      # bivariate relationships (e.g. those with an interaction term)
      
      # - mean climate conditions for peak conditional probability
      x <- clim.x$mean_cond_peak %>% 
        select(starts_with("bio")) %>% 
        round(1)
      mean_peak <- paste0("(", x[1], ", ", x[2], ")")
      
      # - mode of optimal climate conditions for peak cond. prob.
      x <- clim.x$mode %>% select(-density) %>% round(1)
      mode_opt <- paste0("(", x[1], ", ", x[2], ")")
      
      data.frame(mean_peak = mean_peak,
                 n_opt_max_mcmc = nrow(clim.x$opt_cond_mcmc %>% filter(opt_type == "max")),
                 n_opt_maxrst_mcmc = nrow(clim.x$opt_cond_restricted),
                 mode_opt_maxrst = mode_opt,
                 opt_maxrst_90CI = NA)
      
    } else if(grepl("gdd5", clim_name.x)) {
      # univariate relationships with only linear terms
      # - easy here b/c only gdd5 in this category
      data.frame(mean_peak = clim.x$mean_cond_peak %>% 
                   select(clim_name.x) %>% 
                   round(1) %>% 
                   as.character,
                 n_opt_max_mcmc = NA,
                 n_opt_maxrst_mcmc = NA,
                 mode_opt_maxrst = NA,
                 opt_maxrst_90CI = NA)
    } else {
      # univariate relationships with lin & quad terms
      # - all the other non interaction bio* covariates
      data.frame(
        mean_peak = clim.x$mean_cond_peak %>% 
          select(clim_name.x) %>% 
          round(1) %>% 
          as.character,
        n_opt_max_mcmc = nrow(clim.x$opt_cond_mcmc %>% filter(opt_type == "max")),
        n_opt_maxrst_mcmc = nrow(clim.x$opt_cond_restricted),
        mode_opt_maxrst = clim.x$mode %>% 
          round(1) %>% 
          as.character,
        opt_maxrst_90CI = clim.x$opt_cond_restricted %>% 
          select(clim_name.x) %>%
          set_colnames("dummy") %>%
          summarize(ql = quantile(dummy, .05) %>% as.numeric %>% round(1),
                    qu = quantile(dummy, .95) %>% as.numeric %>% round(1)) %>%
          with(., paste0("(", ql, ", ", qu, ")"))
        )
    }
  }) %>%
    bind_rows(.id = "relationship") %>%
    mutate(relationship = 
             factor(relationship,
                    levels = c("gdd5", "bio18", "bio5", "bio6", "bio8",  
                               "bio8_bio13", "bio8_bio16", "bio10_bio18"))) %>%
    arrange(relationship)
}) %>% bind_rows(.id = "spp")

# - Figures 3 & S9-S13: Plots for manuscript & appendix ------------------ ####
# Figures 3 & Appendix S1 Figures S9-S12: mean peak & optimal condition by spp
# - make plots for manuscript
p_body <- cond_plots %>% imap(function(ls.x, sp_name.x) {
  # - pull out the plots we need to arrange ----------------------- ####
  top.x <- ls.x %>% 
    map(~.x$p1_meanpred + 
          theme(axis.title = element_text(size = EJ_axis_title_size),
                axis.text  = element_text(size = EJ_axis_text_size,
                                          colour="black"),
                legend.title = element_text(size = EJ_axis_title_size),
                legend.text  = element_text(size = EJ_strip_text_size),
                legend.key.height = unit(4, "mm"),
                legend.key.width = unit(3, "mm"),
                legend.box.spacing = unit(0.1, "mm"),
                plot.margin = unit(c(4, 0, 1, 0), "mm")))
  
  opt.x <- ls.x %>% 
    map(~.x$p1_opt + 
          theme(axis.title = element_text(size = EJ_axis_title_size),
                axis.text  = element_text(size = EJ_axis_text_size,
                                          colour="black"),
                legend.title = element_text(size = EJ_axis_title_size),
                legend.text  = element_text(size = EJ_strip_text_size),
                legend.key.height = unit(4, "mm"),
                legend.key.width = unit(3, "mm"),
                legend.box.spacing = unit(0.1, "mm"),
                plot.margin = unit(c(4, 0, 1, 0), "mm")))
  
  # - set rel_widths (to handle plots with legends) --------------- ####
  n_inter <- grepl("_b", names(ls.x)) %>% sum
  n_other <- length(ls.x) - n_inter
  rel_widths <- c(rep(500, n_other),
                  rep(600, n_inter))
  rel_widths <- rel_widths / sum(rel_widths)
  
  # - make the figure! -------------------------------------------- ####
  # plot differently for abpr since there is only one plot
  # (all others plot vertically, although small adjustment for 
  #  psmem & quke inside of the code below b/c of gdd5's lack of opt)
  if (length(top.x) == 1) {
    plot_body <- plot_grid(top.x[[1]], opt.x[[1]],
                           align = "h",
                           axis = "ltb",
                           nrow = 1,
                           labels = c("(a)", "(b)"),
                           label_size = EJ_label_size, 
                           hjust = 0, vjust = 1)
  } else {
    # align by type horizontally, then stack
    p_top <- plot_grid(
      plotlist = rev(top.x),
      align = "h",
      axis = "lb",
      nrow = 1,
      rel_widths = rel_widths,
      labels = paste0("(", letters[1:length(rel_widths)], ")"),
      label_size = EJ_label_size, 
      hjust = 0, 
      vjust = 1)
    
    # - adjust for any spp with gdd5 as a covariate (psmem & quke)
    if (any(opt.x %>% map(~is.null(.x)) %>% unlist)) {
      # this is messy to do dynamically, so I've hardcoded it...
      # (luckily these two species only have 3 relationships to consider)
      p_bottom <- plot_grid(opt.x[[2]], opt.x[[1]],
                            align = "h",
                            axis = "lb",
                            nrow = 1,
                            rel_widths = rel_widths[2:3],
                            labels = paste0("(", 
                                            letters[(length(rel_widths) + 1):(length(rel_widths)*2-1)], # changed this....
                                            ")"), 
                            label_size = EJ_label_size, 
                            hjust = 0, 
                            vjust = 1)
      p_bottom <- plot_grid(NULL, p_bottom,
                            rel_widths = c(rel_widths[1], sum(rel_widths[2:3])))
      
      plot_body <- plot_grid(p_top, p_bottom,
                             align = "v",
                             axis = "l",
                             nrow = 2,
                             rel_heights = c(1, 1))
      
    } else {
      p_bottom <- plot_grid(
        plotlist = rev(opt.x),
        align = "h",
        axis = "lb",
        nrow = 1,
        rel_widths = rel_widths,
        labels = paste0("(", 
                        letters[(length(rel_widths) + 1):(length(rel_widths)*2)], 
                        ")"),
        label_size = EJ_label_size, 
        hjust = 0, 
        vjust = 1)
      
      plot_body <- plot_grid(p_top, p_bottom,
                             align = "v",
                             axis = "l",
                             nrow = 2,
                             rel_heights = c(1, 1))
    }
    
  }
  
  return(plot_body)
})
p_body %>% # saving by species
  list(.,
       # Ecosphere sizing: keeping trial-and-error best-sizes for appendix
       # figures & adjusting psmem Fig to manu's max-width
       c(2.03, 
         3.35*(EJ_dim_2col_width_max/(6*2.54)),
         rep(3.35*(6.5/6), 3)),
       c(5.5, 
         EJ_dim_2col_width_max/2.54, 
         rep(6.5, 3)),
       c("S9", "3", paste0("S", 10:12))) %>%
  pmap(function(plot.x, ht.x, width.x, name.x) {
    ggsave(filename = paste0(path_images, 
                             "Fig",
                             name.x,
                             ".tiff"),
           plot = plot.x,
           height = ht.x, 
           width = width.x,
           units = "in", 
           dpi = EJ_dpi,
           compression = "lzw")
  })

# Section S5 Figure S13: Variability in mean peak plots for interaction terms
# - legend would be: Standard Error of the conditional probability (based
#   on 2,000 MCMC samples)
se_plots <- with(cond_plots %>% map(function(sp.x) {
  sp.x[grepl("_b", names(sp.x))] %>% 
    map(~.x$p1_se + 
          theme(axis.title = element_text(size = EJ_axis_title_size),
                axis.text  = element_text(size = EJ_axis_text_size,
                                          colour="black"),
                legend.title = element_text(size = EJ_axis_title_size),
                legend.text  = element_text(size = EJ_strip_text_size),
                legend.key.height = unit(4, "mm"),
                legend.key.width = unit(3, "mm"),
                legend.box.spacing = unit(0.1, "mm"),
                plot.margin = margin(5,0,0,2)) + 
          guides(fill = guide_colorbar(title = "SE")))
}),
plot_grid(abpr[[1]], psmem[[1]], qudo[[1]],
          quga4[[1]], quga4[[2]], quke[[1]],
          nrow = 2, align = "hv", axis = "l",
          labels = paste0("(", letters[1:6], ")"),
          label_size = EJ_label_size,
          label_x = 0, label_y = 1,
          hjust = -0.2, vjust = 1.8))
ggsave(filename = paste0(path_images, "FigS13.tiff"),
       plot = se_plots,
       height = 3.26*(6.5/6), 
       width = 6.5,
       units = "in", 
       dpi = EJ_dpi,
       compression = "lzw")


# - (%) response curve shape --------------------------------------------- ####
mod_spp %>% 
  map(~.x$MCMC_theta %>% select(starts_with("beta_"))) %>%
  list(., 
       beta_crosswalk,
       list(abpr = list(c("bio10", "bio18")),
            psmem = list("gdd5", "bio8", c("bio10", "bio18")),
            qudo = list("bio18", "bio6", c("bio8", "bio16")),
            qudo = list("bio5", c("bio8", "bio16"), c("bio10", "bio18")),
            qudo = list("gdd5", "bio5", c("bio8", "bio13")))) %>%
  pmap(function(theta.x, cw.x, rel.x) {
    # cols of theta.x in same order as rows of cw.x
    names(theta.x) <- cw.x$orig_names
    
    rel.x %>% 
      lapply(function(rel.i) {
        theta.i <- theta.x %>% 
          select(contains(rel.i)) %>%  # only grab vars in this rel
          select(any_of(names(theta.x))) # make sure order is same as it was originally
        
        if ("gdd5" %in% rel.i) {
          # only has L terms
          c("Sd" = mean(theta.i[,1] < 0),
            "Si" = mean(theta.i[,1] > 0),
            "F" = mean(theta.i[,1] == 0))
        } else if (length(rel.i) == 1) {
          # bio vars with L&Q terms, but not involved in an interaction
          c("Ud" = mean(theta.i[,2] < 0),
            "Ui" = mean(theta.i[,2] > 0)) %>%
            c(., "other" = (1-sum(.)))
        } else {
          # bio vars involved through an interaction (each with L&Q terms)
          c("Udd" = mean(theta.i[,3] < 0 & theta.i[,4] < 0),
            "Udi" = mean((theta.i[,3] > 0 & theta.i[,4] < 0) |
                             (theta.i[,3] < 0 & theta.i[,4] > 0))) %>%
            c(., "other" = (1-sum(.)))
        }
      }) %>% 
      setNames(rel.x %>% map_chr(~paste(.x, collapse = "_"))) %>%
      map(~.x[.x > 0] %>% bind_rows %>% melt) %>%
      bind_rows(.id = "clim_var") %>%
      mutate(value = round(value, 3)*100) %>%
      rename(resp_curve_shape = variable,
             percent_MCMC = value)
  })

# plot pairs of bivariate rel quad-terms by MCMC sample
# - ideally density will be in the lower-left quadrant, but the closer
#   to zero a parameter's sample was the steeper the slope of the curve
with((mod_spp %>% 
        map(~.x$MCMC_theta %>% select(starts_with("beta_"))) %>%
        list(., 
             beta_crosswalk,
             list(abpr = list(c("bio10", "bio18")),
                  psmem = list(c("bio10", "bio18")),
                  qudo = list(c("bio8", "bio16")),
                  qudo = list(c("bio8", "bio16"), c("bio10", "bio18")),
                  qudo = list(c("bio8", "bio13")))) %>%
        pmap(function(theta.x, cw.x, rel.x) {
          # cols of theta.x in same order as rows of cw.x
          names(theta.x) <- cw.x$orig_names
          
          rel.x %>% 
            lapply(function(rel.i) {
              var_names <- paste0(rel.i, "_2")
              (theta.x %>% 
                 select(contains("_2")) %>%
                 select(contains(rel.i))) %>%
                ggplot() +
                geom_density_2d(aes_string(x = var_names[1], y = var_names[2]),
                                color = "black",
                                bins = 10) +
                geom_hline(aes(yintercept = 0), color = "red") +
                geom_vline(aes(xintercept = 0), color = "red") +
                scale_fill_grey() +
                theme_bw() + 
                theme(aspect.ratio = 1)
            })
        })),
     plot_grid(abpr[[1]], psmem[[1]], qudo[[1]],
               quga4[[1]], quga4[[2]], quke[[1]],
               nrow = 2, align = "hv", axis = "l",
               labels = c("ABPR", "PSMEM", "QUDO", 
                          "QUGA4 (a)", "QUGA4 (b)", "QUKE"),
               label_fontface = "bold", label_size = 12))
# range of climate at plots where spp was observed (present)
data_spp.std %>% map(function(sp.x) {
    sp.x$dat %>% 
      filter(p_a == 1) %>% 
      select(contains("bio"), plot_id) %>%
      melt(id.vars = "plot_id",
           variable.name = "clim_var") %>%
      inner_join(sp.x$mean_sd, by = "clim_var") %>%
      mutate(value = value * std_sd + std_mean) %>%
      group_by(clim_var) %>%
      summarize(min = min(value),
                q_05 = quantile(value, 0.05) %>% as.numeric,
                q_50 = quantile(value, 0.5) %>% as.numeric,
                q_95 = quantile(value, 0.95) %>% as.numeric,
                max = max(value))
  })

# - Discussion: observed range climate covars ---------------------------- ####
data_spp.std %>% map(function(sp.x) {
  # grab names of covars..
  covars <- sp.x$mean_sd$clim_var %>% as.character
  
  # un-standardize and summarize by p_a grouping
  sp.x$dat %>% 
    filter(plot_obsv == 1,
           p_a == 1) %>%  # this last just kept to simplify, remove if needed
    select(all_of(covars), p_a, plot_id) %>%
    melt(id.vars = c("plot_id", "p_a"),
         variable.name = "clim_var",
         value.name = "value_std") %>% 
    left_join(sp.x$mean_sd, by = "clim_var") %>%
    mutate(value_nonstd = value_std * std_sd + std_mean) %>%
    select(plot_id, p_a, clim_var, value_nonstd) %>%
    # dcast(plot_id + p_a ~ clim_var, value.var = "value_nonstd") %>%
    # select(-plot_id) %>%
    group_by(p_a, clim_var) %>%
    summarize(v_med = median(value_nonstd),
              v_min = min(value_nonstd),
              v_max = max(value_nonstd),
              v_q05 = quantile(value_nonstd, 0.05) %>% as.numeric,
              v_q95 = quantile(value_nonstd, 0.95) %>% as.numeric)
})

# == parameter trace & densities ========================================= ####
# Note: these have 95% blue quantile intervals on them 
#       (atm, q's hard-coded in plotting fns...)
# - plots: trace & density ----------------------------------------------- ####
# trace plots
list(mod_spp, beta_crosswalk) %>% 
  pmap(function(mod.x, crosswalk.x) {
  mod.x %>% modify_at("MCMC_theta", function(theta_df) {
    # get bio* style names...
    update_index <- match(crosswalk.x$mod_name, names(theta_df))
    names(theta_df)[update_index] <- crosswalk.x$orig_names

    return(theta_df)
  })
}) %>% 
  list(., dat_ls) %>% 
  pmap(~fn.plot_trace_theta(.x, NA, plot_actual = F, m_name = .y$sp.x))  

# density plots
list(mod_spp, beta_crosswalk) %>% 
  pmap(function(mod.x, crosswalk.x) {
    mod.x %>% modify_at("MCMC_theta", function(theta_df) {
      # get bio* style names...
      update_index <- match(crosswalk.x$mod_name, names(theta_df))
      names(theta_df)[update_index] <- crosswalk.x$orig_names
      
      return(theta_df)
    })
  }) %>% 
  list(., dat_ls) %>% 
  pmap(~fn.plot_density_theta(.x, NA, plot_actual = F, m_name = .y$sp.x))
# - eq: posterior means for FE model form eq's --------------------------- ####
list(mod_spp %>% map(~.x$MCMC_theta), beta_crosswalk) %>% 
  pmap(function(theta.x, crosswalk.x) {
    # get bio* style names...
    update_index <- match(crosswalk.x$mod_name, names(theta.x))
    names(theta.x)[update_index] <- crosswalk.x$orig_names
    rm(update_index)
    
    # tabled summaries 
    theta.x %>%
      select(-n_sim) %>%
      melt(variable.name = "parameter") %>%
      group_by(parameter) %>%
      summarize(post_mean = round(mean(value), 2)) %>%
      ungroup %>%
      select(parameter, post_mean) %>%
      filter(!parameter %in% c("rho", "sigma_z", "sigma_v")) %>%
      mutate(value = case_when(
        parameter != "intercept" ~ paste0(post_mean, parameter),
        TRUE ~ paste0(post_mean))) %>%
      select(-post_mean) 
  })

# - table: summary of mean & 90% intervals ------------------------------- ####
table_theta <- list(mod_spp %>% map(~.x$MCMC_theta), beta_crosswalk) %>% 
  pmap(function(theta.x, crosswalk.x) {
    # get bio* style names...
    update_index <- match(crosswalk.x$mod_name, names(theta.x))
    names(theta.x)[update_index] <- crosswalk.x$orig_names
    rm(update_index)
    
    # tabled summaries 
    theta.x %>%
      select(-n_sim) %>%
      melt(variable.name = "parameter") %>%
      group_by(parameter) %>%
      summarize(expected = mean(value),
                q_025 = quantile(value, .05) %>% as.numeric,
                q_975 = quantile(value, .95) %>% as.numeric) %>%
      mutate_if(is.numeric, round, 2) %>%
      ungroup %>%
      mutate(exp_quant = paste0(expected, " (", q_025, ", ", q_975, ")")) %>%
      select(parameter, exp_quant)
  }) %>% 
  bind_rows(.id = "spp") %>%
  dcast(parameter ~ spp, value.var = "exp_quant", fill = "") %>%
  mutate(parameter = factor(parameter,
                            levels = c("intercept", 
                                       "gdd5",
                                       grep("bio", levels(parameter), value = T),
                                       "rho",
                                       "sigma_z",
                                       "sigma_v"))) %>%
  arrange(parameter)

# == Model checks - part 1 of 2 for Table 3 ============================== ####
# Table 3: ESS and AUC in this script, for the Moran's I model-check see
# script 12-...R

# - effective sample size ------------------------------------------------ ####
# Checking that the effective sample size is large enough...
# (univariate ESS)
ESS_uni <- mod_spp %>% map(function(mod.x) {
  # What is the univariate ESS for each parameter & species?
  # - note: set imse & verbose `F` b/c `T` gives warnings on lags
  fn.ESS_wrapper <- function(parm.x) {
    batchmeans::ess(parm.x, 
                    # (T) removes corr beyond a certain lag to reduce noise
                    imse = FALSE, 
                    # (T & imse == T) prints lag and warns if too small
                    verbose = FALSE)
  }
  
  univarESS <- list(
    MCMC_theta = apply(mod.x$MCMC_theta %>% select(-n_sim), 2, fn.ESS_wrapper),
    MCMC_z = apply(mod.x$MCMC_z, 2, fn.ESS_wrapper),
    MCMC_v = apply(mod.x$MCMC_v, 2, fn.ESS_wrapper))
  
}) 

ESS_uni %>%  # table summary of min (or med) ESS
  map(function(x) {
    theta.x <- x$MCMC_theta 
    beta_id <- grepl("beta", names(theta.x))
    z.x <- x$MCMC_z
    v.x <- x$MCMC_v
    
    c(sum(beta_id), min(theta.x[beta_id]), min(theta.x[beta_id][-1]), 
      median(theta.x[beta_id]), median(theta.x[beta_id][-1]),
      theta.x[!beta_id],
      length(z.x), min(z.x), median(z.x),
      length(v.x), min(v.x), median(v.x)) %>% 
      setNames(c("n_beta", "beta_min", "beta_min.wo_int", 
                 "beta_median", "beta_median.wo_int", 
                 names(theta.x[!beta_id]),
                 "n_z", "z_min", "z_median", 
                 "n_v", "v_min", "v_median")) %>% 
      bind_rows
  }) %>% 
  bind_rows(.id = "spp")

# - AUC (pROC::roc) ------------------------------------------------------ ####
# need to load preds_PRISM
preds_PRISM <- readRDS(paste0(path_output, "ch2-7-preds-PRISM.rds"))$preds_PRISM

# get AUC for each MCMC sample...
auc_ls <- preds_PRISM %>% map(function(preds_ls) {
  obsv_FIA_plots <- preds_ls$newdat.x %>%
    filter(y_index == 1) %>% # only observed FIA plots
    pull(plot_id)
  
  # vector of observed p_a, sorted by plot_id...
  y_obsv <- preds_ls$newdat.x %>% 
    filter(y_index == 1) %>% # only observed FIA plots
    arrange(plot_id) %>%
    pull(p_a)
  
  # get AUC for each MCMC map
  preds_ls$maps_1k %>%
    map(function(preds.x) {
      # get mu_pred for only observed FIA plots, sorted by plot_id...
      mu_pred <- preds.x %>%
        filter(plot_id %in% obsv_FIA_plots) %>%
        arrange(plot_id) %>%
        pull(mu_hat)
      
      # get AUC for MCMC sample
      (pROC::roc(response = as.numeric(y_obsv), 
                 predictor = mu_pred, 
                 levels = c(0, 1), 
                 direction = "<",
                 quite = TRUE)
      )$auc[1] %>%
        as.numeric
    }) %>%
    unlist
}) 

# summarize these AUC values...
# - min/max/med
auc_ls %>%
  map(~c(min(.x), median(.x), max(.x)) %>% 
        setNames(c("min", "med", "max")) %>% 
        round(2)) %>% 
  bind_rows(.id = "spp")

# clean-up space 
rm(preds_PRISM)

# == z-maps ============================================================== ####
# maps of the spatial random effects... below I'll show the z-maps for:
# - the mean of the MCMC samples, which is used in calculating mean-preds
# - the z's from the last MCMC sample drawn

# taking mean of z's
z_mean_map <- list(mod_spp, hex_polygons_public.sf, names(mod_spp)) %>% 
  pmap(~fn.map_z_sample_sf(..1, ..2, ..3, 
                           study_area.sf,
                           MCMC_iter = 1000, # 
                           expect_T = TRUE,
                           color_mu_T = TRUE))

# sample z's
z_samp_map <- list(mod_spp, hex_polygons_public.sf, names(mod_spp)) %>% 
  pmap(~fn.map_z_sample_sf(..1, ..2, ..3, 
                           study_area.sf,
                           MCMC_iter = 1000,
                           expect_T = FALSE,
                           color_mu_T = TRUE))

# arrange all in one figure
plot_grid(plot_grid(plotlist = z_mean_map, nrow = 1),
          plot_grid(plotlist = z_samp_map, nrow = 1),
          nrow = 2)

# == ES info (Appendix S1, Section S1 & other ref) ======================= ####
# - ES membership table (Appendix S1, Section S1) ------------------------ ####
mod_spp %>% 
  map(~data.frame("ES" = .x$MCMC_v %>% names %>% gsub("ES", "", .))) %>% 
  bind_rows(.id = "spp") %>% 
  dcast(ES ~ spp, fill = "") %>% 
  mutate(ES = factor(ES, levels = c(
    "242A","242B","M242A","M242B", "M242C","M242D","261A","261B",
    "M261A", "M261B","M261C","M261D","M261E","M261F","M261G","262A",
    "M262A","M262B","263A","M322G", "M332G","M333A","341D","341F",
    "342B","342H","342I","322A", "322C", "331A"))) %>% 
  arrange(ES) %>%
  mutate_at(vars(names(mod_spp)), ~ifelse(.x == "", .x, "X"))

# - map of ES by spp for quick reference... ------------------------------ ####
list(plot_polygons_public.sf, dat_ls %>% map(~.x$moddat)) %>%
  pmap(~.x %>% 
         inner_join(.y %>% select(plot_id, ES), by = "plot_id") %>%
         select(ES)) %>%
  imap(~.x %>%
         ggplot(aes(geometry = geometry)) +
         geom_sf(data = study_area.sf, fill = "white") +
         geom_sf(aes(fill = ES), color = NA) +
         theme_bw() +
         ggtitle(.y)) %>%
  pluck("qudo")

# == Conceptual diagrams (Figure 1 & 2) ================================== ####
# - Figure 1: rel. of plots to hex --------------------------------------- ####
# find three neighboring polygons...
hex_sf <- hex_polygons_public.sf$quga4[c(4, 2, 11),] # three adj polygons

# randomly insert 'plots' -- prettiest representation for set.seed
set.seed(9)
points_sf <- hex_sf %>%
  split(.$fiahex_id) %>%
  list(., c(1, 5, 1), list(c("observed"), 
                           c("observed", "observed", "unobserved", "observed", "unobserved"), 
                           c("unobserved"))) %>%
  pmap(function(hex.x, n.x, y.x) {
    points.x <- hex.x %>% 
      st_sample(., size = n.x, type = "random", exact = TRUE) %>%
      st_sf(.)
    points.x$response <- y.x
    
    points.x
  }) %>%
  setNames(hex_sf$fiahex_id) %>%
  bind_rows(.id = "fiahex_id")

# add new label names to EX hexagons and plots...
hex_sf$z <- paste0("z[", c(1,2,3), "]")
points_sf$y <- paste0("y[", c(4,2,1,6,3,5,7), "]")

# Plot example relationship between sample & spatial units
# - sample units (FIA plot) can be observed or unobserved, and multiple
#   sample units may occur within a single spatial unit (hexagon)
p_body <- ggplot() +
  geom_sf(data = hex_sf, fill = "lightgrey", size = 1) +
  geom_sf(data = points_sf, aes(color = response), size = 2) +
  geom_sf_text(data = points_sf, aes(label = y), parse = TRUE,
               nudge_x = c(450, -450, -450, -450, 450, -400, -400),
               nudge_y = c(200, 100, 100, 100, -100, -250, 300),
               # re-scales so multiply by 5/14 to get pt font we desire...
               size = 7 * (5/14)) +
  geom_sf_text(data = hex_sf, aes(label = z), parse = TRUE,
               nudge_x = 1500, nudge_y = 1500,
               size = 9 * (5/14)) +
  scale_color_manual(values = c("#39568cff", "#EB5760FF")) +
  theme_void() + 
  theme(legend.title = element_text(size = EJ_axis_title_size),
        legend.text  = element_text(size = EJ_strip_text_size))

ggsave(filename = paste0(path_images, "Fig1.tiff"),
       plot = p_body,
       height = 2.47*(EJ_dim_1col_width_max/3.7), 
       width = EJ_dim_1col_width_max,
       units = EJ_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")


# - Figure 2: rel. of persistence/expansion/contraction ------------------ ####
p_body <- data.frame(x = c(10, 12),
                     y = c(10, 17),
                     r = c(10, 9),
                     pred = c("current", "future")) %>%
  mutate(geometry = map2(x, y, ~st_point(c(.x, .y)))) %>%
  st_as_sf() %>%
  st_buffer(dist = .$r, nQuadSegs = c(5,2)) %>%
  st_intersection() %>%
  mutate(area = st_area(.) %>% round(-1) %>% as.numeric) %>%
  mutate(area = paste0(area, " km^2")) %>%
  mutate(peca_n = row_number()) %>%
  mutate(peca = case_when(peca_n == 1 ~ "Contract",
                          peca_n == 2 ~ "Persist",
                          peca_n == 3 ~ "Expand")) %>%
  mutate(peca = factor(peca, levels = c("Expand", "Persist", "Contract"))) %>%
  ggplot() +
  geom_sf(aes(fill = peca)) +
  geom_sf_label(aes(label = area), 
                nudge_x = c(1, 0, -2),
                nudge_y = c(-2, 0, 2),
                size = 7 * (5/14)) +
  scale_fill_brewer(palette = "YlGnBu") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text  = element_text(size = EJ_strip_text_size))

ggsave(filename = paste0(path_images, "Fig2.tiff"),
       plot = p_body,
       height = 2.63*(EJ_dim_1col_width_max/3), 
       width = EJ_dim_1col_width_max,
       units = EJ_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")
  
# == How many plot/hex dropped b/c of isolated hex? ====================== ####
list(data_spp.std %>% map(~.x$dat %>% select(plot_id, fiahex_id)),
     dat_ls %>% map(~.x$moddat %>% select(plot_id, fiahex_id))) %>%
  pmap(~.x %>%
         filter(!plot_id %in% .y$plot_id) %>%
         summarize(n_hex = n_distinct(fiahex_id),
                   n_plots = n(),
                   PoT_hex = n_distinct(fiahex_id) / n_distinct(.x$fiahex_id),
                   PoT_plots = n() / nrow(.x))) %>%
  bind_rows(.id = "spp")




