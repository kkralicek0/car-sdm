###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 3-explore-data.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@udsa.gov)
##
## This R script is in support of 'Climate change induced shifts in suitable
## habitat projected for PNW tree species with spatial-Bayesian models' by 
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - explore biologically meaningful climate variables
## - identify important vars to spp
## - create table summaries and plots/maps
##
## About - output:
## - saved figures from this script in path_images for easy inspection outside
##   of R... done via ggsave
###############################################################################
library(plyr)         # for rbind.fill(), load before dplyr to avoid masking
library(dplyr)        # for ...E V E R Y T H I N G...
library(magrittr)     # all the pipes (e.g. %$%)
library(stringr)      # for string help (e.g. str_detect in functions)
library(reshape2)     # melt() and dcast()
library(purrr)        # working with lists
library(ggplot2)      # for plotting
library(gridExtra)    # for grids of plots
library(GGally)       # for ggpairs()
library(viridis)      

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_images <- paste0(path_top, "images/2-explore-data/")

# == Load data =========================================================== ####
load(paste0(path_data, "ch2-2-data.Rdata"))
spp <- names(FIA_spp)

# == FIA_spp: overall summaries (n_plot, prev, hex, obsv/unobsv) ========= ####
# - Overall summaries by spp {n-plots, prevalence(%)} ----
# (all data)
FIA_spp %>% map(function(x) {
  x %>%
    summarize(n_plots = n(),
              n_ES = n_distinct(ES),
              perc_prev = sum(p_a)/n() * 100 %>% round(2))
}) %>% bind_rows(.id = "spp")

# (observed v. unobserved response - this will only change for qudo & quga4)
# - plot-level
FIA_spp %>% map(function(x) {
  x %>% 
    group_by(plot_obsv) %>%
    summarize(n_plots = n())
}) %>% bind_rows(.id = "spp") %>%
  dcast(spp ~ plot_obsv)

# - hex-level
FIA_spp %>% map(function(x) {
  x %>% 
    group_by(fiahex_obsv) %>%
    summarize(n_plots = n())
}) %>% bind_rows(.id = "spp") %>%
  dcast(spp ~ fiahex_obsv)

# - Nonforest plots/hex breakdown (qudo & quga4 only) ----
FIA_spp[c("qudo", "quga4")] %>% map(function(x) {
  x %>% 
    group_by(fiahex_id) %>%
    mutate(n_plots_hex = n()) %>%
    ungroup() %>%
    filter(n_plots_hex > 1) %>%
    group_by(fiahex_id) %>%
    summarise(n_plots_obsv = sum(plot_obsv),
              n_plots_unobsv = sum(abs(plot_obsv - 1))) %>%
    ungroup() %>%
    group_by(n_plots_obsv, n_plots_unobsv) %>%
    summarize(n_hex = n())
}) %>% bind_rows(.id = "spp")

# == Explore: box-density|ggpairs ======================================== ####
# - format data ---------------------------------------------------------- ####
# group variables for easy viewing
vars_clim <- c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio13", "bio16", "bio18")

# create more managable dfs for plotting
vars_group <- c("plot_id", "LAT_ACTUAL", "LON_ACTUAL", 
                "ES", "p_a")
plot_ls <- FIA_spp %>% map(~.x %>% select(all_of(c(vars_group, vars_clim))))
plot_df <- plot_ls %>% bind_rows(.id = "spp")
plot_df_long <- plot_df %>% melt(id.vars = c(vars_group, "spp"))

# - box/density plot ----------------------------------------------------- ####
(plot_df_long %>% 
  ggplot() + 
  geom_boxplot(aes(y = value, x = spp), alpha = .05, width = 0.2) +
  geom_violin(aes(y = value, x = spp, fill = spp), alpha = .5) +
  facet_wrap(~variable, scales = "free") + 
  scale_fill_viridis_d()) %>% 
  ggsave(filename = paste0(path_images, "boxdensity.jpg"),
         plot = .,
         device = "jpeg",
         width = 20,
         height = 10)

# - ggpairs plots -------------------------------------------------------- ####
# Need to do in stages b/c too many variables

# define function: assign dif colors for P & A 
# (this is complex so fn to clean it up)
fn.ggpairs_new_colors <- function(data, group_by_var, 
                                  colors_for_data, alpha_val, ...) {
  ggpairs(data = data,
          mapping =  aes_string(colour = group_by_var),
          lower   = list(
            continuous = function(data, mapping, ...){
              ggally_points(data = data, mapping = mapping, alpha = alpha_val) + 
                scale_color_manual(values = colors_for_data)},
            combo = function(data, mapping, ...){
              ggally_facethist(data = data, mapping = mapping) + 
                scale_fill_manual(values = colors_for_data)}),
          diag = list(
            continuous = function(data, mapping, ...){
              ggally_densityDiag(data = data, mapping = mapping, alpha = alpha_val*10) +
                scale_fill_manual(values = colors_for_data)},
            discrete = function(data, mapping, ...){
              ggally_barDiag(data = data, mapping = mapping) +
                scale_fill_manual(values = colors_for_data)}),
          upper = list(
            continuous = function(data, mapping, ...){
              ggally_cor(data = data, mapping = mapping) +
                scale_color_manual(values = colors_for_data)},
            combo = function(data, mapping, ...){
              ggally_box_no_facet(data = data, mapping = mapping) +
                scale_fill_manual(values = colors_for_data)})) +
    theme(panel.grid.major = element_blank(),
          panel.border = element_rect(linetype = "solid", colour = "black", fill = NA))
}

# specify vars of interest for each species based on physiology:
# - from silvics of north america & Devine et al. (2012) tree profiles,
#   both citations in the manuscript's lit cited
pairs_vars <- list(
  "abpr"   = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio16", "bio18"),
  "psmem"  = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio16", "bio18"),
  "qudo"   = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio13", "bio16", "bio18"),
  "quga4"  = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio16", "bio18"),
  "quke"   = c("hmi", "gdd5", "bio5", "bio6", "bio8", "bio10", "bio12", "bio13", "bio16", "bio18"))

# plot and save with pmap...
list(FIA_spp, pairs_vars, spp) %>% 
  pmap(function(df.x, vars.x, sp.x) {
    # subset to variables of interest for this spp
    df.x <- df.x %>% select(all_of(c(vars.x, "p_a")))
    
    # Change P-A to non-numeric values (to appease new_colors fn)
    df.x$p_a <- as.character(df.x$p_a)
    df.x[df.x$p_a == 0, ]$p_a <- "absent"
    df.x[df.x$p_a == 1, ]$p_a <- "present"
    df.x <- df.x %>% arrange(p_a)
    
    # plot and save
    (fn.ggpairs_new_colors(data = df.x,
                           group_by_var = "p_a",
                           colors_for_data = c("#252525", "#006d2c"),
                           alpha_val = 0.05) +
        ggtitle(sp.x)) %>%
      ggsave(filename = paste0(path_images,
                               paste0("ggpairs_", sp.x, ".jpg")),
             plot = .,
             device = "jpeg",
             width = 20,
             height = 10)
  })

# == Explore: maps of study area ========================================= ####
# a rough take at looking at PRISM data, just to get a feel for things
# - create base_map of OR/WA/CA ----
map_pw <- map_data("state") %>%
  filter(region %in% c("oregon", "washington", "california")) %>%
  fortify

base_map <- ggplot() +
  geom_polygon(data = map_pw, aes(x = long, y = lat, group = group),
               fill = "white",
               color = "black") +
  coord_map() # Mercator projection by default, you can change this though

# - define plot/save fns ----
fn.p_map <- function(df, clim_var, title, savetitle, lnish_TF = FALSE) {
  (if (lnish_TF) {
    base_map  +
      geom_point(data = df %>% mutate(ln25 = log(!!as.name(clim_var) + 25)), 
                 aes_string(x = "LON_ACTUAL", y = "LAT_ACTUAL", color = "ln25"), 
                 alpha = .25) +
      scale_colour_viridis_c() +
      ggtitle(title)
  } else {
    base_map  +
      geom_point(data = df, 
                 aes_string(x = "LON_ACTUAL", y = "LAT_ACTUAL", color = (clim_var)), 
                 alpha = .25) +
      scale_colour_viridis_c() +
      ggtitle(title)
  }) %>%
    ggsave(filename = paste0(path_images, paste0("map_", savetitle, ".jpg")),
           plot = .,
           device = "jpeg",
           height = 10)
  
}
# - map: hmi(s) for all data (FIA) ----
fn.p_map(FIA, "hmi", "annual heat-moisture index", "hmi")
fn.p_map(FIA, "gdd5", "growing degree days > 5C", "gdd5")
fn.p_map(FIA, "bio5", "max temp of warmest month", "bio5")
fn.p_map(FIA, "bio6", "min temp of coldest month", "bio6")
fn.p_map(FIA, "bio8", "mean temp of wettest quarter", "bio8")
fn.p_map(FIA, "bio10", "mean temp of warmest quarter", "bio10")
fn.p_map(FIA, "bio12", "annual precipitation", "bio12")
fn.p_map(FIA, "bio13", "ppt of the wettest month", "bio13")
fn.p_map(FIA, "bio16", "ppt of wettest quarter", "bio16")
fn.p_map(FIA, "bio18", "ppt of warmest quarter", "bio18")

