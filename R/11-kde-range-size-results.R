###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 11-kde-range-size-results.R
## Updated - 01-28-2023
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Spatial-Bayesian models project shifts in 
## suitable habitat for Pacific Northwest tree species under climate change' by
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - get range-size estimates based on BKDE range-estimation method
## - calculate other statistics and compare range estimates between the BKDE 
##   and VP methods
## - make figures & table summaries for paper & appendices:
##   - Manuscript: Figures 6 & 7
##   - Appendix S1: Figure S1 and Tables S3-S5 
##
## About - output:
## - ch2-11-subset-pps-PRISM.rds & ch2-11-subset-pps-CNA.rds
##  - location: path_output
##  - b/c original prediction RDS's were so large, needed to work with subset
##    for BKDE calculations... these objects are just the MCMC samples' dfs
##    subset-ed to just the plot_id & pps_pa info
## - ch2-11-kde-current.rds
##  - location: path_output
##  - BKDE range-size estimates based on 95% kde contour polygon
## - ch2-11-kde-future-diff.rds
##  - location: path_output
##  - BKDE range-size estimates & persist/expand/contract based on current
##    estimate for the same MCMC sample
###############################################################################
# data manipulation
library(magrittr)    # for set_rownames()
library(dplyr)       # for ... E V E R Y T H I N G ...
library(purrr)       # working with lists
library(reshape2)    # melt() and dcast()
library(stringr)     # for string help (e.g. str_detect in functions)

# plotting
library(ggplot2)
library(patchwork)   # overrides align_plots from cowplot pkg...

# voronoi...
library(sf)

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

# == Load data =========================================================== ####
# - (create|commented-out) pps-only RDS ---------------------------------- ####
# # About - 
# # - prediction RDS are LARGE... to reduce the memory taxation, creating a 
# #   subsetted version of this data to work with...
# # - if exists statement to load if already made, if not, make, close, and
# #   then open/load/continue...
# 
# # Need these to create kde rds... if those already made, don't load b/c
# # these still take up quite a lot of space
# 
# if (file.exists(paste0(path_output, "ch2-11-subset-pps-PRISM.rds"))) {
#   pps_PRISM <- readRDS(paste0(path_output, "ch2-11-subset-pps-PRISM.rds"))
# } else {
#   # current climate predictions
#   pps_PRISM <- readRDS(paste0(path_output, "ch2-7-preds-PRISM.rds"))$preds_PRISM
#   
#   # drop plot ids
#   drop_plot_ids <- readRDS(paste0(path_output, "ch2-8-drop_plot_ids.rds"))
#   
#   # remove drop-plot_id's & subset to only the data we need
#   pps_PRISM <- list(pps_PRISM, drop_plot_ids) %>%
#     pmap(function(preds.x, drop_ids.x) {
#       preds.x$maps_1k %>%
#         map(~.x %>% 
#               filter(!plot_id %in% drop_ids.x) %>% 
#               select(plot_id, pps_pa))
#     })
#   
#   # save (as it didn't already exist)
#   saveRDS(pps_PRISM, paste0(path_output, "ch2-11-subset-pps-PRISM.rds"))
#   
#   # clean-up space
#   rm(drop_plot_ids)
# }
# 
# if (file.exists(paste0(path_output, "ch2-11-subset-pps-CNA.rds"))) {
#   pps_CNA <- readRDS(paste0(path_output, "ch2-11-subset-pps-CNA.rds"))
# } else {
#   # future climate predictions
#   # - names of future scenarios
#   future_scenarios <- readRDS(paste0(path_output, "ch2-7-future_scenarios_chr.rds"))
# 
#   # - read-in predictions!
#   pps_CNA <- lapply(future_scenarios, function(x) {
#     readRDS(paste0(path_output, "ch2-7-preds-CNA_", x, ".rds"))$preds_CNA
#   }) %>%
#     setNames(future_scenarios)
#   
#   # drop plot ids
#   drop_plot_ids <- readRDS(paste0(path_output, "ch2-8-drop_plot_ids.rds"))
#   
#   # remove drop-plot_id's & subset to only the data we need
#   pps_CNA <- pps_CNA %>%
#     map(~list(.x, drop_plot_ids) %>%
#           pmap(function(preds.x, drop_ids.x) {
#             preds.x$maps_1k %>%
#               map(~.x %>% 
#                     filter(!plot_id %in% drop_ids.x) %>% 
#                     select(plot_id, pps_pa))
#           }))
#   
#   # save (as it didn't already exist)
#   saveRDS(pps_CNA, paste0(path_output, "ch2-11-subset-pps-CNA.rds"))
#   
#   # clean-up space
#   rm(drop_plot_ids)
#   
# }
# 
# - fns + spatial data --------------------------------------------------- ####
# load fns
source(paste0(path_fn, "0-functions.R"))

# spatial data 
plot_points_true.sf  <- readRDS(paste0(path_output, "ch2-8-plot_points_true_sf.rds"))
hex_polygons_true.sf <- readRDS(paste0(path_output, "ch2-8-hex_polygons_true_sf.rds"))
study_area.sf <- readRDS(paste0(path_output, "ch2-8-study_area_sf.rds"))

# == Define projections (crs) ============================================ ####
# Primary planar projection:
# - [distance preserving] UTM Zone 11N (EPSG: 32611)
crs_UTM11N <- 32611

# Area planar projection:
# - [area preserving] Albers Equal Area Conic (EPSG: 5070)
# - only used temporarily for calculations, data projected back to primary
#   projection post.
crs_Albers <- 5070

# == Arguments for figures =============================================== ####
scen_colPal <- data.frame(
  scenario = c("current",
               "CCSM4_rcp45_m2085",
               "HadGEM2_ES_rcp45_m2085",
               "CCSM4_rcp85_m2085",
               "HadGEM2_ES_rcp85_m2085"),
  # colorblind safe (cb) diverging palettes 
  # - reds (H) and blues (C)
  hex_cb_br = c("#999999", 
                "#56B4E9",  "#E69F00",  "#0072B2",  "#D55E00"),
  # - purple (H) and greens (C)
  # This is the one we are really going with...
  hex_cb_gp = c("#999999", 
                "#af8dc3",  "#7fbf7b", "#762a83",  "#1b7837"),
  hex_cb_mix = c("#999999",
                 "#56B4E9", "#E69F00", "#009E73", "#F0E442"),
  # - blue (C) and green (H)
  hex_cb_gb = c("#999999",
                "#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")
)

# common names for species
spp_common <- c("noble fir",
                "coastal Douglas-fir",
                "blue oak",
                "white oak",
                "black oak") %>%
  setNames(names(plot_points_true.sf))

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

# == BKDE: range size polygons & estimates =============================== ####
# - (create|load) kde range size estimates ------------------------------- ####
if (file.exists(paste0(path_output, "ch2-11-kde-current.rds")) &
    file.exists(paste0(path_output, "ch2-11-kde-future-diff.rds"))) {
  # kde_current <- readRDS(paste0(path_output, "ch2-11-kde-current.rds"))
  kde_future_diff <- readRDS(paste0(path_output, "ch2-11-kde-future-diff.rds"))
} else {
  # - define functions --------------------------------------------------- ####
  # fn.kde_args: gets sp-specific global arguments for the kde fn
  fn.kde_args <- function(ppoints_sf.x, hex_sf.x,
                          study_area.sf = study_area.sf,
                          approx_cell_area_ha = 2500/5,
                          mbuffer = 40) {
    # This fn returns: ----
    # - H: bandwidth matrix for a species. Choice of H is determined based on
    #   observed presences (original data set), aggregated to a hexagon-level to
    #   avoid over-weighting areas where FIA plots were spatially intensified.
    #   The same H is used with each of the 2k MCMC maps for each scenario.
    # - boundary.x: The species-specific study area boundary. Intersection of 
    #   the states boundary with the union of in-sample hexagons for a given spp
    # - bbox_buff: bounding box of boundary.x plus a buffer. A four element 
    #   vector describing the area over which the kde grid should be evaluated
    # - n_xy: A 2 element vector with the dimensions of the kde grid-size
    
    # == Define sp-specific study area & grid-size for kde ================= ####
    # boundary for this species
    # - combination of in-sample hexagons and state's boundary
    boundary.x <- hex_sf.x %>% st_union %>% st_intersection(., study_area.sf)
    
    # grid-size & min/max values for the grid with which to est the density
    # (default is a 151 x 151 grid... since our bounding boxes are not square
    #  we probably want points more evenly spaced N-S and E-W)
    # - want to add a buffer to our bounding box so that full 95% density polygon
    approx_cell_len_m <- approx_cell_area_ha^.5 * 100
    bbox_buff <- st_bbox(boundary.x) %>% as.numeric
    bbox_buff <- c(bbox_buff[1:2] - approx_cell_len_m * mbuffer,
                   bbox_buff[3:4] + approx_cell_len_m * mbuffer)
    len_xy <- c(bbox_buff[3:4] - bbox_buff[1:2])
    n_xy   <- ceiling(len_xy / approx_cell_len_m)
    
    # clean-up space
    rm(approx_cell_len_m, len_xy)
    
    # == Choose bandwidth matrix (H) for all scenarios/samples ============= ####
    # Filter to presence data & aggregate to hexagon level 
    # - Working with observed presences here so that the H matrix is consistent
    #   across scenarios & MCMC samples (this way, the only uncertainty between
    #   range size polygons will be due to variability between MCMC maps & 
    #   scenarios)
    # - Because FIA plots are not spatially intensified consistently across all 
    #   states and ownership, aggregating any p-observations to the hexagon level  
    #   and associating them with the hex centroid.
    ppsp_coords <- hex_sf.x %>%
      st_set_agr(., "constant") %>% # specify attributes as spatially constant
      st_centroid() %>%
      mutate(x = unlist(map(.$geometry,1)),
             y = unlist(map(.$geometry,2)),
             fiahex_id = fiahex_id %>% as.character) %>%
      st_drop_geometry() %>%
      inner_join(ppoints_sf.x %>% st_drop_geometry,
                 by = "fiahex_id") %>%
      group_by(fiahex_id) %>%
      mutate(hex_p = as.numeric(sum(p_a) > 0)) %>%
      filter(hex_p == 1) %>%
      select(-p_a, -plot_id) %>%
      ungroup
    
    # Choose bandwidth matrix
    # - Using a 2-stage plug-in bandwidth selector with unconstrained pilot 
    #   bandwidth based on methods in Chacon & Duong (2010). From what I've
    # - From what I've read (Gitzen et al., 2006, C&D 2010, and others), a
    #   two-stage approach is preferable at reducing bias, but does tend to
    #   result in larger bandwidth estimates and more smoothing, but is still
    #   preferable to a 0-stage plug-in selector (like the ks::Hns, normal 
    #   scale selector). Also, the unconstrained pilot bandwidth method has
    #   shown in C&D2010's simulation study to perform better with densities
    #   similar to my species' distributions...
    
    # (originally I was specifying deriv.order = 2... but I think what I really
    #  wanted to specify was nstage, as I have done below...)
    H <- ks::Hpi(ppsp_coords %>% select(x, y), 
                 nstage = 2,
                 pilot = "unconstr")
    
    # == Return items ====================================================== ####
    list(H = H,
         boundary.x = boundary.x,
         bbox_buff = bbox_buff,
         n_xy = n_xy)
  }
  
  fn.kde_current <- function(preds.x, kde_args.x,
                             hex_plot_walk = hex_plot_walk) {
    # This fn is intended to work with the current scenario data and returns:
    # - kde95_sf: 95% kde contour polygon as a sf object for each MCMC sample
    # - kde95_area: estimate of the range size (ha) for each MCMC sample
    
    # map by MCMC sample to get kde, kde_95_poly, and then kde_95_area...
    preds.x %>%
      map(function(pred.i) {
        # == get baseline data for kde ===================================== ####
        # (atm dropping scenario & MCMC_sample_n)
        ppsp_coords <- pred.i %>% 
          left_join(hex_plot_walk, by = "plot_id") %>%
          group_by(fiahex_id) %>%
          mutate(hex_p = as.numeric(sum(pps_pa) > 0)) %>%
          filter(hex_p == 1) %>%
          select(-pps_pa, -plot_id) %>%
          ungroup
        
        # == Compute kde =================================================== ####
        # compute kernel density estimate based on pps presence and fixed H
        kde <- ks::kde(x = ppsp_coords %>% select(x, y), 
                       H = kde_args.x$H,
                       gridsize = kde_args.x$n_xy,
                       xmin = kde_args.x$bbox_buff[1:2],
                       xmax = kde_args.x$bbox_buff[3:4])
        
        # == Extract 95% kde polygon ======================================= ####
        # extract a contour containing 95% of the estimated density
        # - note: ks::kde() seems to label these backwards, so pulling '5%' below
        kde95_sf <- with(kde,
                         contourLines(x = eval.points[[1]],
                                      y = eval.points[[2]],
                                      z = estimate,
                                      levels = cont["5%"]))
        
        # remove potential holes from polygons & make sf object
        # - convert to sf object & intersect with sp-specific boundary area
        kde95_sf <- kde95_sf %>% 
          map(~with(.x, matrix(c(x, y), ncol = 2, byrow = FALSE))) %>%
          st_polygon %>% # this step to remove potential holes
          st_sfc %>%     # make into sf so we can specify correct crs
          st_set_crs(crs_UTM11N) %>% # ... this crs!
          st_make_valid %>% # (will freak out next step otherwise)
          st_intersection(., kde_args.x$boundary.x)
        
        # == Calculate area of range size estimate (ha) ==================== ####
        # here using Albers Equal Area Conic (planar to preserve area)
        kde95_area <- kde95_sf %>%
          st_transform(crs = crs_Albers) %>% 
          st_area(.) %>% 
          units::set_units(., ha) %>% 
          as.numeric 
        
        # == Return the items we need ====================================== ####
        list(kde95_sf = kde95_sf,
             kde95_area = kde95_area)
      })
  }
  
  # changed for the moment to return list...
  # - make into df post... there was an error suggesting an issue with one of 
  #   the arguments being length 0... checked kde_current and that is all good
  #   also checked some of the abpr future runs where the kde95_area is 0, and
  #   that code worked... so I think it is kde95_area_diff when there might
  #   not be an overlap between the current and future distributions.
  fn.kde_future_diff <- function(preds.x, kde_current.x, kde_args.x,
                                 hex_plot_walk = hex_plot_walk) {
    # This fn is intended to work with the future scenario data and current
    # scenario kde results to return a df where rows are MCMC samples and 5
    # columns are area (ha) summaries for:
    #   (1) current kde95_area
    #   (2) future kde95_area
    #   (3) persistence: suitable currently and suitable in the future
    #       {overlap between kde95_sf of current & future}
    #   (4) expansion: unsuitable currently, but suitable in the future
    #       {(2) - (3)}
    #   (5) contraction: suitable currently, but unsuitable in the future
    #       {(1) - (3)}
    # Special consideration:
    # - if all pps_pa == 0 for a scenario + MCMC sample combination, as has
    #   happened for some of the abpr samples before the truly final preds,
    #   then kde will throw error as ppsp_coords will be empty. So when this
    #   is the case, the range size estimate is zero and we can skip the sf
    #   steps
    
    # map by MCMC sample to get kde, kde_95_poly, and then kde_95_area...
    list(preds.x, kde_current.x) %>%
      pmap(function(pred.i, kde_current.i) {
        if (any(pred.i$pps_pa > 0)) {
          # - get baseline data for kde ---------------------------------- ####
          # (atm dropping scenario & MCMC_sample_n)
          ppsp_coords <- pred.i %>% 
            left_join(hex_plot_walk, by = "plot_id") %>%
            group_by(fiahex_id) %>%
            mutate(hex_p = as.numeric(sum(pps_pa) > 0)) %>%
            filter(hex_p == 1) %>%
            select(-pps_pa, -plot_id) %>%
            ungroup
          
          # - compute kde ------------------------------------------------ ####
          # compute kernel density estimate based on pps presence and fixed H
          kde <- ks::kde(x = ppsp_coords %>% select(x, y), 
                         H = kde_args.x$H,
                         gridsize = kde_args.x$n_xy,
                         xmin = kde_args.x$bbox_buff[1:2],
                         xmax = kde_args.x$bbox_buff[3:4])
          
          # - extract 95% kde polygon ------------------------------------ ####
          # extract a contour containing 95% of the estimated density
          # - note: ks::kde() seems to label these backwards, so pulling '5%' below
          kde95_sf <- with(kde,
                           contourLines(x = eval.points[[1]],
                                        y = eval.points[[2]],
                                        z = estimate,
                                        levels = cont["5%"]))
          
          # remove potential holes from polygons & make sf object
          # - convert to sf object & intersect with sp-specific boundary area
          kde95_sf <- kde95_sf %>% 
            map(~with(.x, matrix(c(x, y), ncol = 2, byrow = FALSE))) %>%
            map(~rbind(.x, .x[1,])) %>% # close polygons (just in case)
            st_polygon %>% # this step to remove potential holes
            st_sfc %>%     # make into sf so we can specify correct crs
            st_set_crs(crs_UTM11N) %>% # ... this crs!
            st_make_valid %>% # (will freak out next step otherwise)
            st_intersection(., kde_args.x$boundary.x)
          
          # - calc area estimates (ha) ----------------------------------- ####
          # here using Albers Equal Area Conic (planar to preserve area)
          # - future range size estimate (ha)
          kde95_area <- kde95_sf %>%
            st_transform(crs = crs_Albers) %>% 
            st_area(.) %>% 
            units::set_units(., ha) %>% 
            as.numeric 
          
          # - area estimated to persist
          kde_95_intersection <- kde95_sf %>%
            st_intersection(., kde_current.i$kde95_sf) %>% 
            st_transform(crs = crs_Albers) %>% 
            st_area(.) %>% 
            units::set_units(., ha) %>% 
            as.numeric 
        } else {
          # the case when all of the pps were absences...
          kde95_area <- 0
          kde_95_intersection <- 0
        }
        
        # # Return area summary df
        # data.frame(current_est = kde_current.i$kde95_area,
        #            future_est  = kde95_area,
        #            persistence = kde_95_intersection) %>%
        #   mutate(expansion = future_est - persistence,
        #          contraction = current_est - persistence)
        
        # return a list instead... 
        list(current_est = kde_current.i$kde95_area,
             future_est  = kde95_area,
             persistence = kde_95_intersection,
             expansion = kde95_area - kde_95_intersection,
             contraction = kde_current.i$kde95_area - kde_95_intersection)
      }) 
    # %>% 
    #   bind_rows %>%
    #   mutate(MCMC_sample_n = 1:n()) # because these are in order
  }
  
  # - make crosswalk between plots and their containing hexagons --------- ####
  # Crosswalk plot_id's to their containing hex id's & centroid coords
  hex_plot_walk <- list(plot_points_true.sf, hex_polygons_true.sf) %>%
    pmap(function(ppoints_sf.x, hex_sf.x) {
      hex_sf.x %>%
        st_set_agr(., "constant") %>% # specify attributes as spatially constant
        st_centroid() %>%
        mutate(x = unlist(map(.$geometry,1)),
               y = unlist(map(.$geometry,2)),
               fiahex_id = fiahex_id %>% as.character) %>%
        st_drop_geometry() %>%
        inner_join(ppoints_sf.x %>% st_drop_geometry,
                   by = "fiahex_id")
    }) %>%
    bind_rows() %>%
    distinct()
  
  # - choose H & get kde_95 for current ---------------------------------- ####
  kde_args <- list(plot_points_true.sf, hex_polygons_true.sf) %>%
    pmap(~fn.kde_args(..1, ..2, study_area.sf,
                      approx_cell_area_ha = 2500/5,
                      mbuffer = 40))
  
  # - will take ~2.2 hrs for all spp x 2k MCMC samples
  start_time <- Sys.time()
  kde_current <- list(pps_PRISM, kde_args) %>%
    pmap(~fn.kde_current(.x, .y, hex_plot_walk))
  run_time <- Sys.time() - start_time
  
  saveRDS(kde_current, 
          file = paste0(path_output, "ch2-11-kde-current.rds"))
  
  # - get kde_95 for future & difference from current -------------------- ####
  # Note: pps_CNA is organized by scenario, grab names with bind_rows
  # - will take ~20.3 hrs for all scenarios x spp x 2k MCMC samples
  start_time <- Sys.time()
  kde_future_diff <- pps_CNA %>%
    map(~list(.x, kde_current, kde_args) %>%
          pmap(~fn.kde_future_diff(..1, ..2, ..3, hex_plot_walk)))
  run_time2 <- Sys.time() - start_time
  
  saveRDS(kde_future_diff, 
          file = paste0(path_output, "ch2-11-kde-future-diff.rds"))
  
  
}

# clean up two case for future scenario predictions
# (1) some pps_pa were `1`, but still no area is allocated to a 95% polygon
#     under the future scenario. fix by setting future_est to 0 & calc the
#     change vars (only one case: MCMC_iter 1002 of abpr & rcp85_CCSM4)
# (2) fix where there is no predicted overlap between the current and future 
#     scenario's range estimates. This happened for abpr and prior to this fix 
#     appears as entries with estimates for the current and future 
#     distributions, but gives `numeric(0)` for persistence, expansion, 
#     and contraction --- which occurs if the intersection area is empty 
#     (e.g. no polygon to calculate the area from...)
kde_future_diff <- kde_future_diff %>%
  map_depth(3, function(x) {
    # here x is a list of results with elements:
    # {current_est, future_est, persistence, expansion, contraction}
    
    # case 1
    if (length(x$future_est) == 0) {
      x$future_est <- 0 # this is what was not filled
    }
    
    # case 2 & resolving rest of case 1
    # - note, if no area in future, the intersection will also be zero
    if (length(x$persistence) == 0) {
      kde_95_intersection <- 0 # this is what was not filled
      
      # fix p/e/c
      x$persistence <- kde_95_intersection
      x$expansion   <- x$future_est - kde_95_intersection
      x$contraction <- x$current_est - kde_95_intersection
    }
    
    # return updated/fixed list of results
    return(x)
  })

# into a list by spp of dfs that are easy to summarize/plot
kde_future_diff <- kde_future_diff %>% 
  transpose %>%
  # add MCMC iteration number to results & make df
  map_depth(2, function(scen.x) {
    list(scen.x, 1:2000) %>%
      pmap(function(results.x, MCMC_iter.x) {
        results.x$MCMC_iter <- MCMC_iter.x
        results.x %>% bind_rows
      }) %>% bind_rows 
  }) %>%
  # bind results for a single spp and add scenario variable
  map(~.x %>% bind_rows(.id = "scenario")) 

# - BKDE's PEC table (Appendix S1: Table S5) ----------------------------- ####
# means and 90% intervals...
# - (thousand km^2)
kde_future_diff %>%
  map(function(sp.x, round_digits = 3) {
    sp.x %>% 
      select(scenario, persistence, expansion, contraction) %>%
      melt(id.vars = c("scenario")) %>%
      mutate(value = (value/100)/1000,
             variable = factor(variable, 
                               levels = c("expansion", "persistence", "contraction"))) %>%
      group_by(scenario, variable) %>%
      summarize(mean = mean(value),
                q_05 = quantile(value, .05) %>% unname,
                q_95 = quantile(value, .95) %>% unname) %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      mutate(scenario = 
               case_when(scenario == "CCSM4_rcp45_m2085" ~ "RCP4.5 CCSM4",
                         scenario == "HadGEM2_ES_rcp45_m2085" ~ "RCP4.5 HadGEM2-ES",
                         scenario == "CCSM4_rcp85_m2085" ~ "RCP8.5 CCSM4",
                         scenario == "HadGEM2_ES_rcp85_m2085" ~ "RCP8.5 HadGEM2-ES")) %>%
      mutate(scenario =factor(scenario,
                              levels = c("RCP4.5 CCSM4",
                                         "RCP4.5 HadGEM2-ES",
                                         "RCP8.5 CCSM4",
                                         "RCP8.5 HadGEM2-ES"))) %>%
      mutate(value = paste0(mean, " (", q_05, ", ", q_95, ")")) %>%
      select(scenario, variable, value) %>%
      ungroup %>%
      dcast(variable ~ scenario, value.var = "value")
  }) %>%
  bind_rows(.id = "spp")

# - (as percent relative to the current range estimate)
kde_future_diff %>%
  map(function(sp.x, round_digits = 3) {
    sp.x %>% 
      mutate(persistence = persistence/current_est * 100,
             expansion = expansion/current_est * 100,
             contraction = contraction/current_est * 100) %>%
      select(scenario, persistence, expansion, contraction) %>%
      melt(id.vars = c("scenario")) %>%
      mutate(variable = factor(variable, 
                               levels = c("expansion", "persistence", "contraction"))) %>%
      group_by(scenario, variable) %>%
      summarize(mean = mean(value),
                q_05 = quantile(value, .05) %>% unname,
                q_95 = quantile(value, .95) %>% unname) %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      mutate(scenario = 
               case_when(scenario == "CCSM4_rcp45_m2085" ~ "RCP4.5 CCSM4",
                         scenario == "HadGEM2_ES_rcp45_m2085" ~ "RCP4.5 HadGEM2-ES",
                         scenario == "CCSM4_rcp85_m2085" ~ "RCP8.5 CCSM4",
                         scenario == "HadGEM2_ES_rcp85_m2085" ~ "RCP8.5 HadGEM2-ES")) %>%
      mutate(scenario =factor(scenario,
                              levels = c("RCP4.5 CCSM4",
                                         "RCP4.5 HadGEM2-ES",
                                         "RCP8.5 CCSM4",
                                         "RCP8.5 HadGEM2-ES"))) %>%
      mutate(value = paste0(mean, " (", q_05, ", ", q_95, ")")) %>%
      select(scenario, variable, value) %>%
      ungroup %>%
      dcast(variable ~ scenario, value.var = "value")
  }) %>%
  bind_rows(.id = "spp")

# - make summary figure (Figure 7) --------------------------------------- ####
# Figure that is 5 (spp) x 2, where 2 cols are:
# - posterior distributions of range-size estimates for the current and four
#   future climate scenarios (BKDE method)
# - scatter plot of area of contraction vs. expansion
p_body <- kde_future_diff %>% 
  imap(function(df.x, sp_name.x) {
    # - make new legend labels & create base data set ---------------------- ####
    # colors and old names...
    legend_labels <- with(scen_colPal, hex_cb_gp %>% setNames(scenario))
    future_lab_levels <- c("RCP 4.5 CCSM4",
                           "RCP 4.5 HadGEM2-ES",
                           "RCP 8.5 CCSM4",
                           "RCP 8.5 HadGEM2-ES")
    names(legend_labels) <- case_when(
      names(legend_labels) == "current" ~ "current",
      names(legend_labels) == "CCSM4_rcp45_m2085" ~ future_lab_levels[1],
      names(legend_labels) == "HadGEM2_ES_rcp45_m2085" ~ future_lab_levels[2],
      names(legend_labels) == "CCSM4_rcp85_m2085" ~ future_lab_levels[3],
      names(legend_labels) == "HadGEM2_ES_rcp85_m2085" ~ future_lab_levels[4])
    
    df.x <- df.x %>% 
      rowwise %>%
      mutate(scenario = names(legend_labels)[grep(scenario, scen_colPal$scenario)]) %>% 
      mutate(scenario = factor(scenario, levels = future_lab_levels)) %>%
      mutate_if(is.numeric, ~ (.x/100)/1000) %>% # convert from ha to 1k km^2
      arrange(MCMC_iter)
    
    # - make the base figures ---------------------------------------------- ####
    # get end points for 1-1 line segment
    line_1to1 <- min(max(df.x$contraction), max(df.x$expansion))
    std_margin <- margin(0, 0, 5.5, 0)
    abpr_margin <- margin(0, 0, 7.5, 0)
    
    # plot it!
    # - do a little differently for abpr (many near zero future ranges est)
    if (sp_name.x == "abpr") {
      p_range_est <- (df.x %>% 
                        select(current_est, MCMC_iter) %>%
                        mutate(scenario = "current") %>%
                        rename(range_est = current_est) %>%
                        select(scenario, range_est, MCMC_iter) %>%
                        rbind.data.frame(df.x %>% 
                                           rename(range_est = future_est) %>%
                                           select(scenario, range_est, MCMC_iter)) %>%
                        ungroup %>%
                        mutate(spp = recode(sp_name.x, !!!spp_common),
                               spp = factor(spp, levels = unname(spp_common)))) %>%
        ggplot() +
        geom_density(aes(x = log1p(range_est), 
                         y = log1p(..density..), 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ .) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_fill_manual(values = legend_labels,
                          labels = sub("^\\S+\\s+\\S+\\K\\s+", "\n", 
                                       names(legend_labels), perl=TRUE))  +
        xlab("log1p(range-size)") +
        ylab("log1p(density)") +
        theme_bw() +
        theme(plot.margin = abpr_margin,
              axis.title.x = element_text(size = EJ_axis_title_size,
                                          vjust = 3),
              strip.text = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size))
      
      p_exp_v_contract <- df.x %>%
        ggplot() +
        geom_segment(aes(xend = log1p(line_1to1), 
                         yend = log1p(line_1to1)), 
                     x = 0, y = 0) +
        geom_point(aes(x = log1p(contraction), 
                       y = log1p(expansion), 
                       color = scenario),
                   alpha = .25, shape = 20) + 
        expand_limits(x = 0, y = 0) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_y_continuous(position = "right") +
        ylab("log1p(expan.)") +
        theme_bw() +
        theme(plot.margin = abpr_margin,
              axis.title.x = element_text(size = EJ_axis_title_size,
                                          vjust = 3),
              strip.text = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size)) 
      
    } else {
      p_range_est <- (df.x %>% 
                        select(current_est, MCMC_iter) %>%
                        mutate(scenario = "current") %>%
                        rename(range_est = current_est) %>%
                        select(scenario, range_est, MCMC_iter) %>%
                        rbind.data.frame(df.x %>% 
                                           rename(range_est = future_est) %>%
                                           select(scenario, range_est, MCMC_iter)) %>%
                        ungroup %>%
                        mutate(spp = recode(sp_name.x, !!!spp_common),
                               spp = factor(spp, levels = unname(spp_common)))) %>%
        ggplot() +
        geom_density(aes(x = range_est, fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ .) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_fill_manual(values = legend_labels,
                          labels = sub("^\\S+\\s+\\S+\\K\\s+", "\n", 
                                       names(legend_labels), perl=TRUE))  +
        xlab("range-size") +
        theme_bw() +
        theme(plot.margin = std_margin,
              axis.title.x = element_text(size = EJ_axis_title_size),
              strip.text = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size))
      
      p_exp_v_contract <- df.x %>%
        ggplot() +
        geom_segment(aes(xend = line_1to1, yend = line_1to1), x = 0, y = 0) +
        geom_point(aes(x = contraction, y = expansion, color = scenario),
                   alpha = .25, shape = 20)  + 
        expand_limits(x = 0, y = 0) +
        scale_color_manual(values = legend_labels, guide = FALSE) +
        scale_y_continuous(position = "right") +
        theme_bw() +
        theme(plot.margin = std_margin,
              axis.title.x = element_text(size = EJ_axis_title_size),
              strip.text = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text  = element_text(size = EJ_axis_text_size,
                                        color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size))
    }
    
    # - return figs & base legend for the last ----------------------------- ####
    # do differently for quke
    # - quke will be the bottom figure, so the only one that needs x title 
    # - we'll also harvest the legend from this species figure
    if (sp_name.x == "quke") {
      # add some space to the left of the legend
      p_range_est <- p_range_est + 
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(legend.position = "bottom",
              legend.box.margin = margin(3, 0, 0, 0),
              legend.title = element_blank())
      
      list(p_range_est = p_range_est, 
           p_exp_v_contract = p_exp_v_contract + 
             guides(color = "none"))
    } else if (sp_name.x == "abpr") {
      list(p_range_est = p_range_est  + 
             guides(color = "none", fill = "none"), 
           p_exp_v_contract = p_exp_v_contract  + 
             guides(color = "none"))
    } else {
      list(p_range_est = p_range_est + 
             theme(axis.title.x = element_blank()) + 
             guides(color = "none", fill = "none"), 
           p_exp_v_contract = p_exp_v_contract + 
             theme(axis.title.x = element_blank()) + 
             guides(color = "none"))
    }
  }) %>%
  map(~ .x$p_range_est + .x$p_exp_v_contract) %>%
  reduce(., `/`)

# view it!
p_body + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Bivariate KDE range size estimates",
                  subtitle = "Based on 2k MCMC samples",
                  caption = "Note: abpr visualized on logp1 scale (ln(x + 1))") & 
  theme(legend.position = "bottom")

# save it
ggsave(filename = paste0(path_images, "Fig7.tiff"),
       plot = p_body + 
         plot_layout(guides = "collect") & 
         theme(legend.position = "bottom"),
       height = 5.7 * (EJ_dim_2col_width_max/6), 
       width = EJ_dim_2col_width_max,
       units = EJ_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")

# - Alternative figures for interpreting results ------------------------- ####
# persistence v contraction (with one-to-one line)
kde_future_diff %>% 
  imap(function(df.x, sp_name.x) {
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
    
    # - create base data set ----------------------------------------------- ####
    df.x <- df.x %>% 
      rowwise %>%
      mutate(scenario = legend_labels[grep(scenario, scen_colPal$scenario)]) %>% 
      mutate(scenario = factor(scenario, levels = legend_labels)) %>%
      mutate_if(is.numeric, ~ (.x/100)/1000) %>% # convert from ha to 1k km^2
      arrange(MCMC_iter)
    
    # - make the base figures ---------------------------------------------- ####
    # get end points for 1-1 line segment
    line_1to1 <- min(max(df.x$contraction), max(df.x$persistence))
    std_margin <- margin(0, 0, 5.5, 0)
    
    # plot it!
    # - do a little differently for abpr (many near zero future ranges est)
    if (sp_name.x == "abpr") {
      p <- df.x %>%
        ggplot() +
        geom_segment(aes(xend = log1p(line_1to1), yend = log1p(line_1to1)), x = 0, y = 0) +
        geom_point(aes(x = log1p(contraction), y = log1p(persistence), color = scenario),
                   alpha = .25, shape = 20) + 
        expand_limits(x = 0, y = 0) +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_y_continuous(position = "right") +
        ylab("log1p(pers.)") +
        theme_bw() +
        theme(plot.margin = std_margin)
      
    } else {
      p <- df.x %>%
        ggplot() +
        geom_segment(aes(xend = line_1to1, yend = line_1to1), x = 0, y = 0) +
        geom_point(aes(x = contraction, y = persistence, color = scenario),
                   alpha = .25, shape = 20)  + 
        expand_limits(x = 0, y = 0) +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_y_continuous(position = "right") +
        theme_bw() +
        theme(plot.margin = std_margin)
    }
    
    # - return figs & base legend for the last ----------------------------- ####
    # do differently for quke
    # - quke will be the bottom figure, so the only one that needs x title 
    # - we'll also harvest the legend from this species figure
    if (sp_name.x == "quke") {
      p + guides(color = "none")
    } else if (sp_name.x == "abpr") {
      p + guides(color = "none")
    } else {
      p + theme(axis.title.x = element_blank()) + 
        guides(color = "none")
    }
  }) %>%
  reduce(., `/`)

# stacked bar-chart with (exp, persist, contract & current est)
(kde_future_diff %>%
   map(function(sp.x) {
     sp.x %>% 
       select(scenario, persistence, expansion, contraction, MCMC_iter, current_est) %>%
       mutate(scenario = 
                case_when(scenario == "CCSM4_rcp45_m2085" ~ "RCP4.5 CCSM4",
                          scenario == "HadGEM2_ES_rcp45_m2085" ~ "RCP4.5 HadGEM2-ES",
                          scenario == "CCSM4_rcp85_m2085" ~ "RCP8.5 CCSM4",
                          scenario == "HadGEM2_ES_rcp85_m2085" ~ "RCP8.5 HadGEM2-ES")) %>%
       melt(id.vars = c("scenario", "MCMC_iter", "current_est")) %>%
       mutate(area_kkm2 = (value/100)/1000,
              current_est = (current_est/100)/1000,
              variable = factor(variable, 
                                levels = rev(c("expansion", "persistence", "contraction"))),
              scenario = factor(scenario,
                                levels = c("RCP4.5 CCSM4",
                                           "RCP4.5 HadGEM2-ES",
                                           "RCP8.5 CCSM4",
                                           "RCP8.5 HadGEM2-ES"))) %>%
       select(-value)
   }) %>% 
   bind_rows(.id = "spp") %>%
   mutate(spp = toupper(spp))) %>% 
  ggplot() + 
  geom_line(aes(x = MCMC_iter, y = current_est),
            size = .25,
            alpha = .5) +
  geom_bar(aes(x = MCMC_iter, y = area_kkm2, fill = variable),
           stat = "identity",
           alpha = .75) + 
  facet_grid(spp ~ scenario, scales = "free_y") +
  theme_bw() +
  scale_fill_viridis_d(direction = 1, option = "viridis") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  ylab("area (thousand km^2)") +
  xlab("MCMC sample")

# contraction/persistence/expansion as boxplots
(kde_future_diff %>%
  map(function(sp.x) {
    sp.x %>% 
      mutate(persistence = persistence/current_est * 100,
             expansion = expansion/current_est * 100,
             contraction = contraction/current_est * 100) %>%
      select(scenario, persistence, expansion, contraction) %>%
      melt(id.vars = c("scenario")) %>%
      mutate(variable = factor(variable, 
                               levels = c("expansion", "persistence", "contraction"))) %>%
      group_by(scenario, variable) %>%
      summarize(ymean = mean(value),
                q_05 = quantile(value, .05) %>% unname,
                q_95 = quantile(value, .95) %>% unname,
                ymin = min(value),
                ymax = max(value)) %>%
      mutate(scenario = 
               case_when(scenario == "CCSM4_rcp45_m2085" ~ "RCP4.5 CCSM4",
                         scenario == "HadGEM2_ES_rcp45_m2085" ~ "RCP4.5 HadGEM2-ES",
                         scenario == "CCSM4_rcp85_m2085" ~ "RCP8.5 CCSM4",
                         scenario == "HadGEM2_ES_rcp85_m2085" ~ "RCP8.5 HadGEM2-ES")) %>%
      mutate(scenario =factor(scenario,
                              levels = c("RCP4.5 CCSM4",
                                         "RCP4.5 HadGEM2-ES",
                                         "RCP8.5 CCSM4",
                                         "RCP8.5 HadGEM2-ES")))
    }) %>% 
  bind_rows(.id = "spp") %>%
  mutate(spp = toupper(spp))) %>%
  ggplot() +
  geom_boxplot(aes(ymin = ymin, ymax = ymax,
                   lower = q_05, upper = q_95,
                   middle = ymean,
                   x = variable, fill = scenario),
               stat = "identity",
               width = .5,
               alpha = .9,
               position = position_dodge(0.75)) +
  facet_grid(spp ~ ., scales = "free_y") +
  scale_fill_manual(values = unname(scen_colPal$hex_cb_gp[2:5]))  +
  ylab("area relative to current range-size (%)") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(plot.margin = margin(0, 0, 0, 0),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggtitle("Projected area relative to current range-size estimate")

# == BKDE v. Voronoi polygon estimates =================================== ####
# - load voronoi polygon estimates --------------------------------------- ####
range_size_key <- readRDS(paste0(path_output, "ch2-10-range_size_key.rds"))

# - make new legend labels & create base data set ------------------------ ####
# colors and old names...
legend_labels <- with(scen_colPal, hex_cb_gp %>% setNames(scenario))
future_lab_levels <- c("RCP 4.5 CCSM4",
                       "RCP 4.5 HadGEM2-ES",
                       "RCP 8.5 CCSM4",
                       "RCP 8.5 HadGEM2-ES")
names(legend_labels) <- case_when(
  names(legend_labels) == "current" ~ "current",
  names(legend_labels) == "CCSM4_rcp45_m2085" ~ future_lab_levels[1],
  names(legend_labels) == "HadGEM2_ES_rcp45_m2085" ~ future_lab_levels[2],
  names(legend_labels) == "CCSM4_rcp85_m2085" ~ future_lab_levels[3],
  names(legend_labels) == "HadGEM2_ES_rcp85_m2085" ~ future_lab_levels[4])

# - BKDE: change & rel-change (Figure 6) --------------------------------- ####
bkde_boxplot.df <- kde_future_diff %>%
  # create base data set
  imap(function(df.x, sp_name.x) {
    df.x %>% 
      rowwise %>%
      mutate(scenario = names(legend_labels)[grep(scenario, scen_colPal$scenario)]) %>% 
      mutate(scenario = factor(scenario, levels = future_lab_levels)) %>%
      ungroup %>%
      mutate_if(is.numeric, ~ (.x/100)/1000) %>% # convert from ha to 1k km^2
      # get change (thousand km^2) and relative-change (%)
      mutate(change = (future_est - current_est),
             relative_change = ((future_est - current_est)/current_est)*100) %>%
      select(scenario, change, relative_change) %>%
      melt(id.vars = "scenario") %>%
      group_by(scenario, variable) %>%
      summarize(q_05 = quantile(value, .05) %>% as.numeric,
                q_95 = quantile(value, .95) %>% as.numeric,
                ymean = mean(value),
                ymin = min(value),
                ymax = max(value)) %>%
      ungroup() %>%
      mutate(spp = recode(sp_name.x, !!!spp_common),
             spp = factor(spp, levels = unname(spp_common))) %>%
      mutate(yint = ifelse(sp_name.x == "abpr", NA, 0))
  }) %>%
  bind_rows()

p_body <- bkde_boxplot.df %>%
  split(.$variable) %>%
  imap(function(df.x, var_name.x) {
    # set up...
    std_margin <- margin(0, 0, 5.5, 0)
    legend_labels.w_return <- sub("^\\S+\\s+\\S+\\K\\s+", "\n", 
                                  names(legend_labels), perl=TRUE)
    
    # make the figure - depends on if it will be on the LHS or RHS...
    if (var_name.x == "change") {
      df.x$variable <- "Change in range-size"
      # figure on the LHS
      ggplot(df.x) +
        geom_boxplot(aes(ymin = ymin, ymax = ymax,
                         lower = q_05, upper = q_95,
                         middle = ymean,
                         x = scenario, fill = scenario),
                     position=position_dodge(width = 0),
                     stat = "identity",
                     width = .5,
                     alpha = .9) + 
        geom_hline(aes(yintercept = yint), color = "black") +
        facet_grid(spp ~ variable, scales = "free_y") +
        scale_fill_manual(values = legend_labels[2:5],
                          labels = legend_labels.w_return[2:5]) +
        ylab(expression(paste("thousand", ~km^2))) +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(plot.margin = margin(0, 0, 0, 0),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              strip.text.x = element_text(size = EJ_strip_text_size),
              strip.text.y = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text.y  = element_text(size = EJ_axis_text_size,
                                          color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size)) 
    } else {
      # figure on the RHS
      df.x$variable <- "Relative-change in range-size"
      ggplot(df.x) +
        geom_boxplot(aes(ymin = ymin, ymax = ymax,
                         lower = q_05, upper = q_95,
                         middle = ymean,
                         x = scenario, fill = scenario),
                     position=position_dodge(width = 0),
                     stat = "identity",
                     width = .5,
                     alpha = .9) + 
        geom_hline(aes(yintercept = yint), color = "black") +
        facet_grid(spp ~ variable, scales = "free_y") +
        scale_color_manual(values = legend_labels[2:5], guide = FALSE)  +
        scale_fill_manual(values = legend_labels[2:5],
                          labels = legend_labels.w_return[2:5]) +
        scale_y_continuous(position = "right") +
        ylab("percent") +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(plot.margin = margin(0, 0, 0, 0),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              strip.background.y = element_blank(),
              strip.text.y = element_blank(),
              strip.text.x = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text.y  = element_text(size = EJ_axis_text_size,
                                          color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size)) 
    }
  }) %>%
  reduce(., `+`)

# View it
p_body + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Change and relative-change in range-size",
                  subtitle = "Based on 2k MCMC samples & BKDE method") & 
  theme(legend.position = "bottom")

# save it
ggsave(filename = paste0(path_images, "Fig6.tiff"),
       plot = p_body + 
         plot_layout(guides = "collect") & 
         theme(legend.position = "bottom"),
       height = 4.6 * (EJ_dim_2col_width_max / 5.5),
       width = EJ_dim_2col_width_max,
       units = EJ_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")

  
# - VP: change & rel-change (Appendix S1: Figure S1) --------------------- ####
vp_boxplot.df <- range_size_key %>%
  # create base data set
  imap(function(df.x, sp_name.x) {
    df.x %>%
      filter(map_group == "maps_1k") %>%
      # add new scenario legend titles for plotting 
      rowwise %>%
      mutate(scenario = names(legend_labels)[grep(scenario, scen_colPal$scenario)]) %>% 
      ungroup %>%
      # summarize the difference in range size estimate by MCMC_iter...
      split(.$MCMC_iter) %>%
      map(function(x) {
        current_est <- x[x$scenario == "current", ]$range_size_ha
        current_est <- (current_est/100)/1000 # convert from ha to 1k km^2
        x %>% 
          filter(scenario != "current") %>%
          # convert from ha to 1k km^2
          mutate(range_size_ha = (range_size_ha/100)/1000) %>%
          # get change (thousand km^2) and relative-change (%)
          mutate(change = (range_size_ha - current_est),
                 relative_change = ((range_size_ha - current_est)/current_est)*100) %>%
          select(scenario, change, relative_change)
      }) %>%
      bind_rows() %>%
      melt(id.vars = "scenario") %>%
      group_by(scenario, variable) %>%
      summarize(q_05 = quantile(value, .05) %>% as.numeric,
                q_95 = quantile(value, .95) %>% as.numeric,
                ymean = mean(value),
                ymin = min(value),
                ymax = max(value)) %>%
      ungroup() %>%
      mutate(scenario = factor(scenario, levels = future_lab_levels),
             spp = recode(sp_name.x, !!!spp_common),
             spp = factor(spp, levels = unname(spp_common)),
             yint = ifelse(sp_name.x == "abpr", NA, 0))
  }) %>%
  bind_rows()

p_body <- vp_boxplot.df %>%
  split(.$variable) %>%
  imap(function(df.x, var_name.x) {
    # set up...
    std_margin <- margin(0, 0, 5.5, 0)
    legend_labels.w_return <- sub("^\\S+\\s+\\S+\\K\\s+", "\n", 
                                  names(legend_labels), perl=TRUE)
    
    # make the figure - depends on if it will be on the LHS or RHS...
    if (var_name.x == "change") {
      df.x$variable <- "Change in range-size"
      # figure on the LHS
      ggplot(df.x) +
        geom_boxplot(aes(ymin = ymin, ymax = ymax,
                         lower = q_05, upper = q_95,
                         middle = ymean,
                         x = scenario, fill = scenario),
                     position=position_dodge(width = 0),
                     stat = "identity",
                     width = .5,
                     alpha = .9) + 
        geom_hline(aes(yintercept = yint), color = "black") +
        facet_grid(spp ~ variable, scales = "free_y") +
        scale_fill_manual(values = legend_labels[2:5],
                          labels = legend_labels.w_return[2:5]) +
        ylab(expression(paste("thousand", ~km^2))) +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(plot.margin = margin(0, 0, 0, 0),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              strip.text.x = element_text(size = EJ_strip_text_size),
              strip.text.y = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text.y  = element_text(size = EJ_axis_text_size,
                                          color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size)) 
    } else {
      # figure on the RHS
      df.x$variable <- "Relative-change in range-size"
      ggplot(df.x) +
        geom_boxplot(aes(ymin = ymin, ymax = ymax,
                         lower = q_05, upper = q_95,
                         middle = ymean,
                         x = scenario, fill = scenario),
                     position=position_dodge(width = 0),
                     stat = "identity",
                     width = .5,
                     alpha = .9) + 
        geom_hline(aes(yintercept = yint), color = "black") +
        facet_grid(spp ~ variable, scales = "free_y") +
        scale_color_manual(values = legend_labels[2:5], guide = FALSE)  +
        scale_fill_manual(values = legend_labels[2:5],
                          labels = legend_labels.w_return[2:5]) +
        scale_y_continuous(position = "right") +
        ylab("percent") +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(plot.margin = margin(0, 0, 0, 0),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              strip.background.y= element_blank(),
              strip.text.y = element_blank(),
              strip.text.x = element_text(size = EJ_strip_text_size),
              axis.title.y = element_text(size = EJ_axis_title_size),
              axis.text.y  = element_text(size = EJ_axis_text_size,
                                          color = "black"),
              legend.text  = element_text(size = EJ_strip_text_size)) 
    }
  }) %>%
  reduce(., `+`)

# view it!
p_body + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Change and relative-change in range-size",
                  subtitle = "Based on 2k MCMC samples & VP method") & 
  theme(legend.position = "bottom")

# save it
ggsave(filename = paste0(path_images, "FigS1.tiff"),
       plot = p_body + 
         plot_layout(guides = "collect") &
         theme(legend.position = "bottom"),
       height = 4.6 * (EJ_app_dim_width_max/5.5),
       width = EJ_app_dim_width_max,
       units = EJ_app_dim_units, 
       dpi = EJ_dpi,
       compression = "lzw")

# - Range-size estimates (Appendix S1: Table S4) ------------------------- ####
list(range_size_key, kde_future_diff) %>%
  pmap(function(vp.x, bkde.x, round_digits = 3) {
    quant_range <- c(0.05,0.95)
    vp.x <- vp.x %>%
      filter(map_group == "maps_1k") %>%
      select(-map_group) %>%
      mutate(range_ksqkm = (range_size_ha/100)/1000) %>% # convert from ha to 1k km^2
      select(-range_size_ha) %>%
      group_by(scenario) %>%
      summarize(ql = quantile(range_ksqkm, min(quant_range)) %>% unname,
                qu = quantile(range_ksqkm, max(quant_range)) %>% unname,
                m  = mean(range_ksqkm))  %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      group_by(scenario) %>%
      summarize(vp = paste0(m, " (", ql, ", ", qu, ")"))
    
    bkde.x <- bkde.x %>% 
      select(MCMC_iter, current_est) %>% 
      distinct() %>%
      mutate(scenario = "current",
             range_ksqkm = (current_est/100)/1000) %>%
      select(scenario, range_ksqkm) %>%
      rbind.data.frame(bkde.x %>%
                         mutate(range_ksqkm = (future_est/100)/1000) %>%
                         select(scenario, range_ksqkm)) %>%
      group_by(scenario) %>%
      summarize(ql = quantile(range_ksqkm, min(quant_range)) %>% unname,
                qu = quantile(range_ksqkm, max(quant_range)) %>% unname,
                m  = mean(range_ksqkm))  %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      group_by(scenario) %>%
      summarize(bkde = paste0(m, " (", ql, ", ", qu, ")"))
    
    inner_join(vp.x, bkde.x, by = "scenario") %>%
      mutate(scenario = factor(scenario, levels = scen_colPal$scenario)) %>%
      arrange(scenario)
  })

# - change & relative-change in range-size (Appendix S1: Table S3) ------- ####
# 90% equal-tailed credible intervals for expansion (+), contraction (-), 
# or neither (include zero)

# Change: (future - current)
diff_90IC <- list(range_size_key, kde_future_diff) %>%
  pmap(function(vp.x, bkde.x, round_digits = 3) {
    vp.x <- vp.x %>%
      filter(map_group == "maps_1k") %>%
      split(.$MCMC_iter) %>%
      map(function(x) {
        current_est <- x[x$scenario == "current", ]$range_size_ha
        x %>% 
          filter(scenario != "current") %>%
          mutate(range_diff_ha = range_size_ha - current_est) %>%
          mutate(range_diff_ksqkm = (range_diff_ha/100)/1000) %>%
          select(scenario, range_diff_ksqkm, MCMC_iter)
      }) %>%
      bind_rows() %>%
      group_by(scenario) %>%
      summarize(ql = quantile(range_diff_ksqkm, 0.05) %>% unname,
                qu = quantile(range_diff_ksqkm, 0.95) %>% unname,
                m  = mean(range_diff_ksqkm))  %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      group_by(scenario) %>%
      summarize(vp = paste0(m, " (", ql, ", ", qu, ")"))
    
    bkde.x <- bkde.x %>% 
      mutate(diff_ksqkm = ((future_est - current_est)/100)/1000) %>%
      group_by(scenario) %>%
      summarize(ql = quantile(diff_ksqkm, 0.05) %>% unname,
                qu = quantile(diff_ksqkm, 0.95) %>% unname,
                m  = mean(diff_ksqkm))  %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      group_by(scenario) %>%
      summarize(bkde = paste0(m, " (", ql, ", ", qu, ")"))
    
    inner_join(vp.x, bkde.x, by = "scenario") %>%
      mutate(scenario = factor(scenario, levels = scen_colPal$scenario)) %>%
      arrange(scenario)
  })

# relative-change: (future - current)/ current
pc_90CI <- list(range_size_key, kde_future_diff) %>%
  pmap(function(vp.x, bkde.x, round_digits = 3) {
    vp.x <- vp.x %>%
      filter(map_group == "maps_1k") %>%
      split(.$MCMC_iter) %>%
      map(function(x) {
        current_est <- x[x$scenario == "current", ]$range_size_ha
        x %>% 
          filter(scenario != "current") %>%
          mutate(range_pc = ((range_size_ha - current_est)/current_est)*100) %>%
          select(scenario, range_pc, MCMC_iter)
      }) %>%
      bind_rows() %>%
      group_by(scenario) %>%
      summarize(ql = quantile(range_pc, 0.05) %>% unname,
                qu = quantile(range_pc, 0.95) %>% unname,
                m  = mean(range_pc))  %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      group_by(scenario) %>%
      summarize(vp = paste0(m, " (", ql, ", ", qu, ")"))
    
    bkde.x <- bkde.x %>% 
      mutate(range_pc = ((future_est - current_est)/current_est)*100) %>%
      group_by(scenario) %>%
      summarize(ql = quantile(range_pc, 0.05) %>% unname,
                qu = quantile(range_pc, 0.95) %>% unname,
                m  = mean(range_pc))  %>%
      mutate(across(where(is.numeric), ~ .x %>% round(round_digits))) %>%
      group_by(scenario) %>%
      summarize(bkde = paste0(m, " (", ql, ", ", qu, ")"))
    
    inner_join(vp.x, bkde.x, by = "scenario") %>%
      mutate(scenario = factor(scenario, levels = scen_colPal$scenario)) %>%
      arrange(scenario)
  })

# - Alternative figure for interpreting results -------------------------- ####
# VP & BKDE comparison figure (change in range-size)
p_diff_vp <- range_size_key %>%
  imap(function(sp.x, sp_name) {
    # new scenario legend titles for plotting
    legend_values <- with(scen_colPal, hex_cb_gp %>% setNames(scenario))
    legend_labels <- legend_values %>%
      names %>%
      gsub("_m2085", "", .) %>% 
      gsub("_rcp", " rcp", .) %>% 
      str_split(., " ") %>% 
      lapply(function(x) paste0(x[2], " ", x[1])) %>% 
      unlist %>% gsub("NA ", "", .)
    names(legend_values) <- legend_labels
    
    sp.x <- sp.x %>%
      filter(map_group == "maps_1k") %>%
      # add new scenario legend titles for plotting
      rowwise %>%
      mutate(scenario = legend_labels[grep(scenario, scen_colPal$scenario)]) %>% 
      mutate(scenario = factor(scenario, levels = legend_labels)) %>%
      ungroup %>%
      # summarize the difference in range size estimate by MCMC_iter...
      split(.$MCMC_iter) %>%
      map(function(x) {
        current_est <- x[x$scenario == "current", ]$range_size_ha
        x %>% 
          filter(scenario != "current") %>%
          mutate(range_diff_ha = range_size_ha - current_est) %>%
          mutate(range_diff_ksqkm = (range_diff_ha/100)/1000) %>%
          select(scenario, range_diff_ksqkm, MCMC_iter) %>%
          mutate(spp = toupper(sp_name))
      }) %>%
      bind_rows() %>%
      mutate(vp_bkde = "Voronoi polygon")
    
    legend_labels <- legend_labels[2:5]
    legend_values <- legend_values[2:5]
    
    # make the plots
    if (sp_name == "quke") {
      ggplot(sp.x %>% filter(spp != "abpr")) +
        geom_density(aes(x = range_diff_ksqkm, 
                         y = ..density.., 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ ., scales = "free_y") +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_fill_manual(values = legend_values,
                          labels = gsub(" ", "\n", legend_labels))  +
        xlab("change") +
        ylab("density") +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(plot.margin = margin(0, 0, 5.5, 0),
              legend.title = element_blank())
    } else if (sp_name == "abpr") {
      ggplot(sp.x) +
        geom_density(aes(x = range_diff_ksqkm, 
                         y = ..density.., 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ vp_bkde, scales = "free_y") +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_fill_manual(values = legend_values,
                          labels = gsub(" ", "\n", legend_labels))  +
        xlab("change") +
        ylab("density") +
        theme_bw() +
        theme(plot.margin = margin(0, 0, 5.5, 0),
              axis.title.x = element_blank()) + 
        guides(color = "none", fill = "none")
    } else {
      ggplot(sp.x) +
        geom_density(aes(x = range_diff_ksqkm, 
                         y = ..density.., 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(spp ~ ., scales = "free_y") +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_fill_manual(values = legend_values,
                          labels = gsub(" ", "\n", legend_labels))  +
        xlab("change") +
        ylab("density") +
        theme_bw() +
        theme(plot.margin = margin(0, 0, 5.5, 0),
              axis.title.x = element_blank()) + 
        guides(color = "none", fill = "none")
    }
  }) %>%
  map(~.x + geom_vline(aes(xintercept = 0), color = "black"))

p_diff_bkde <- kde_future_diff %>%
  imap(function(sp.x, sp_name) {
    # new scenario legend titles for plotting
    legend_values <- with(scen_colPal, hex_cb_gp %>% setNames(scenario))
    legend_labels <- legend_values %>%
      names %>%
      gsub("_m2085", "", .) %>% 
      gsub("_rcp", " rcp", .) %>% 
      str_split(., " ") %>% 
      lapply(function(x) paste0(x[2], " ", x[1])) %>% 
      unlist %>% gsub("NA ", "", .)
    names(legend_values) <- legend_labels
    
    sp.x <- sp.x %>%
      # add new scenario legend titles for plotting
      rowwise %>%
      mutate(scenario = legend_labels[grep(scenario, scen_colPal$scenario)]) %>% 
      mutate(scenario = factor(scenario, levels = legend_labels)) %>%
      ungroup  %>%
      # get difference and convert to thousand km^2
      mutate(range_diff_ksqkm = ((future_est - current_est)/100)/1000) %>%
      select(scenario, range_diff_ksqkm, MCMC_iter) %>%
      mutate(spp = toupper(sp_name)) %>%
      mutate(vp_bkde = "BKDE")

    legend_labels <- legend_labels[2:5]
    legend_values <- legend_values[2:5]
    
    # make the plots
    if (sp_name == "quke") {
      ggplot(sp.x) +
        geom_density(aes(x = range_diff_ksqkm, 
                         y = ..density.., 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_fill_manual(values = legend_values,
                          labels = gsub(" ", "\n", legend_labels))  +
        scale_y_continuous(position = "right") +
        xlab("change") +
        ylab("density") +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(plot.margin = margin(0, 0, 5.5, 0),
              legend.title = element_blank())
    } else if (sp_name == "abpr") {
      ggplot(sp.x) +
        geom_density(aes(x = range_diff_ksqkm, 
                         y = ..density.., 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        facet_grid(. ~ vp_bkde, scales = "free_y") +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_fill_manual(values = legend_values,
                          labels = gsub(" ", "\n", legend_labels))  +
        xlab("change") +
        ylab("density") +
        theme_bw() +
        theme(plot.margin = margin(0, 0, 5.5, 0),
              axis.title.x = element_blank()) + 
        scale_y_continuous(position = "right") +
        guides(color = "none", fill = "none")
    } else {
      ggplot(sp.x) +
        geom_density(aes(x = range_diff_ksqkm, 
                         y = ..density.., 
                         fill = scenario),
                     alpha = .5,
                     trim = TRUE) +
        scale_color_manual(values = legend_values, guide = FALSE) +
        scale_fill_manual(values = legend_values,
                          labels = gsub(" ", "\n", legend_labels))  +
        xlab("change") +
        ylab("density") +
        scale_y_continuous(position = "right") +
        theme_bw() +
        theme(plot.margin = margin(0, 0, 5.5, 0),
              axis.title.x = element_blank()) + 
        guides(color = "none", fill = "none")
    }
  }) %>%
  map(~.x + geom_vline(aes(xintercept = 0), color = "black"))

p_body <- list(p_diff_vp, p_diff_bkde, names(p_diff_vp)) %>%
  pmap(function(vp.x, bkde.x, sp_name.x) {
    # do differently for quke
    # - quke will be the bottom figure, so the only one that needs x title 
    # - we'll also harvest the legend from this species figure
    if (sp_name.x == "quke") {
      # add some space to the left of the legend
      vp.x <- vp.x + 
        theme(legend.position = "bottom",
              legend.box.margin = margin(3, 0, 0, 0),
              legend.title = element_blank())
      
      list(vp.x = vp.x, 
           bkde.x = bkde.x  + 
             guides(color = "none", fill = "none"))
      
    } else {
      list(vp.x = vp.x + 
             theme(axis.title.x = element_blank()) + 
             guides(color = "none", fill = "none"), 
           bkde.x = bkde.x + 
             theme(axis.title.x = element_blank()) + 
             guides(color = "none"))
    }
  }) %>%
  map(~ .x$vp.x + .x$bkde.x) %>%
  reduce(., `/`)

p_body + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Posterior densities for range size change",
                  subtitle = "Based on 2k MCMC samples; change is (thousand km^2)",
                  caption = "Note: (+) for future expansion and (-) for future contraction of total range size") & 
  theme(legend.position = "bottom")