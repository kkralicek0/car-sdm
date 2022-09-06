###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 8-voronoi-polygons.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Climate change induced shifts in suitable
## habitat projected for PNW tree species with spatial-Bayesian models' by 
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - creates spatial data.frames (sf) for use with plotting and area-based 
##   calculations for both plot-level and hexagon-level results.
## - Note on backward use:
##   - Script 6-fit-CAR-models.R uses both hex_polygon_*.sf objects to examine
##     organization of the spatial random effects (z-maps, both true & public)
##
## About - output:
## - ch2-8-plot_polygons_true.rds; ch2-8-plot_polygons_public.rds
##   - location: path_output
##   - .rds for plot_polygons_true.sf & plot_polygons_public.sf
##   - Voronoi polygons based on true plot locations (confidential) and the
##     public plot locations (fuzzed and swapped; only these can be used for
##     manuscript figures / sharing / etc.)
## - ch2-8-hex_polygons_true.rds; ch2-8-hex_polygons_public.rds
##   - location: path_output
##   - .rds for hex_polygons_true.sf & hex_polygons_public.sf
##   - hexagon layers based on true locations (confidential) and the
##     public locations (perturbed non-confidential version; only these can be
##     used for manuscript figures / sharing / etc.)
## - ch2-8-drop_plot_ids.rds
##   - location: path_output
##   - Plots to drop for prediction b/c not covered by CNA data set or b/c 
##     polygons created in this script occur entirely within ocean
###############################################################################
# data manipulation
library(dplyr)       # for ... E V E R Y T H I N G ...
library(purrr)       # working with lists

# plotting
library(ggplot2)

# spatial data (st_voronoi, etc.)
library(sf)       # for st_voronoi(), etc.

# == Paths =============================================================== ####
# path to Rdata
path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/ch2 - code output/")

path_data   <- paste0(path_top, "data/")
path_output <- paste0(path_top, "output/")
path_public_crosswalk <- paste0(path_USB, 
                                "Karin phd thesis files/r - raw master data/",
                                "crosswalk plots - true and public/")

# path to hex spatial data
path_hex <- paste0(path_USB, 
                   "Karin phd thesis files/r - raw master data/",
                   "FIA hex layer data - true version")
path_states <- paste0(path_data, "study area polygon")

# == Define projections (crs) ============================================ ####
# Primary planar projection:
# - [distance preserving] UTM Zone 11N (EPSG: 32611)
crs_UTM11N <- 32611

# Area planar projection:
# - [area preserving] Albers Equal Area Conic (EPSG: 5070)
# - only used temporarily for calculations, data projected back to primary
#   projection post.
crs_Albers <- 5070

# == Load data =========================================================== ####
# - Study area boundaries: load & project -------------------------------- ####
# state boundaries as sf
# - Note on data: OR/WA/CA boundary based on 2017 Cartographic Boundary 
#   of the US (US Census Bureau); clipped to OR/WA/CA extent by Ty Nietupski
#   (emailed on 8/11/2020)
study_area.sf <- read_sf(dsn = path_states, layer = "ORWACA_StateBoundaries") %>%
  st_transform(crs = crs_UTM11N)

# - Hex data: load & project (true & public) ----------------------------- ####
# hex polygons: true location (confidential)
hex_polygons_true.sf <- read_sf(dsn = path_hex, layer = "pnw_p2_hex") %>%
  st_transform(crs = crs_UTM11N)

# hex polygons: public location (non-confidential)
hex_polygons_public.sf <- read_sf(dsn = path_hex, layer = "FIA_pres_hex") %>%
  st_transform(crs = crs_UTM11N)

# - Plot data: load (true & public crosswalk & drop-ids) ----------------- ####
# True location
# - use: confidential, for manuscript summaries
# - version loaded is training data used to fit Ch2 models (formatted & std)
spp <- c("abpr", "psmem", "qudo", "quga4", "quke")
plot_true_spp <- lapply(spp, function(x) {
  readRDS(paste0(path_output, "ch2-5-sampler-ready-data_", x, ".rds"))
}) %>% setNames(spp)

# Public location crosswalk to true 
# (contains confidential info)
# - use: non-confidential, fuzzed & swapped, for manuscript maps
if (file.exists(paste0(path_public_crosswalk, 
                       "crosswalk_true_public_PLT_CN.rds"))) {
  public_crosswalk <- readRDS(paste0(path_public_crosswalk, 
                                    "crosswalk_true_public_PLT_CN.rds"))
} else {
  # if it doesn't already exist, make it!
  # plot data -- true & public
  # (will give warnings about coercion... okay though...)
  df_true <- read_xlsx(
    paste0(path_USB, "Karin phd thesis files/r - raw master data/",
           "the data - true version/",
           "TheData_True.xlsx"),
    col_types = c("text",             # PLT_CN
                  rep("skip", 8),     # ...
                  rep("numeric", 2),  # Lat, Lon
                  rep("skip", 131)))  # ...
  df_true$PLT_CN <- as.factor(df_true$PLT_CN)
  df_public <- read_xlsx(
    paste0(path_USB, "Karin phd thesis files/r - raw master data/",
           "the data - public version/", 
           "TheData_Public.xlsx"),
    col_types = c("text",             # PLT_CN
                  rep("skip", 3),     # ...
                  rep("numeric", 2),  # Lat, Lon
                  rep("skip", 130)))  # ...
  df_true$PLT_CN <- as.factor(df_true$PLT_CN)
  
  # make & save crosswalk
  public_crosswalk <- df_true %>% inner_join(df_public, by = "PLT_CN")
  saveRDS(public_crosswalk, 
          file = paste0(path_public_crosswalk, 
                        "crosswalk_true_public_PLT_CN.rds"))
  
  # clean-up space
  rm(df_true, df_public)
}

# Drop plot ids
# - use: plot_id's to drop for prediction
# - these plots were outside of ClimateNA's extent for future climate data
drop_plot_ids <- readRDS(file = paste0(path_output, "ch2-7-drop_plot_ids.rds"))

# == Drop non-pred & project ============================================= ####
# - Plots (true & public) ------------------------------------------------ ####
plot_points_true.sf <- plot_true_spp %>% map(function(plot.x) {
  plot.x$moddat %>% 
    # only keep data for matching to polygons (& p_a for code development)
    select(plot_id, fiahex_id, LAT_ACTUAL, LON_ACTUAL, p_a) %>% 
    # subset training data to prediction locations (i.e. no drop plots/hexs)
    filter(!plot_id %in% drop_plot_ids) %>%
    # convert to sf object for plot points
    st_as_sf(coords = c("LON_ACTUAL", "LAT_ACTUAL"),
             crs = "NAD83") %>%
    st_transform(crs = crs_UTM11N)
})

plot_points_public.sf <- plot_true_spp %>% map(function(plot.x) {
  plot.x$moddat %>%
    # only keep data for matching to polygons (& p_a for code development)
    select(plot_id, fiahex_id, p_a, PLT_CN) %>%
    # subset training data to prediction locations (i.e. no drop plots/hexs)
    filter(!plot_id %in% drop_plot_ids) %>%
    # add Public coordinate data
    left_join(public_crosswalk, by = "PLT_CN") %>%
    select(-PLT_CN, -LAT_ACTUAL, -LON_ACTUAL) %>%
    # convert to sf object for plot points
    st_as_sf(coords = c("LON_FUZZ", "LAT_FUZZ"),
             crs = "NAD83") %>%
    st_transform(crs = crs_UTM11N)
})

# - Hex (true) ----------------------------------------------------------- ####
# True location: hex_boundary
# - use: clipping voronoi polygons (based on true plot locations)
# - hex_boundary (multipolygon) by spp
hex_boundary_true.sf <- plot_points_true.sf %>% map(function(plot.x) {
  hex_polygons_true.sf %>% 
    filter(FIAHEX_ID %in% plot.x$fiahex_id) %>%
    st_union
})

# True location: hex_polygons
# - use: z-maps
# - list by spp of hex polygons
hex_polygons_true.sf <- plot_points_true.sf %>% map(function(plot.x) {
  hex_polygons_true.sf %>% 
    filter(FIAHEX_ID %in% plot.x$fiahex_id) %>%
    select(FIAHEX_ID, geometry) %>%
    rename(fiahex_id = FIAHEX_ID)
})

# - Hex (public) --------------------------------------------------------- ####
# Crosswalk true to public hex data

# Public location: hex_polygons
# - use: z-maps
# - list by spp of hex polygons
# - below, matching fiahex_id to public hex locations (perturbed) 
hex_polygons_public.sf <- hex_polygons_true.sf %>% map(function(hex_true.x) {
  hex_polygons_public.sf %>% 
    select(ROW3, PATH3) %>% 
    st_set_agr(., "constant") %>%
    st_join(hex_true.x %>% 
              st_set_agr(., "constant") %>%
              st_centroid,
            left = FALSE) %>%
    select(fiahex_id)
})

# Public location: hex_boundary
# - use: clipping voronoi polygons (based on public plot locations)
# - hex_boundary (multipolygon) by spp
hex_boundary_public.sf <- hex_polygons_public.sf %>% map(~.x %>% st_union)

# == Find Voronoi polygons + area ======================================== ####
# get (rough) Voronoi polygons & intersect w. study-area & hex boundaries
# - why intersect with these features?
#   - study_area: b/c some on coast have hex's that include water
#   - hex boundaries: b/c only want to comment on these units...
#     (we don't want to predict/comment on distribution in: developed areas
#     aka towns; areas where the spp will not grow (non-for for some...); or
#     areas that are outside study extent (200km from observed p-plot))

# - True locations (polygons + area) ------------------------------------- ####
# get Vornoi polygons
plot_polygons_true.sf <- list(plot_points_true.sf, hex_boundary_true.sf) %>% 
  pmap(function(point.x, hex.x) {
    point.x %>% 
      # get Voronoi tessellation 
      st_union %>%    # make MULTIPOINT geometry for st_voronoi (i.e. 1 row)
      st_voronoi %>%  # returns GEOMETRYCOLLECTION...
      st_cast %>%     # ... split geometries to get individual POLYGONs
      # hex boundary intersection
      st_intersection(., hex.x) %>%
      st_cast
  })

# add plot_id attribute to features
plot_polygons_true.sf <- list(plot_polygons_true.sf, plot_points_true.sf) %>% 
  pmap(function(ploy.x, point.x) {
    # add plot_id attribute to features & return object (out_sf)
    # - note: st_join by default does a left join, but in our case this is okay 
    #   b/c every point (plot) has an associated polygon
    ploy.x %>%
      data.frame(geometry = .) %>%
      st_sf(.) %>%
      st_join(., point.x %>% select(plot_id, p_a)) %>%
      st_set_agr(., "constant") # specify attributes as spatially constant
  })

# intersect with study_area boundaries
plot_polygons_true.sf <- plot_polygons_true.sf %>% map(function(poly.x) {
  poly.x %>%
    # state boundary intersection
    st_intersection(., study_area.sf) %>%
    st_cast
})

# add area attribute to features
plot_polygons_true.sf <- plot_polygons_true.sf %>% map(function(poly.x) {
  poly.x %>%
    st_transform(crs = crs_Albers) %>% # Albers Equal Area Conic (planar to preserve area)
    mutate(area_ha = st_area(.) %>% 
             units::set_units(., ha) %>% 
             as.numeric) %>%
    st_transform(crs = crs_UTM11N) %>% # back to same crs as hex & study_area..
    st_set_agr(., "constant") # specify attributes as spatially constant
})

# create non-sf df with area info
plot_polygons_true.df <- plot_polygons_true.sf %>% 
  map(~as.data.frame(.x) %>% select(plot_id, area_ha, p_a))

plot_polygons_true.df$qudo %>% 
  group_by(p_a) %>% 
  summarize(area_ha = sum(area_ha)) # check... del post

# - Public locations (polygons only, for plotting only) ------------------ ####
# get Vornoi polygons
plot_polygons_public.sf <- list(plot_points_public.sf, hex_boundary_public.sf) %>% 
  pmap(function(point.x, hex.x) {
    point.x %>% 
      # get Voronoi tessellation 
      st_union %>%    # make MULTIPOINT geometry for st_voronoi (i.e. 1 row)
      st_voronoi %>%  # returns GEOMETRYCOLLECTION...
      st_cast %>%     # ... split geometries to get individual POLYGONs
      # hex boundary intersection
      st_intersection(., hex.x) %>%
      st_cast
  })

# add plot_id attribute to features
plot_polygons_public.sf <- list(plot_polygons_public.sf, plot_points_public.sf) %>% 
  pmap(function(ploy.x, point.x) {
    # add plot_id attribute to features & return object (out_sf)
    # - note: st_join by default does a left join, but in our case this is okay 
    #   b/c every point (plot) has an associated polygon
    ploy.x %>%
      data.frame(geometry = .) %>%
      st_sf(.) %>%
      st_join(., point.x %>% select(plot_id, p_a)) %>%
      st_set_agr(., "constant") # specify attributes as spatially constant
  })

# intersect with study_area boundaries
plot_polygons_public.sf <- plot_polygons_public.sf %>% map(function(poly.x) {
  poly.x %>%
    # state boundary intersection
    st_intersection(., study_area.sf) %>%
    st_cast
})

# == Update & save drop_plot_ids ========================================= ####
# Drop plot ids
# - use: plot_id's to drop for prediction
# - add new plots dropped for having polygons that occur entirely within ocean
#   to those plots outside of ClimateNA's extent for future climate data
drop_plot_ids2 <- list(plot_points_true.sf, plot_polygons_true.sf) %>% 
  pmap(function(points.x, poly.x) {
    c((points.x %>% filter(!plot_id %in% poly.x$plot_id))$plot_id,
      drop_plot_ids)
  })

saveRDS(drop_plot_ids2,
        file = paste0(path_output, "ch2-8-drop_plot_ids.rds"))

# # == Save: polygon data ================================================== ####
# # True location / area files...
# saveRDS(study_area.sf,
#         file = paste0(path_output, "ch2-8-study_area_sf.rds"))
# saveRDS(hex_polygons_true.sf,
#         file = paste0(path_output, "ch2-8-hex_polygons_true_sf.rds"))
# saveRDS(plot_polygons_true.sf,
#         file = paste0(path_output, "ch2-8-plot_polygons_true_sf.rds"))
# saveRDS(plot_polygons_true.df,
#         file = paste0(path_output, "ch2-8-plot_polygons_true_df.rds"))
# saveRDS(plot_points_true.sf,
#         file = paste0(path_output, "ch2-8-plot_points_true_sf.rds"))
# 
# # Public locations (associated with true names - so confidential!)...
# saveRDS(hex_polygons_public.sf,
#         file = paste0(path_output, "ch2-8-hex_polygons_public_sf.rds"))
# saveRDS(plot_polygons_public.sf,
#         file = paste0(path_output, "ch2-8-plot_polygons_public_sf.rds"))
# saveRDS(plot_points_public.sf,
#         file = paste0(path_output, "ch2-8-plot_points_public_sf.rds"))
