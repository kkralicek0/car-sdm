###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 1-load-manip-raw-data.R
## Updated - 09-02-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Climate change induced shifts in suitable
## habitat projected for PNW tree species with spatial-Bayesian models' by 
## Karin Kralicek, Jay M. Ver Hoef, Tara M. Barrett, and Temesgen Hailemariam.
## 
## About - this script:
## - loads and formats data to be used in subsequent scripts, including data
##   from (prev. queried) FIA plot data & associated Norm81m data, spp codes,
##   hex-membership of plots, and ES.
## 
## About - output:
## - ch2-1-all-vars.Rdata
##   - location: path_data_TRUEloc
##   - nearly cleaned data(?)
## - ch2-1-elev-lat-plot-key.rds
##   - location: path_data_TRUEloc
##   - plot key with elevation & latitude (this for script 10-...R)
###############################################################################
## Notes & Reminders to self:::                                            ----
## Note - on floats and equality
##  - beware of [((i = .1 + .05) == 0.15) FALSE]; instead for floats use
##    [isTRUE(all.equal(i, 0.15)) TRUE]
##  - all.equal() compares two objects using a numeric tolerance of
##    .Machine$double.eps ^ 0.5; for greater tol need to look deeper
##  - Should be ok with PLT_CN since 14 char integer (not float), but! issues
##    exist when comparing Lat/Lon -- one solution when merging is to change
##    lat/lon into character vector for merge.
##
## Note - NAs "not allowed in subscripted assignments of data frames"
## - Suppose you want to do something like: 
##     df[df$col1 == 'a', ]$col2 <- 'b'
## - Then the above will throw an error if col1 has any NAs, so you actually
##   need to do:
##     df[df$col1 == 'a' & !is.na(df$col1), ]$col2 <- 'b'
## - You can see how R sees this by:
##     nrow(df %>% filter(col1 == 'a'))
##     nrow(df[df$col1 == 'a' & !is.na(df$col1), ])
##     nrow(df[df$col1 == 'a', ])
##   the first two will return the same nrow, but the last will give this number
##   PLUS then number of rows with NAs
###############################################################################
library(dplyr)
library(purrr)
library(stringr)
library(readxl)
library(sf)
library(ggplot2)
library(reshape2)
library(tidyr)

# == Paths =============================================================== ####
# IN - paths to data this script requires
# - note for path_hex: read_sf only works w/o trailing '/' at end of path...

path_USB <- "G:/"
path_top <- paste0(path_USB, "Karin phd thesis files/")
path_top_raw    <- paste0(path_top, "r - raw master data/")
path_top_proces <- paste0(path_top, "r - processed master data/")

path_FIAclim    <- paste0(path_top_raw, "the data - true version/")
path_keepPLT    <- paste0(path_top_raw, "key to cleaned data - true version/")
path_hex        <- paste0(path_top_raw, "FIA hex layer data - true version")
path_psmem_mask <- paste0(path_top_proces, 
                          "spp clean-buffer masks/processed data/")
path_ES         <- paste0(path_top_proces, 
                          "ECOMAP Ecological Sections/processed data/")
path_ESname     <- paste0(path_top_proces,
                          "ECOMAP Ecological Sections/",
                          "raw data - ecological sections shp/")

# OUT - destination of the true location .Rdata output from this script
path_data_TRUEloc <- paste0(path_top, "ch2 - code output/data/")


# == Load data: sppcd & specify spp of interest ========================== ####
# Specify spp of interest
keep_spp <- c("abpr", "psme", "qudo", "quga4", "quke")

# Load spp code key translating 4-letter abrev to numeric codes
sppcd <- read_xlsx(paste0(path_FIAclim, "Key_to_ThePublicData.xlsx"),
                   sheet = "SpeciesKey")
sppcd <- sppcd %>% rename(spp_abrev = SPECIES_SYMBOL) %>%
  filter(spp_abrev %in% str_to_upper(keep_spp)) %>%
  mutate(spp_abrev = str_to_lower(spp_abrev))

# Clean-up space
rm(keep_spp)


# == Load data: FIA (Y) & create short plot_id var ======================= ####
# True location FIA spp P/A and other plot data
# - (different format than public version of this file)
# - PRISM data in this file is based on Norm81m dataset

# load data - will give warnings about coercion... okay though...
FIA <- read_xlsx(paste0(path_FIAclim, "TheData_True.xlsx"),
                 col_types = c("text",             # PLT_CN
                               rep("skip", 8),     # PLT info: INVYR, etc.
                               rep("numeric", 2),  # Lat, Lon
                               rep("skip", 5),     # Elev/Aspect/Slope
                               rep("numeric", 70),    # Norm81m vars
                               rep("skip", 2),     # TotBAm, TotTPH
                               rep("numeric", 27), # bam-spp
                               rep("skip", 27)))   # tph-spp
FIA$PLT_CN <- as.factor(FIA$PLT_CN)

# Filter cols to spp of interest & rename cols
keep_cols <- FIA %>% select(-starts_with('bam')) %>% names
keep_cols <- c(keep_cols, paste0("bam", sppcd$SPCD))
FIA <- FIA[names(FIA) %in% keep_cols]
names(FIA)[(ncol(FIA) - nrow(sppcd) + 1):ncol(FIA)] <- sppcd$spp_abrev

# Create chr_LAT/chr_LON vars
# - these variables will be used to subset other datasets later on...
FIA <- FIA %>% mutate(chr_LAT = as.character(LAT_ACTUAL), 
                      chr_LON = as.character(LON_ACTUAL))

# Create plot_id variable (shorter version of PLT_CN)
FIA$plot_id <- 1:nrow(FIA)
FIA <- FIA %>% select(plot_id, everything()) # re order so plot_id 1st var.


# == Load data: hex & to FIA ============================================= ####
# Extract centroids to define CAR neighbor relations
# - here, grabbing centroids and associated plots with hex-polygons 
# - we'll use these centroids to determine neighbor relations
# - at the end everything returned in lat/lon (NAD83)

# read in shp file as sf object (tibble)
hex <- read_sf(dsn = path_hex, layer = "pnw_p2_hex")

# create a points collection from the FIA data
fia_sf <- st_as_sf(FIA %>% select(plot_id, LAT_ACTUAL, LON_ACTUAL),
                   coords = c("LON_ACTUAL", "LAT_ACTUAL"),
                   crs = as.numeric(st_crs(hex)$input))

# transform lat/lon to planar
# - sf assumes planar objects for these operations (intersect and centroids)
# - planar projection used here is UTM Zone 11N (EPSG: 32611)
hex_pl <- st_transform(hex, crs = 32611)     # hex-polygons
fia_pl <- st_transform(fia_sf, crs = 32611)  # point data

# intersect hex-polygons with point collection of FIA plot-level data
# - specifying `sparse = FALSE` returns a matrix where rows correspond to arg1
#   (hex_pl) and cols to arg2 (fia_pl)
hex_walk <- st_intersection(hex_pl, fia_pl)

# associate FIA plots w. FIAHEX_ID
st_geometry(hex_walk) <- NULL # remove geometry to turn back into a regular df
hex_walk <- hex_walk %>% select(FIAHEX_ID, plot_id)

# get hex centroids for these FIA plots
# - will throw a warning, but it is okay (just says now associating area-attr
#   with point data, so be careful)
centroid_pl <- st_centroid(hex_pl) %>% 
  select(FIAHEX_ID, geometry) %>%
  filter(FIAHEX_ID %in% hex_walk$FIAHEX_ID)

# associate FIA plots w. centroid locations
centroid_df <- st_coordinates(centroid_pl) %>% as.data.frame
centroid_df$FIAHEX_ID <- centroid_pl$FIAHEX_ID
hex_walk <- left_join(hex_walk, centroid_df, 
                      by = "FIAHEX_ID") %>%
  rename(hex_centroid.y_UTM11N_ACTUAL = Y,
         hex_centroid.x_UTM11N_ACTUAL = X,
         fiahex_id = FIAHEX_ID)

# merge with FIA
FIA <- left_join(FIA, hex_walk, by = "plot_id")

# clean-up space
rm(hex, fia_sf, hex_pl, fia_pl, 
   centroid_pl, centroid_df,
   hex_walk)

# at this point, how many hexagons have multiple plots within them?
FIA %>% 
  group_by(fiahex_id) %>% 
  summarise(n.plots = n()) %>% 
  filter(n.plots > 1) %>% 
  ungroup %>% 
  group_by(n.plots) %>% 
  summarize(n.hex = n()) %>% 
  ungroup

# == Load data: nonfor & to FIA ========================================== ####
# Nonforest reason code for excluding problematic nonforest plots.
# - correct PLT_CN forest/fulforst coding (raw true-coord data wrong)
# - TB-created 'NonForNat' variable; indicates PRESNFCD cond where trees could 
#   still naturally grow in theory. Kept PRESNFCDs in our data set, are:
#   (10) Ag land; (12) Pasture; (13) idle farmland; (20) rangeland; (40) other;
#   (41) nonvegetated; (42) wetland; (45) nonforest-chaparral
# - for CAR (z = (z_o, z_u)), create plot_obsv and fiahex_obsv vars...

# Load data
# - will give warnings about coercion... okay though...
nonfor <- read_xlsx(paste0(path_FIAclim, "TheData_NonForestReason.xlsx"),
                    col_types = c("text",              # PLT_CN
                                  rep("numeric", 5)))  # condition status vars
nonfor$PLT_CN <- as.factor(nonfor$PLT_CN)

# Create a keep_plot variable
nonfor$keep_plot <- nonfor$NonForNat

# - fix typo: there isn't a PRESNFCD = 0, so fix and drop these plots 
typo_PRESNFCD_0 <- (!is.na(nonfor$COND_PRESNFCD)) & (nonfor$COND_PRESNFCD == 0)
sum(typo_PRESNFCD_0) # number of such plots

nonfor[typo_PRESNFCD_0, ]$COND_PRESNFCD <- NA
nonfor[typo_PRESNFCD_0, ]$keep_plot <- 0

# - drop un-sampled plots from sample
#   note: some nonforested plots are unsampled but PRESNFCDs are identified 
#   via areal photos etc.
sum(is.na(nonfor$Sampled)) # number of such plots total
sum(is.na((nonfor %>% filter(NonForNat == 1))$Sampled)) # unsampled nonfor w. keep presnfcd

nonfor[is.na(nonfor$Sampled), ]$keep_plot <- 0

# Peak at this data structure
nonfor %>% 
  group_by(NonForNat, Sampled, Forest, FulForst) %>% 
  summarize (n_plts = n())

with(nonfor, table(keep_plot, COND_PRESNFCD, useNA = "always"))

# Create 'plot_obsv' variable from 'Forest'
# - this indicates if the plot is forested or nonforested 
#   (e.g. if the plot has an observed or unobserved response)
nonfor <- nonfor %>% mutate(plot_obsv = ifelse(is.na(Forest), 0, Forest))

# Filter FIA plots to keep_plots
FIA <- FIA %>% inner_join(nonfor %>% 
                            filter(keep_plot == 1) %>% 
                            select(PLT_CN, plot_obsv),
                          by = "PLT_CN")

# Create 'fiahex_obsv' variable from 'Forest'
# - this indicates whether the hex is composed of entirely observed plots 
#   (e.g. forested), unobserved plots (e.g. nonforested), or if it is a mix of
#   obsv and unobserved
FIA <- FIA %>%
  group_by(fiahex_id) %>%
  mutate(n_plots_hex = n(),
         n_plots_obsv = sum(plot_obsv),
         n_plots_unobsv = sum(abs(plot_obsv - 1))) %>%
  mutate(fiahex_obsv = ifelse(n_plots_obsv > 0 & n_plots_unobsv > 0,
                              "mix",
                              ifelse(n_plots_obsv > 0,
                                     "obsv",
                                     "unobsv"))) %>%
  select(-n_plots_obsv, -n_plots_unobsv) %>%
  ungroup

# Clean-up space
rm(nonfor, typo_PRESNFCD_0)

# == Load data: ES (re) & to FIA; Create ES_key ========================== ####
# Load ES crosswalk to PLT_CN
# - ECOMAP Ecological Sections are subregions
# - adding these now to potentially use as random effects later or subset spp
# - will give warnings about coercion... okay though...
ES_crosswalk <- read_xlsx(paste0(path_ES, "EcoSections.xlsx"),
                          col_types = c("text",           # PLTCN
                                        rep("numeric",2), # lat/lon
                                        "text"))          # ES code
ES_crosswalk$PLT_CN <- as.factor(ES_crosswalk$PLT_CN)

# Add ES variable to FIA df
FIA <- FIA %>% inner_join(ES_crosswalk %>% select(PLT_CN, ES), 
                          by = "PLT_CN")

# Create ES_key
ES_key <- read_xls(paste0(path_ESname, "ES descriptions for OR WA CA.xls"),
                   sheet = "ES - OR WA CA",
                   col_types = c(rep("skip", 7), # higher level info
                                 rep("text", 2), # ES code and ES name
                                 "skip"))        # in depth description

# Clean-up space
rm(ES_crosswalk)


# == Load data: distance buffer and cleaning keys ======================== ####
# - Load: Distance buffer mask ------------------------------------------- ####
# - drops plots outside a 200km+ radius from plots where spp P (False A)
# - Created by Bianca Eskelson (BE) to remove all absence-plots > 200km from 
#   a plot on which the spp was present. The intent is to include a measure 
#   of unsuitable habitat for modeling, but not an unnecessary number of plots.
# - Here '1' is keep and 'NA' is drop/mask
mask_distance <- read.csv(paste0(path_keepPLT, "DataToKarin.txt"),
                          colClasses = c("factor",                 # PLT_CN
                                         rep("numeric", 27)))      # spp 1/NA

# - Load: Cleaning key --------------------------------------------------- ####
# - removes potential spp mis-IDs (identified by BE in ArcMap)
# - The plots this key removes were suspected to be non-valid observations of 
#   a spp presence (e.g. data entry error, but as we could not check the FIA  
#   plot jackets, we remove these plots here just to be safe). 
# - These plots are noted in BE's SAS script as being removed to create the
#   mask_distance, but we still need to manually remove these plots here
#   (e.g. although these plots didn't influence the dist buffer for A-plots, 
#   they will still show up as '1' if within 200km of a valid presence plot)

# load data - will give warnings about coercion... okay though...
cleaning_key <- read_xlsx(
  paste0(path_keepPLT, "suspected nonvalid presence plots.xlsx"),
  col_types = c(rep("text", 2),     # spp, PLT_CN
                rep("skip", 2),     # PLOT_ACTUAL, STATECD
                rep("numeric", 2),  # Lat, Lon
                rep("text", 2)))    # who suspected, why
cleaning_key$PLT_CN <- as.factor(cleaning_key$PLT_CN)

# - Load: PSME variety mask (to identify PSMEM sample plots) ------------- ####
# - Two varieties of PSME in study area, PSMEM (coastal; v. menziesii) and 
#   PSMEG (Rocky Mountain (interior); v. glauca) -- we're focusing on PSMEM
# - Divisions between varieties based on work by Grugger et al. (2010)
# - This dataset masks out three things:
#   - plots in PSMEG (var glauca)'s region 
#   - non-valid presence plots (suspected FIA data entry errors)
#   - absence-plots > 200 km (geodesic) distance from valid PSMEM presence-plots

# Note -- to keep this script clean, we'll change spp codes at the end,
# but for now know that 'psme' once masked is 'psmem'
mask_psmem_clean <- read_xlsx(
  paste0(path_psmem_mask, "psmem_clean_mask.xlsx"),
  col_types = c("text", #           PLT_CN
                rep("numeric", 2), # lat/lon
                "numeric"))        # 0(drop)/1(keep)
mask_psmem_clean$PLT_CN <- as.factor(mask_psmem_clean$PLT_CN)
mask_psmem_clean <- mask_psmem_clean %>% rename(psme = psmem_keep)
mask_psmem_clean[mask_psmem_clean$psme == 0, ]$psme <- NA

# - Create keepPLT & keepPLT_spp ----------------------------------------- ####
# Add chr_LAT and chr_LON vars from FIA to keepPLT by PLT_CN
keepPLT <- mask_distance %>% 
  inner_join(FIA %>% select(plot_id, PLT_CN, chr_LAT, chr_LON), 
             by = c("PLT_CN"))

# Filter cols to spp of interest & rename cols
# - note: cols in alphabeta-(spp_abrev)-order, so okay to rename this way
keep_cols <- paste0("Plot_", sppcd$spp_abrev)
keep_cols <- c("plot_id", "PLT_CN", "chr_LAT", "chr_LON", keep_cols)
keepPLT   <- keepPLT %>% select(!!keep_cols)

names(keepPLT)[5:ncol(keepPLT)] <- str_sort(sppcd$spp_abrev)

# Drop old distance buffer and replace with new psme variety mask for PSMEM
keepPLT <- keepPLT %>% select(-psme) %>%
  inner_join(mask_psmem_clean %>% select(PLT_CN, psme),
             by = "PLT_CN") %>%
  select(plot_id, PLT_CN, chr_LAT, chr_LON, !!sppcd$spp_abrev)

# Change entries in keepPLT to 'NA' for spp x PLT_CN combos in cleaning_key
keepPLT[(keepPLT$PLT_CN %in% 
           (cleaning_key %>% filter(spp == "abpr"))$PLT_CN),]$abpr <- NA
keepPLT[(keepPLT$PLT_CN %in% 
           (cleaning_key %>% filter(spp == "psme"))$PLT_CN),]$psme <- NA
keepPLT[(keepPLT$PLT_CN %in% 
           (cleaning_key %>% filter(spp == "qudo"))$PLT_CN),]$qudo <- NA

# Drop any plots that are 'NA' for all the keep-spp
row_sum <- apply(keepPLT[5:ncol(keepPLT)], 1, sum, na.rm = TRUE)
keepPLT <- keepPLT[row_sum > 0, ]

# Create keepPLT_spp: a list of dfs by spp
# - each df is sp-specific and only has in-plots for that spp
keepPLT_spp <- sppcd$spp_abrev %>% 
  lapply(function(x) (keepPLT %>% select(PLT_CN, plot_id, x))) %>% 
  map(function(df) df[complete.cases(df), ])

# - Subset FIA to keepPLT ------------------------------------------------ ####
FIA <- FIA %>% inner_join(keepPLT %>% select(plot_id), by = "plot_id")

# Clean-up space
rm(keep_cols, row_sum, mask_distance, cleaning_key, mask_psmem_clean)

# == Create P/A vars from bam** in FIA =================================== ####
# Pre-queried dataset we're working with has TPA and BAM for each species on
# a plot... using BAM to filter to plots on which the spp was present (checked
# and as would be expected/hoped, same plots returned if TPA was used).

# Change all 'NA' to zero & '> 0' to one
FIA <- FIA %>% 
  mutate_at(sppcd$spp_abrev, ~ifelse(is.na(.), 0, .)) %>%
  mutate_at(sppcd$spp_abrev, ~ifelse(. > 0, 1, .))

# == Rename psme for what it is... psmem ================================= ####
FIA <- FIA %>% rename(psmem = psme)

# == Reorganize FIA ====================================================== ####
FIA <- FIA %>% 
  select(-chr_LAT, -chr_LON) %>%
  select(plot_id, fiahex_id, 
         LAT_ACTUAL, LON_ACTUAL,
         hex_centroid.x_UTM11N_ACTUAL, hex_centroid.y_UTM11N_ACTUAL,
         ES, plot_obsv, fiahex_obsv, n_plots_hex,
         everything())

# == Split FIA into sp-specific list of dfs ============================== ####
# Split FIA by sp into a list of dfs
spp_names <- c(sppcd$spp_abrev[c(1,3:5)], "psmem")[c(1,5,2:4)]
nonspp_cols <- FIA %>% select(-!!(spp_names)) %>% names
FIA_spp <- lapply(spp_names, function(x) {
  df.x <- FIA %>% 
    select(!!nonspp_cols, x) %>% 
    mutate(spp = x) %>%
    rename(p_a = !!x) %>%
    select(plot_id, spp, p_a, ES, everything())
})

names(FIA_spp) <- spp_names

# filter FIA to keep-plots for each spp
FIA_spp <- pmap(list(FIA_spp, keepPLT_spp), 
                function(df.x, df.y) df.x %>% filter(plot_id %in% df.y$plot_id))

# clean-up space
rm(nonspp_cols)

# == Create plot_key with lat lon and plot_id ============================ ####
plot_key <- FIA_spp %>% map(function(df) {
  df <- df %>% select(plot_id, PLT_CN, LAT_ACTUAL, LON_ACTUAL)
})

# == Export as Rdata ===================================================== ####
save(sppcd, spp_names, ES_key, plot_key,  # reference dfs
     FIA, FIA_spp,                        # data (Y, r.e., W, U)
     file = paste0(path_data_TRUEloc, "ch2-1-all-vars.Rdata"))

# == save a plot key with elevation & latitude for script 10..R ========== ####
# load data - will give warnings about coercion... okay though...
elev <- read_xlsx(paste0(path_FIAclim, "TheData_True.xlsx"),
                 col_types = c("text",            # PLT_CN
                               rep("skip", 10),   # PLT info, Lat, Lon
                               "numeric",         # Elev
                               rep("skip", 130))) # etc.
elev$PLT_CN <- as.factor(elev$PLT_CN)

tmp <- plot_key %>% 
  bind_rows %>%
  distinct(plot_id, .keep_all = TRUE) %>%
  left_join(elev, by = "PLT_CN") %>%
  select(-PLT_CN)

saveRDS(tmp,
        file = paste0(path_data_TRUEloc, "ch2-1-elev-lat-plot-key.rds"))




