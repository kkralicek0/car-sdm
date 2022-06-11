###############################################################################
## Project - CAR-SDM & climate change                                      ----
## Script  - 5-ready-data-rholookup-for-sampler.R
## Updated - 06-10-2022
## Author  - Karin Kralicek (karin.kralicek@usda.gov)
## 
## This R script is in support of 'Climate change induced shifts in suitable
## habitat for five tree species in the Pacific Northwest projected with 
## spatial-Bayesian hierarchical models' by Karin Kralicek, Jay Ver Hoef, 
## Tara Barrett, and Temesgen Hailemariam.
##
## About - this script:
## - format data for sampler and create rho look-up tables...
##
## About - output:
## - ch2-5-sampler-ready-data_*.rds; ch2-5-rholookup-table_*.rds
##   - location: path_output (USB)
##     - naming convention (*) == (spp abbreviation) 
##       (e.g. ch2-5-sampler-ready-data_abpr.Rdata, etc.)
##   - MCMC sampler ready data and rho look-up tables for each species
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

# path to functions script
path_Ndrive <- paste0("C:/Users/karin/Desktop/phd thesis - sample/",
                      "Ch2 - CAR prediction for species/code/")
path_fn   <- paste0(path_Ndrive, "scripts - real data/")

# == Load fn ============================================================= ####
load(paste0(path_data, "ch2-2-data.Rdata")) # FIA, FIA_spp, etc.
load(paste0(path_data, "ch2-4-modeling-data.Rdata")) # subset modeling data
spp <- names(FIA_spp)

# == Format data ========================================================= ####
# - create: moddat, W_ls, U, Q ------------------------------------------- ####
# modeling data prep (part 1)
moddat_spp <- list(data_spp.std, clim_vars) %>%
  pmap(function(df.x, clim.x) {
    # - grab data ---------------------------------------------------------- ####
    df.x <- df.x$dat %>% as.data.frame()
    
    # - identify variables to manipulate ----------------------------------- ####
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
    
    # - format data: explicitly add new variable forms --------------------- ####
    # add quadratic terms
    for(i in 1:length(var_quad)) {
      varname <- var_quad.name[i]
      var_lin <- var_quad[i]
      
      df.x <- df.x %>% 
        mutate(!!varname := (!!as.name(var_lin))^2)
    }
    
    # add interaction terms
    for(i in 1:length(var_inter)) {
      varname  <- var_inter.name[[i]]
      var_lin1 <- var_inter[[i]][1]
      var_lin2 <- var_inter[[i]][2]
      
      df.x <- df.x %>% 
        mutate(!!varname := (!!as.name(var_lin1)) * (!!as.name(var_lin2)))
    }
    
    # - add ob/un indexes & NA nonfor plots (p_a) -------------------------- ####
    df.x <- df.x %>%
      rename(y_index = plot_obsv) %>%
      # make sure nonfor plots 'NA' for the response
      mutate(p_a = ifelse(y_index == 0, NA, p_a))
    
    # - order variables in df ---------------------------------------------- ####
    # so starting values are associated with the right variables...
    df.x <- df.x %>% select(starts_with(c("gdd5", "hmi", "bio")),
                            everything())
    # - add rownames & make fiahex_id chr ---------------------------------- ####
    df.x$fiahex_id <- as.character(df.x$fiahex_id)
    rownames(df.x) <- df.x$plot_id
    
    # - return df.x -------------------------------------------------------- ####
    return(df.x)
  })

# sparse matrix of neighbor relations (for distinct fiahex_id's)
W_spp <- moddat_spp %>% map(function(df.x) {
  # About: ----
  # Specifying spatial random effects for hex-polygons, not plots; b/c of  
  # intensification on NFS lands, multiple plots can exist within the same
  # FIA hexagon.
  # - Do this in a few steps:
  #   - identify neighboring hex for each fiahex_id (not plot_id)
  #   - identify isolated hex's (i.e. those without any neighbors)
  #   - remove isolated hex's (and their associated plots) from data set
  #   - re-identify hex neighbors with this clean dataset & create sparse W
  #     matrix of neighbor relations from that...
  # - Need to remove isolated hex's so that W (and therefore Sigma of the 
  #   spatial RE) is positive definite.
  # - After creating W_spp, use U_spp to match plot_id's up with their 
  #   respective fiahex_id's
  
  # Part 0: grab relevant data & set n_nb and nn_max_dist ---------------- ####
  # grab coords for fiahex_id's (e.g. not plot_id's) 
  df.x <- df.x %>% 
    select(fiahex_id,
           fiahex_obsv,
           hex_centroid.x_UTM11N_ACTUAL,
           hex_centroid.y_UTM11N_ACTUAL) %>%
    distinct()
  
  # specify max dist for a nearest neighbor
  # - this done through looking at the data and hex calculations (blk nb)
  # - all nn fall within ~5200-5400 m, with next dist ~9km, so chop at 6km
  nn_max_dist <- 6000
  
  # hexagons - so 6 neighbors, plus one for self
  n_nb <- 6 + 1 
  
  # Part 1: identify isolated hex's to drop ------------------------------ ####
  # - grab coord data with rownames -------------------------------------- ####
  xycoord <- data.frame(xcoord = df.x$hex_centroid.x_UTM11N_ACTUAL,
                        ycoord = df.x$hex_centroid.y_UTM11N_ACTUAL)
  rownames(xycoord) <- df.x$fiahex_id  # keep track of rownames
  
  # - get neighbor relations --------------------------------------------- ####
  # hex-grid and add one b/c nabor::knn() counts self
  knnout <- knn(xycoord, k = n_nb)
  
  # only keep the immediately adjacent neighbors sharing a side (and ditch self)
  nb_list <- lapply(1:nrow(df.x), function(i) {
    knnout$nn.idx[i, 2:n_nb][knnout$nn.dist[i, 2:n_nb] <= nn_max_dist]
  })
  
  names(nb_list) <- rownames(xycoord) # keep track of rownames
  
  # - identify isolated hexs & drop from df.x ---------------------------- ####
  # keep track isolated hexagon's fiahex_id's and remove their corresponding
  # plots from df.x 
  hex_isolated <- names(nb_list)[nb_list %>% map(~length(.x) == 0) %>% unlist]
  
  # re-subset list to keep hexagons
  df.x <- df.x[!df.x$fiahex_id %in% hex_isolated, ]
  
  # Part 2: create W sparse matrix --------------------------------------- ####
  # - add z_index -------------------------------------------------------- ####
  # create z_index (to help w. prop-llike calc):
  # - 0: all plots unobserved in hex
  # - 1: all plots observed or mixed in hex
  df.x <- df.x %>% mutate(z_index = ifelse(fiahex_obsv == "unobsv", 0, 1))
  z_index <- df.x$z_index %>% setNames(df.x$fiahex_id)
  
  # - grab coord data with rownames -------------------------------------- ####
  xycoord <- data.frame(xcoord = df.x$hex_centroid.x_UTM11N_ACTUAL,
                        ycoord = df.x$hex_centroid.y_UTM11N_ACTUAL)
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
  
  # - return W and hex_isolated (to drop from moddat) -------------------- ####
  return(list(W = W,
              hex_isolated = hex_isolated,
              z_index = z_index))
})

# modeling data prep (part 2)
np_pre <- moddat_spp %>% map(~nrow(.x)) %>% unlist

# - remove plots occurring in isolated hex (problem for CAR models...)
moddat_spp <- list(moddat_spp, W_spp) %>%
  pmap(function(df.x, hex.x) {
    df.x[!df.x$fiahex_id %in% hex.x$hex_isolated, ]
  })

# - how many hex's was this an issue for?
W_spp %>% map(~length(.x$hex_isolated)) %>% unlist

# - how many plots did we drop?
np_pre - (moddat_spp %>% map(~nrow(.x)) %>% unlist)

# design matrix for the spatial random effects
# - crosswalk between fiahex_id's and the (pot. multiple) plot_id's
U_spp <- list(moddat_spp, W_spp) %>% pmap(function(df.x, W_ls.x) {
  sparse.model.matrix(
    object = plot_id ~ fiahex_id - 1, 
    data = df.x %>%
      mutate(fiahex_id = factor(fiahex_id, levels = names(W_ls.x$z_index))) %>%
      set_rownames(.$plot_id)) # b/c mutate drops them...
  
})

# design matrix for the random effects (ES)
Q_spp <- moddat_spp %>% map(function(df.x) {
  sparse.model.matrix(
    object = plot_id ~ ES - 1, 
    data = df.x %>% set_rownames(.$plot_id)) # b/c mutate drops them...
})

# - create parms & tune_start (shell) ------------------------------------ ####
# About `parms_spp`: 
# - create a list of parms_names, parms_start, and the original names for the
#   parameters (e.g. 'bio1' instead of beta_1, etc...., b/c parms_names are 
#   spp dependent)
# - Note: for all spp, setting a naive start for rho (0.9) & sigma_z (1)
rho_start <- 0.9
sigma_z_start <- 1

parms_spp <- list(parms_start_spp, moddat_spp) %>%
  pmap(function(est.x, dat.x) {
    # - parms_names (chr vector) ----------------------------------------- ####
    covars <- dat.x %>% select(starts_with(c("gdd5", "hmi", "bio"))) %>% names
    betas  <- paste0("beta_", 0:length(covars))
    parms_names <- c(betas, "rho", "sigma_z", "sigma_v", "z_ob", "z_un", "v")
    parms_names_orig <- c("intercept", covars,
                          "rho", "sigma_z", "sigma_v", "z_ob", "z_un", "v")
    
    # - parms_start (names ls) ------------------------------------------- ####
    # get starting values for betas, in the correct order...
    # - sigma_v is always in last slot
    sigma_v_start <- est.x[length(est.x)] %>% unname 
    
    # - simplify est.x down to the betas
    est.x <- est.x[1:(length(est.x) -1)] 
    
    # - reorder if ggd5 exists s.t. ggd5 comes after the fixed intercept
    gdd5_loc   <- grep("gdd5", names(est.x)) # slot for ggd5
    betas_start <- if (length(gdd5_loc) > 0) {
      # re-order: {int, gdd5, all other betas}
      c(est.x[1], est.x[gdd5_loc], est.x[-c(1, gdd5_loc)])
    } else {
      # no re-ordering needed!
      est.x
    }
    
    # create the parms_start list!
    parms_start <- vector("list", length(parms_names)) %>% setNames(parms_names)
    parms_start[betas] <- betas_start
    parms_start[["rho"]] <- rho_start           # defined outside this `pmap`
    parms_start[["sigma_z"]] <- sigma_z_start   # defined outside this `pmap`
    parms_start[["sigma_v"]] <- sigma_v_start
    parms_start[c("z_ob", "z_un", "v")] <- rep(NA, 3)
    
    # - drop "z_un" if all plots were observed --------------------------- ####
    if(!any(dat.x$y_index == 0)) {
      parms_names <- parms_names[!grepl("z_un", parms_names)]
      parms_start <- parms_start[!grepl("z_un", names(parms_start))]
    }
    
    # - return obj as list ----------------------------------------------- ####
    list(parms_names = parms_names,
         parms_start = parms_start,
         parms_names_orig = parms_names_orig)
  })

# About `tune_spp`:
# - create a shell for tuning parameters (will fill when running the sampler)
tune_spp <- parms_spp %>% map(function(parms.x) {
  parms_names <- parms.x$parms_names
  tune_start <- rep(NA, length(parms_names)) %>% setNames(parms_names)
})

# == Create rho_lookup =================================================== ####
# About:
# - creating a look-up table of R and R.det_log for different values of rho
# - why? to speed things up!
# - takes a while to make the look-up table, so: make -> save 
# Requires:
# - W matrix of neighbor relations (data-set dependent!)
# - set rho values

# - get rho values ------------------------------------------------------- ####
# see notes (blk nb): here may want to use (-6:6) or (-8:8) because most data
# sets are closer to 0.99 ... and we don't want to be butting up too close to
# the edge of this space (b/c rho can only take on values between 0 & 1).
fn.inv_logit <- function(x) exp(x) / (1 + exp(x)) # this also in 0-..R script
rho_values <- fn.inv_logit(seq(-8, 8, length.out = 101))

# - create rho_lookup ---------------------------------------------------- ####
rho_lookup_spp <- W_spp %>% map(function(W) {
  W <- W$W

  lapply(rho_values, function(x) {
    R <- Diagonal(x = as.vector(W %*% rep(1, nrow(W)))) - x * W
    R.det_log <- determinant(R, logarithm = TRUE)$modulus # issue w. large W o/w ...

    return(list(rho = x,
                R = R,
                R.det_log = R.det_log))
  })
})

# == Save data & rho_lookup for each spp (rds) =========================== ####
# Why by spp?
# - saving separately for each species to save on space; that is, if we want
#   to work with species one-at-a-time later on, we can w/o the extra stuff.

# save sampler ready data
lapply(spp, function(sp.x) {
  saveRDS(list(sp.x = sp.x,
               moddat = moddat_spp[[sp.x]],
               W_ls = W_spp[[sp.x]],
               U = U_spp[[sp.x]],
               Q = Q_spp[[sp.x]],
               parms_names = parms_spp[[sp.x]]$parms_names,
               parms_start = parms_spp[[sp.x]]$parms_start,
               tune_start = tune_spp[[sp.x]],
               parms_names_orig = parms_spp[[sp.x]]$parms_names_orig),
          file = paste0(path_output,
                        "ch2-5-sampler-ready-data_", sp.x, ".rds"))
})

# save rho_lookup
# - rho_lookup's are big... so save separately from other sampler data...
lapply(spp, function(sp.x) {
  saveRDS(list(sp.x = sp.x,
               rho_lookup = rho_lookup_spp[[sp.x]]),
          file = paste0(path_output,
                        "ch2-5-rholookup-table_", sp.x, ".rds"))
})


