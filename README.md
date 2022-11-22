### These R scripts are in support of 'Spatial-Bayesian models project shifts in suitable habitat for Pacific Northwest tree species under climate change.'
#### Karin Kralicek<sup>1,2</sup>, Jay M. Ver Hoef<sup>3</sup>, Tara M. Barrett<sup>4</sup>, and Temesgen Hailemariam<sup>2</sup> 
<sup>1</sup>Forest Inventory and Analysis Program, Rocky Mountain Research Station, USDA Forest Service, Fort Collins, CO, United States; <sup>2</sup>Forest Measurements and Biometrics Laboratory, Department of Forest Engineering, Resources, and Management, Oregon State University, Corvallis, OR, United States; <sup>3</sup>Marine Mammal Laboratory, Alaska Fisheries Science Center, NOAA National Marine Fisheries Service, Seattle, WA, United States; <sup>4</sup>Pacific Northwest Research Station, USDA Forest Service, Wenatchee, WA, United States

The manuscript has been accepted by Ecosphere.

---
### There are 13 scripts in this project, which begin with a number indicating their order in the workflow. 

#### Definition of scripts:

* 0-functions.R: contains function definitions and is sourced in other project scripts in order to keep those other scripts tidy. 

* 1-load-manip-raw-data.R: loads and formats data to be used in subsequent scripts, including FIA plot data & associated PRISM-Norm81m data, species codes, association of FIA plots with their containing FIA hexagons and membership to Ecological Sections (ES's).

* 2-add-bioclimatic-vars.R: calculates bio-style climate variables (from PRISM-Norm81m data) for each FIA plot.

* 3-explore-data.R: preliminary exploration of the data to identify biologically meaningful climate variables to consider for modeling.

* 4-model-forms.R: identifies the final species-specific model forms for which spatial-Bayesian hierarchical models will be fit.

* 5-ready-data-rholookup-for-sampler.R: formats the data to be sampler-ready for script 6; also creates the rho look-up tables for each species.

* 6-fit-CAR-models.R: tune the CAR SDM sampler for each species; includes all adjustments from the very first run through the final burn-in and 2000 MCMC sample from each species' model. Code to check models while tuning is included towards the end of this script.

* 7-predict-maps.R: generates predictions maps for the current and four future climate scenarios. For each climate data set, 2000 maps are produced which correspond to the 2000 MCMC samples kept after burn-in for each species-specific model.

* 8-voronoi-polygons.R: creates spatial data to use in subsequent scripts for plotting and for area-based calculations, for both plot-level and hexagon-level results. Note on backward use: script 6 reference objects created in this script (8) to examine the organization of the spatial random effects while tuning the models.

* 9-analyze-results.R: Calculates summary statistics based on prediction maps and base data sets. This includes table summaries for plot hexagon or ES membership, conceptual diagrams in Figures 1 & 2, posterior distribution summaries for parameters, species-response curve results, model checks of effective sample size and AUC, etc. 

* 10-maps-and-range-size-results.R: creates maps and plots for mean-predictions and prediction uncertainty for the current and four future climate scenarios, as well as maps associated with the VP (Voronoi polygon) method of range estimation.

* 11-kde-range-size-results.R: obtain range-size estimates based on the BKDE (bivariate kernel density estimation) method of range estimation. Calculates other statistics and compares range estimates between the BKDE and VP methods.

* 12-model-check-morans-i.R: performs a type of Bayesian model-check by evaluating model performance with Moran's I at different distance classes.


#### Note on data: 

The raw data supporting the conclusions of this manuscript are based on confidential, precise-location coordinates for FIA field data. Imprecise-location (fuzzed) versions of raw FIA field data are available at https://apps.fs.usda.gov/fia/datamart/datamart.html. For more information on access to precise-location FIA data see Burrill et al. (2021). Requests to access these datasets should be directed to https://www.fia.fs.usda.gov/about/about_us/#questions-requests.

#### Citations:

Burrill, E. A., DiTommaso, A. M., Turner, J. A., Pugh, S. A., Menlove, J., Christiansen, G., et al. (2021). The Forest Inventory and Analysis Database: database description and user guide version 9.0.1 for Phase 2. U.S. Department of Agriculture, Forest Service. Available online at: http://www.fia.fs.usda.gov/library/database-documentation/
