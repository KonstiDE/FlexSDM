# load and install required packages
install.packages("remotes")
library(terra)
library(dplyr)

remotes::install_github("sjevelazco/flexsdm")
library(flexsdm)



########### Data Preparation ###########

# Load raster with environmantal data
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)
names(somevar) <- c("aet", "cwd", "tmx", "tmn")

# Species occurence data (presence-only)
data(hespero)
hespero <- hespero %>% dplyr::select(-id)

# Load ecoregions of California 
regions <- system.file("external/regions.tif", package = "flexsdm")
regions <- terra::rast(regions)
regions <- as.polygons(regions)

# Create subset only containing the ecoregion where the plant can be found
sp_region <- terra::subset(regions, regions$category == "SCR") # ecoregion where *Hesperocyparis stephensonii* is found

# visualize the ecoregion and the species occurrences
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  main = "Hesperocyparis stephensonii occurrences"
)
points(hespero[, c("x", "y")], col = "black", pch = 16)
cols <- rep("gray80", 8)
cols[regions$category == "SCR"] <- "yellow"
terra::inset(
  regions,
  loc = "bottomleft",
  scale = .3,
  col = cols
)



########### Pre-Modelling ##########

# Define the model´s calibrated area
# 25km buffer around  previously loaded points
ca <- calib_area(
  data = hespero,
  x = "x",
  y = "y",
  method = c('buffer', width=25000),
  crs = crs(somevar)
)

# Visualize the species occurrences & calibration area
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  main = "Calibration area and occurrences")
plot(ca, add=TRUE)
points(hespero[,c("x", "y")], col = "black", pch = 16)


# Create pseudo-absence data
# Use the calibration area to sample the equal number of
# pseudo-absences as there are presence points
psa <- sample_pseudoabs(
  data = hespero,
  x = "x",
  y = "y",
  n = sum(hespero$pr_ab), # selecting number of pseudo-absence points that is equal to number of presences
  method = "random",
  rlayer = somevar,
  calibarea = ca
)

# Visualize species presences and pseudo-absences
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  xlim = c(289347, 353284),
  ylim = c(-598052,  -520709),
  main = "Presence (yellow) and Pseudo-absence (black)")
plot(ca, add=TRUE)
points(psa[,c("x", "y")], cex=0.8, pch=16, col = "black") # Pseudo-absences
points(hespero[,c("x", "y")], col = "yellow", pch = 16, cex = 1.5) # Presences

# Bind the dataframes presences and pseudo-absences by rows into a single dataframe
hespero_pa <- bind_rows(hespero, psa)
hespero_pa


# Partition data for evaluating the models (testing and training data)
set.seed(10)
# Partition Method: repeated K-fold cross validation method
hespero_pa2 <- part_random(
  data = hespero_pa,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 5, replicates = 10)
)


# Extract environmental predictors for the presence and pseudo-absence locations
hespero_pa3 <- sdm_extract(
  data = hespero_pa2,
  x = 'x',
  y = 'y',
  env_layer = somevar,
  variables = c('aet', 'cwd', 'tmx', 'tmn')
)



######### Modelling (Standard Models) ##########

# Generalized Linear Model
mglm <- fit_glm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part', # column name with training and validation partition groups
  thr = 'max_sens_spec' # threshold to get binary suitability values (for perfromance metrics)
)

# Generalized Boosted Regression Model
mgbm <- fit_gbm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

# Support Vector Machine Model
msvm <-  fit_svm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

# Spatial prediction from individual models 
mpred <- sdm_predict(
  models = list(mglm, mgbm, msvm),
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)



########## Modelling (Ensemble of Small Models (ESM)) ##########

# When predicting and esm, it is only possible to process one at a time

# Construct Generalized Linear Model using the ESM approach
eglm  <- esm_glm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

# Construct Generalized Boosted Regression Model using the ESM approach
egbm <- esm_gbm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

# Construct Support Vector Machine Model using the ESM approach
esvm <-  esm_svm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

# Predicting ESMs
eglm_pred <- sdm_predict(
  models = eglm ,
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

egbm_pred <- sdm_predict(
  models = egbm ,
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

esvm_pred <- sdm_predict(
  models = esvm,
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)



########## Comparing the Models ###########

# spatial outputs suggest that the standard models tend to predict broader areas
# with high suitability values than the ESM´s
par(mfrow = c(3, 2))
plot(mpred$glm, main = 'Standard GLM')
plot(eglm_pred[[1]], main = 'ESM GLM')

plot(mpred$gbm, main = 'Standard GBM')
plot(egbm_pred[[1]], main = 'ESM GBM')

plot(mpred$svm, main = 'Standard SVM')
plot(esvm_pred[[1]], main = 'ESM SVM')


# Performance metrics

# Merge model performance tables (combined model performance table for all input models)
merge_df <- sdm_summarize(models = list(mglm, mgbm, msvm, eglm, egbm, esvm))

# Create table
knitr::kable(
  merge_df %>% dplyr::select(
    model,
    AUC = AUC_mean,
    TSS = TSS_mean,
    JACCARD = JACCARD_mean,
    BOYCE = BOYCE_mean,
    IMAE = IMAE_mean
  )
)

# AUC, TSS, and Jaccard index are higher for the ESM´s
# Boyce index and the Inverse Mean Absolute Error are slightly higher for the standard models