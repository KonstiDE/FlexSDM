install.packages("remotes")
library(terra)
library(dplyr)

remotes::install_github("sjevelazco/flexsdm")
library(flexsdm)


# Load points and locate them in a visualization plot
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)
names(somevar) <- c("aet", "cwd", "tmx", "tmn")

# species occurence data (presence-only)
data(hespero)
hespero <- hespero %>% dplyr::select(-id)

# California ecoregions
regions <- system.file("external/regions.tif", package = "flexsdm")
regions <- terra::rast(regions)
regions <- as.polygons(regions)
sp_region <- terra::subset(regions, regions$category == "SCR") # ecoregion where *Hesperocyparis stephensonii* is found

# visualize the species occurrences
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


# Buffer the previously loaded points with pre-modeling calib_area() method
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



# Sample the same number of species presences
psa <- sample_pseudoabs(
  data = hespero,
  x = "x",
  y = "y",
  n = sum(hespero$pr_ab) * 3, # selecting number of pseudo-absence points that is equal to number of presences
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


hespero_pa <- bind_rows(hespero, psa)
hespero_pa



set.seed(10)
# Repeated K-fold method
hespero_pa2 <- part_random(
  data = hespero_pa,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 5, replicates = 10)
)

hespero_pa3 <- sdm_extract(
  data = hespero_pa2,
  x = 'x',
  y = 'y',
  env_layer = somevar,
  variables = c('aet', 'cwd', 'tmx', 'tmn')
)

mglm <- fit_glm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)
mgbm <- fit_gbm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)
msvm <-  fit_svm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)

mpred <- sdm_predict(
  models = list(mglm, mgbm, msvm),
  pred = somevar,
  con_thr = TRUE,
  predict_area = ca
)

eglm  <- esm_glm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)
egbm <- esm_gbm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)
esvm <-  esm_svm(
  data = hespero_pa3,
  response = 'pr_ab',
  predictors = c('aet', 'cwd', 'tmx', 'tmn'),
  partition = '.part',
  thr = 'max_sens_spec'
)


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

par(mfrow = c(3, 2))
plot(mpred$glm, main = 'Standard GLM')
plot(eglm_pred[[1]], main = 'ESM GLM')

plot(mpred$gbm, main = 'Standard GBM')
plot(egbm_pred[[1]], main = 'ESM GBM')

plot(mpred$svm, main = 'Standard SVM')
plot(esvm_pred[[1]], main = 'ESM SVM')



merge_df <- sdm_summarize(models = list(mglm, mgbm, msvm, eglm, egbm, esvm))

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



