# Install flexsdm
remotes::install_github("sjevelazco/flexsdm")
library(flexsdm)

# Load the other required libs
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(geodata)
library(raster)
library(purrr)

# Load regional gulf of mexcio data for a nicer picture
usa <- read_sf("data/countries/USA_shp/gadm36_USA_0.shp")
mex <- read_sf("data/countries/MEX_shp/gadm36_MEX_0.shp")

# Load the predictor variables file
varfile <- rast("data/clim/varfile_gulf.tif")

# Load the shapefile, where the fish theoretically could appear
shape <- read_sf("data/occs/data_0.shp")

# Load the occurences of the fish
points <- read.csv("data/occs/points_data.csv")

# Crop the USA and MEX shapefile for better plotting
usa <- st_crop(usa, shape)
mex <- st_crop(mex, shape)



# !!!Coding task!!!
# ----------------------------------------------------------------------------------
# | Convert the points csv to an sf object an plot everything together in ggplot() |
# ----------------------------------------------------------------------------------
points_sf <- st_as_sf(points, coords = 19:20, crs = crs(shape))

ggplot() +
  geom_sf(data = usa) +
  geom_sf(data = mex) +
  geom_sf(data = shape, aes(fill=1), alpha = 0.5) +
  geom_sf(data = points_sf, col = "red")

# -----------------------------------------------------------------------------------


# Sampling Pseudo-absense
psa <- sample_pseudoabs(
  data = points,
  x = "longitude",
  y = "latitude",
  n = nrow(points),
  method = "random",
  rlayer = varfile,
  calibarea = shape
)

psa_sf <- psa %>%
  as_tibble() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = crs(shape))

ggplot() +
  geom_sf(data = usa) +
  geom_sf(data = mex) +
  geom_sf(data = shape, aes(fill=1), alpha = 0.5) +
  geom_sf(data = points_sf, col = "green") +
  geom_sf(data = psa_sf, col = "yellow")


# !!!Coding task!!!
# ------------------------------------------------------------------------------------------------------------
# | Select the presnece and geometry variable from presenses and pseudo absences and combine them in one df. |
# | Hint: The outcome should be a dataframe with two columns and 14 rows.                                    |
# ------------------------------------------------------------------------------------------------------------

points_sf_selected <- points_sf[,c("presence", "geometry")]
colnames(psa_sf) <- c("presence", "geometry")

points_complete <- rbind(points_sf_selected, psa_sf)

# -----------------------------------------------------------------------------------------------------------
# Optional: Split the location into two seperate x and y columns:

separated_coord <- points_complete %>%
    mutate(x = unlist(map(points_complete$geometry,1)),
           y = unlist(map(points_complete$geometry,2)))

# -----------------------------------------------------------------------------------------------------------


# Random Partitioning of both presence and pseudo-absense
random_partitioning <- part_random(
  data = separated_coord,
  pr_ab = "presence",
  method = c(method = "rep_kfold", folds = 5, replicates = 10)
)

# Remove NA values from the predictors, replace them with zero
varfile <- ifel(is.nan(varfile), 0, varfile)
plot(varfile[[1]])

sdm_extract <- sdm_extract(
  data = random_partitioning,
  x = 'x',
  y = 'y',
  env_layer = varfile,
  variables = c("sst", "chlor_a", "nflh", "poc"),
  filter_na = FALSE
)

# Fitting model
msvm <-  fit_svm(
  data = sdm_extract,
  response = 'presence',
  predictors = c("sst", "chlor_a", "nflh", "poc"),
  partition = '.part',
  thr = 'max_sens_spec'
)
mpred <- sdm_predict(
  models = list(msvm),
  pred = varfile,
  con_thr = TRUE,
  predict_area = vect(shape)
)

# Essamble model
esvm <-  esm_svm(
  data = sdm_extract,
  response = 'presence',
  predictors = c("sst", "chlor_a", "nflh", "poc"),
  partition = '.part',
  thr = 'max_sens_spec'
)
esvm_pred <- sdm_predict(
  models = list(esvm),
  pred = varfile,
  con_thr = TRUE,
  predict_area = vect(shape)
)

# SpatRaster to sf for nicer plotting
coords_svm <- xyFromCell(mpred$svm, cell = 1:ncell(mpred$svm))
raster_df_svm <- data.frame(x = coords_svm[, 1], y = coords_svm[, 2], value = as.vector(mpred$svm))

coords_esvm <- xyFromCell(esvm_pred[[1]], cell = 1:ncell(esvm_pred[[1]]))
raster_df_esvm <- data.frame(x = coords_esvm[, 1], y = coords_esvm[, 2], value = as.vector(esvm_pred[[1]]))


ggplot() +
  geom_sf(data = usa) +
  geom_sf(data = mex) +
  geom_sf(data = shape, alpha = 0.5) +
  geom_raster(data = raster_df_svm, mapping = aes(x, y, fill=value)) +
  scale_fill_gradientn(
    colours = viridis::viridis(10),
    na.value = "transparent",
    guide_colorbar(title = "Probability")
  ) +
  xlab("") +
  ylab("")

merge_df <- sdm_summarize(models = list(esvm))

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





