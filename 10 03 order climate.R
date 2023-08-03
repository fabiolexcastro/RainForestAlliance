

## Order the climate data 
## By Fabio Castro-Llanos
## June 26th 2023

# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(terra, fs, sf, tiydverse, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

zone <- terra::vect('../data/gpkg/eastafrica.gpkg')
 
# =========================================================================
# Functions ---------------------------------------------------------------
# =========================================================================
read.bsln.mnth <- function(dir = '../data/tif/climate_monthly/baseline/eaf/prec', var = 'prec'){
  rst <- dir %>% 
    dir_ls(., regexp = var) %>% 
    as.character() %>% 
    rast()
  rst <- map(.x = 1:12, .f = function(m){
    cat('To process: ', month.abb[m], '\n')
    m <- ifelse(m < 10, paste0('0', m), as.character(m))
    r <- rst[[grep(paste0('-', m), names(rst), value = F)]]
    r <- app(r, mean)
    names(r) <- glue('{var}_{m}_bsl')
    return(r)
  }) %>% 
    reduce(., c)
  return(rst)
}

rclf.binr <- function(stk, thr){
  rsl <- map(.x = 1:12, .f = function(i){
    rcl <- stk[[i]]
    rcl <- terra::classify(rcl, matrix(c(-Inf, -1*thr, -1, -1*thr, sd.tmax, 0, thr, Inf, 1), byrow = T, ncol = 3))
    return(rcl)
  }) %>% 
    do.call('c', .)
  return(rsl)
}
chgnm <- function(stk, var, prd){names(stk) <- glue('{var}_{1:12}_{prd}'); return(stk)}

# =========================================================================
# Baseline ----------------------------------------------------------------
# =========================================================================

# To read
prec.bsln <- read.bsln.mnth(dir = '../data/tif/climate_monthly/baseline/eaf/prec', var = 'prec')
tmax.bsln <- read.bsln.mnth(dir = '../data/tif/climate_monthly/baseline/eaf/tmax', var = 'tmax')
tmin.bsln <- read.bsln.mnth(dir = '../data/tif/climate_monthly/baseline/eaf/tmin', var = 'tmin')

mask <- prec.bsln[[1]] * 0 + 1

# Kelvin to celcius
tmax.bsln <- tmax.bsln - 273.15
tmin.bsln <- tmin.bsln - 273.15

# To resample 
tmax.bsln <- terra::resample(tmax.bsln, mask, method = 'bilinear')
tmin.bsln <- terra::resample(tmin.bsln, mask, method = 'bilinear')
tavg.bsln <- (tmax.bsln + tmin.bsln) / 2

tavg.bsln <- chgnm(stk = tavg.bsln, var = 'tavg', 'bsl')

# Solar radiation 
srad <- '//catalogue/workspace-cluster9/DATA/ET_SolRad' %>% 
  dir_ls(., regexp = 'et_solrad') %>% 
  gtools::mixedsort() %>% 
  rast() %>% 
  crop(., zone) %>% 
  mask(., zone) %>% 
  resample(., mask, method = 'bilinear')

# To calculate ETP 
etps.bsln <- 0.0013 * 0.408 * srad * (tavg.bsln + 17) * (tmax.bsln - tmin.bsln - 0.0123 * prec.bsln) ^ 0.76
etps.bsln <- etps.bsln * c(31,29,31,30,31,30,31,31,30,31,30,31)
for(i in 1:12){etps.bsln[[i]][which.lyr(is.na(etps.bsln[[i]]))] <- 0}
etps.bsln <- terra::crop(etps.bsln, zone) %>% terra::mask(., zone)
names(etps.bsln) <- glue('etps_{c(paste0(0, 1:9), 10:12)}_bsl')

# To calculate the balance
baln.bsln <- prec.bsln - etps.bsln
mtrx.baln <- matrix(c(-Inf, -10, -1, -10, 10, 0, 0, Inf, 1), ncol = 3, nrow = 3, byrow = T)

# Precipitation: 100
# Potential evapotranspiration: 20
# Difference 100 - 20 = 80 (reclassify value: 1 (wet))

baln.bsln.rclf <- classify(baln.bsln, mtrx.baln)
names(baln.bsln.rclf) <- glue('baln_{c(paste0(0, 1:9), 10:12)}_bsl')

# To write the results 
terra::writeRaster(x = prec.bsln, filename = '../data/tif/climate_monthly/baseline/eaf/prec_bsln.tif', overwrite = T)
terra::writeRaster(x = tmin.bsln, filename = '../data/tif/climate_monthly/baseline/eaf/tmin_bsln.tif', overwrite = T)
terra::writeRaster(x = tavg.bsln, filename = '../data/tif/climate_monthly/baseline/eaf/tavg_bsln.tif', overwrite = T)
terra::writeRaster(x = tmax.bsln, filename = '../data/tif/climate_monthly/baseline/eaf/tmax_bsln.tif', overwrite = T)
terra::writeRaster(x = etps.bsln, filename = '../data/tif/climate_monthly/baseline/eaf/etps_bsln.tif', overwrite = T)
terra::writeRaster(x = baln.bsln.rclf, filename = '../data/tif/climate_monthly/baseline/eaf/baln-rclf_bsln.tif', overwrite = TRUE)

# =========================================================================
# Forecast ----------------------------------------------------------------
# =========================================================================
prec.frcs <- rast('../data/tif/climate_monthly/forecast/eaf/prec_frcs.tif')
tmin.frcs <- rast('../data/tif/climate_monthly/forecast/eaf/tmin_frcs.tif')
tmax.frcs <- rast('../data/tif/climate_monthly/forecast/eaf/tmax_frcs.tif')
tavg.frcs <- (tmax.frcs + tmin.frcs) / 2

# To change the names for these rastes
prec.frcs <- chgnm(stk = prec.frcs, var = 'prec', 'frc')
tmin.frcs <- chgnm(stk = tmin.frcs, var = 'tmin', 'frc')
tmax.frcs <- chgnm(stk = tmax.frcs, var = 'tmax', 'frc')
tavg.frcs <- chgnm(stk = tavg.frcs, var = 'tavg', 'frc')

# Convert Kelvin to Celcius
tmin.frcs <- tmin.frcs - 273.15
tmax.frcs <- tmax.frcs - 273.15
tavg.frcs <- tavg.frcs - 273.15

# To resample 
tmin.frcs <- terra::resample(tmin.frcs, mask, method = 'bilinear')
tmax.frcs <- terra::resample(tmax.frcs, mask, method = 'bilinear')
tavg.frcs <- terra::resample(tavg.frcs, mask, method = 'bilinear')

# To calculate ETP 
etps.frcs <- 0.0013 * 0.408 * srad * (tavg.frcs + 17) * (tmax.frcs - tmin.frcs - 0.0123 * prec.frcs) ^ 0.76
etps.frcs <- etps.frcs * c(31,29,31,30,31,30,31,31,30,31,30,31)
for(i in 1:12){etps.frcs[[i]][which.lyr(is.na(etps.frcs[[i]]))] <- 0}
etps.frcs <- terra::crop(etps.frcs, zone) %>% terra::mask(., zone)
names(etps.frcs) <- glue('etps_{c(paste0(0, 1:9), 10:12)}_frc')

# To calculate the balance
baln.frcs <- prec.frcs - etps.frcs
mtrx.baln <- matrix(c(-Inf, -10, -1, -10, 10, 0, 0, Inf, 1), ncol = 3, nrow = 3, byrow = T)

# Precipitation: 100
# Potential evapotranspiration: 20
# Difference 100 - 20 = 80 (reclassify value: 1 (wet))

baln.frcs.rclf <- classify(baln.frcs, mtrx.baln)
names(baln.frcs.rclf) <- glue('baln_{c(paste0(0, 1:9), 10:12)}_frc')

# To write the results 
terra::writeRaster(x = prec.frcs, filename = '../data/tif/climate_monthly/forecast/eaf/prec_frcs.tif', overwrite = T)
terra::writeRaster(x = tmin.frcs, filename = '../data/tif/climate_monthly/forecast/eaf/tmin_frcs.tif', overwrite = T)
terra::writeRaster(x = tavg.frcs, filename = '../data/tif/climate_monthly/forecast/eaf/tavg_frcs.tif', overwrite = T)
terra::writeRaster(x = tmax.frcs, filename = '../data/tif/climate_monthly/forecast/eaf/tmax_frcs.tif', overwrite = T)
terra::writeRaster(x = etps.frcs, filename = '../data/tif/climate_monthly/forecast/eaf/etps_frcs.tif', overwrite = T)
terra::writeRaster(x = baln.frcs.rclf, filename = '../data/tif/climate_monthly/forecast/eaf/baln-rclf_frcs.tif', overwrite = TRUE)

# =========================================================================
# Future ------------------------------------------------------------------
# =========================================================================
prec.ftre <- rast('../data/tif/climate_monthly/future/eaf/prec_avg.tif')
tmax.ftre <- rast('../data/tif/climate_monthly/future/eaf/tmax_avg.tif')
tmin.ftre <- rast('../data/tif/climate_monthly/future/eaf/tmin_avg.tif')
tavg.ftre <- (tmax.ftre + tmin.ftre) / 2

names(prec.ftre) <- paste0(names(prec.ftre), '-ftr')
names(tmin.ftre) <- glue('tmin_{c(paste0(0, 1:9), 10:12)}_ftr')
names(tavg.ftre) <- glue('tavg_{c(paste0(0, 1:9), 10:12)}_ftr')
names(tmax.ftre) <- glue('tmax_{c(paste0(0, 1:9), 10:12)}_ftr')

# To calculate Potential Evapotranspiration
prec.ftre.avrg <- prec.ftre[[grep('_avg', names(prec.ftre), value = FALSE)]]
etps.ftre <- 0.0013 * 0.408 * srad * (tavg.ftre + 17) * (tmax.ftre - tmin.ftre - 0.0123 * prec.ftre.avrg) ^ 0.76
etps.ftre <- etps.ftre * c(31,29,31,30,31,30,31,31,30,31,30,31)
for(i in 1:12){etps.ftre[[i]][which.lyr(is.na(etps.ftre[[i]]))] <- 0}
etps.ftre <- terra::crop(etps.ftre, zone) %>% terra::mask(., zone)
names(etps.ftre) <- glue('etps_{c(paste0(0, 1:9), 10:12)}_ftr')

# To calculate the balance
baln.ftre <- prec.ftre.avrg - etps.ftre
mtrx.baln <- matrix(c(-Inf, -10, -1, -10, 10, 0, 0, Inf, 1), ncol = 3, nrow = 3, byrow = T)

baln.ftre.rclf <- classify(baln.ftre, mtrx.baln)
names(baln.ftre.rclf) <- glue('baln_{c(paste0(0, 1:9), 10:12)}_ftr')

# To write the results 
terra::writeRaster(x = tavg.ftre, filename = '../data/tif/climate_monthly/future/eaf/tavg_avg.tif', overwrite = T)
terra::writeRaster(x = etps.ftre, filename = '../data/tif/climate_monthly/future/eaf/etps_ftr.tif', overwrite = T)
terra::writeRaster(x = baln.ftre.rclf, filename = '../data/tif/climate_monthly/future/eaf/baln-rclf_ftr.tif', overwrite = TRUE)

# To make the stack 
prec.bsln
tmax.ftre
tavg.ftre
tmin.ftre
baln.bsln

prec.ftre
tmax.ftre
tmin.ftre
tavg.ftre
baln.ftre


# =========================================================================
# Temperature difference --------------------------------------------------
# =========================================================================
tmax.dfrn <- chgnm(tmax.ftre - tmax.bsln, 'tmax', 'ftr.dfr')
tmin.dfrn <- chgnm(tmin.ftre - tmin.bsln, 'tmin', 'ftr.dfr')
tavg.dfrn <- chgnm(tavg.ftre - tavg.bsln, 'tmin', 'ftr.dfr')

sd.tmax <- as.numeric(global(app(tmax.dfrn, 'sd'), 'min', na.rm = T))
sd.tmin <- as.numeric(global(app(tmin.dfrn, 'sd'), 'min', na.rm = T))
sd.tavg <- as.numeric(global(app(tavg.dfrn, 'sd'), 'min', na.rm = T))

tmax.binr <- rclf.binr(stk = tmax.dfrn, thr = sd.tmax)
tavg.binr <- rclf.binr(stk = tavg.dfrn, thr = sd.tavg)
tmin.binr <- rclf.binr(stk = tmin.dfrn, thr = sd.tmin)

tavg.binr
tmax.binr <- chgnm(tmax.binr, 'tmax', 'ftr.dfr')
tavg.binr <- chgnm(tavg.binr, 'tavg', 'ftr.dfr')
tmin.binr <- chgnm(tmin.binr, 'tmin', 'ftr.dfr')

# =========================================================================
# To make the stack -------------------------------------------------------
# =========================================================================

stck <- c(prec.bsln, tmin.bsln, tavg.bsln, tmax.bsln, etps.bsln, baln.bsln.rclf,
          prec.frcs, tmin.frcs, tavg.frcs, tmax.frcs, etps.frcs, baln.frcs.rclf,
          prec.ftre, tmin.ftre, tavg.ftre, tmax.ftre, etps.ftre, baln.ftre.rclf, 
          tmin.binr, tavg.binr, tmax.binr)

stck <- terra::crop(stck, zone)
stck <- terra::mask(stck, zone)

# Raster to table 
tble <- terra::as.data.frame(stck, xy = T) %>% as_tibble()
drop <- drop_na(tble)
rstr <- rast(drop, type = 'xyz')

# Now to write the raster of climate dataset ------------------------------
terra::writeRaster(x = rstr, filename = '../data/tif/climate_monthly/stack_climate_eaf.tif', overwrite = TRUE)
rstr <- rast('../data/tif/climate_monthly/stack_climate_eaf.tif')

# =========================================================================
# To add the index rasters ------------------------------------------------
# =========================================================================

indx <- rast('../data/tif/indices/eaf/stck_indices_eaf-tea_v2.tif')
indx <- rast('../data/tif/indices/eaf/stck_indices_eaf-coffee.tif')

# Extract by mask - Raster
mask <- indx[[1]] * 0 + 1
mask <- as.polygons(mask)

rstr <- crop(rstr, mask) %>% mask(., mask)
stck <- terra::resample(stck, rstr, method = 'near')

rstr <- c(stck, indx)
tble <- terra::as.data.frame(rstr, xy = T)
tble <- as_tibble(tble)
tble <- drop_na(tble)
rstr <- terra::rast(tble, type = 'xyz')

# Raster to table again 
tble <- as_tibble(terra::as.data.frame(rstr, xy = TRUE))
colnames(tble)
write.csv(tble, '../data/tbl/climate_index/climate-index_eaf-tea.csv', row.names = F)

tble <- read_csv('../data/tbl/climate_index/climate-index_CIV.csv')
colnames(tble)
