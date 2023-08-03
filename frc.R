
## Forecast analysis - project to 2023:2025
## By: Fabio Castro-Llanos

# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(terra, forecast, fs, sf, tidyverse, crayon, furrr, future, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================

root <- '../data/tif'
dirs <- dir_ls('../data/tif/climate_monthly/baseline/waf', type = 'directory')
dirs <- as.character(dirs)
vars <- c('prec', 'tmax', 'tmea', 'tmin')
dirs <- grep(paste0(vars, collapse = '|'), dirs, value = T)

zone <- terra::vect('../data/gpkg/westafrica.gpkg')

# =========================================================================
# Forecast function -------------------------------------------------------
# =========================================================================
make.forecast <- function(dir, iso){
  
  # dir <- dirs[3]
  # iso <- 'GHA'
  
  var <- basename(dir)
  fls <- dir_ls(dir, regexp = '.tif$') %>% as.character()
  
  # Read as raster file
  rst <- map(.x = fls, .f = rast)
  rst <- reduce(rst, c)
  zne <- zone[zone$GID_0 == iso,]
  
  # Extract by mask 
  rst <- terra::crop(rst, zne)
  rst <- terra::mask(rst, zne)
  
  # Raster to table
  rst <- rst * 1
  tbl <- terra::as.data.frame(rst, xy = T)
  tbl <- as_tibble(tbl)
  tbl <- mutate(tbl, gid = 1:nrow(tbl))
  tbl <- gather(tbl, var, value, -gid, -x, -y)
  tbl <- separate(tbl, col = 'var', into = c('variable', 'date'), sep = '_')
  gds <- pull(tbl, gid) %>% unique()
  crd <- distinct(tbl, gid, x, y)
  
  options(future.globals.maxSize = 8000 * 1024^2)
  plan(cluster, workers = 18, gc = TRUE)
  rsl <- furrr::future_map(.x = gds, .f = function(i){
    cat(i, '\n')
    tb <- filter(tbl, gid == i)
    vl <- pull(tb, value)
    ts <- ts(vl, start = c(1990, 1), end = c(2022, 12), frequency = 12)
    ft <- auto.arima(ts)
    fc <- forecast(ft, h = 36)
    rs <- as.vector(fc$upper[,2])
    rs <- tibble(gid = i, date = c(glue('2023-{1:12}'), glue('2024-{1:12}'), glue('2025-{1:12}')), value = rs)
    return(rs)
  })
  future:::ClusterRegistry('stop')
  rsl <- bind_rows(rsl)
  rsl <- inner_join(rsl, crd, by = 'gid') %>% dplyr::select(x, y, gid, date, value) %>% mutate(variable = basename(dir), .after = gid) 
  write.csv(rsl, file = glue('../data/tbl/forecast/{iso}/{basename(dir)}_frc_raw.csv'), row.names = F)
  
}

# To apply the forecast ---------------------------------------------------
make.forecast(dir = dirs[1], iso = 'GHA')
make.forecast(dir = dirs[2], iso = 'GHA')
make.forecast(dir = dirs[3], iso = 'GHA')
make.forecast(dir = dirs[4], iso = 'GHA')

# To make the summarise ---------------------------------------------------
fles <- dir_ls('../data/tbl/forecast/GHA', regexp = '.csv$') %>% as.character()
tble <- map(fles, read_csv)
tble <- bind_rows(tble)
tble <- mutate(tble, month = as.numeric(str_sub(date, 6, nchar(date))))

# Prec / Tmax / Tmin
prec <- tble %>% filter(variable == 'prec') %>% group_by(x, y, gid, month) %>% summarise(value = mean(value, na.rm = T)) %>% ungroup() %>% arrange(gid) %>% dplyr::select(x, y, month, value) %>% spread(month, value) %>% rast(., type = 'xyz')
tmax <- tble %>% filter(variable == 'tmax') %>% group_by(x, y, gid, month) %>% summarise(value = mean(value, na.rm = T)) %>% ungroup() %>% arrange(gid) %>% dplyr::select(x, y, month, value) %>% spread(month, value) %>% rast(., type = 'xyz')
tmin <- tble %>% filter(variable == 'tmin') %>% group_by(x, y, gid, month) %>% summarise(value = mean(value, na.rm = T)) %>% ungroup() %>% arrange(gid) %>% dplyr::select(x, y, month, value) %>% spread(month, value) %>% rast(., type = 'xyz')

names(prec) <- glue('prec_{1:12}')
names(tmax) <- glue('tmax_{1:12}')
names(tmin) <- glue('tmin_{1:12}')

terra::writeRaster(x = prec, filename = glue('../data/tif/climate_monthly/forecast/GHA/prec_mnt.tif'), overwrite = TRUE)
terra::writeRaster(x = tmin, filename = glue('../data/tif/climate_monthly/forecast/GHA/tmin_mnt.tif'), overwrite = TRUE)
terra::writeRaster(x = tmax, filename = glue('../data/tif/climate_monthly/forecast/GHA/tmax_mnt.tif'), overwrite = TRUE)
