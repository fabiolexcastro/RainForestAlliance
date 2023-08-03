
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
dirs <- dir_ls('../data/tif/climate_monthly/baseline/eaf', type = 'directory')
dirs <- as.character(dirs)
vars <- c('prec', 'tmax', 'tmea', 'tmin')
dirs <- grep(paste0(vars, collapse = '|'), dirs, value = T)

wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- wrld[wrld$sov_a3 %in% c('KEN', 'UGA'),]

# =========================================================================
# Forecast function -------------------------------------------------------
# =========================================================================
make.forecast <- function(dir){
  
  dir <- dirs[1]
  
  var <- basename(dir)
  fls <- dir_ls(dir, regexp = '.tif$') %>% as.character()
  
  # Read as raster file
  rst <- map(.x = fls, .f = rast)
  rst <- reduce(rst, c)
  
  # Raster to table
  tbl <- terra::as.data.frame(rst, xy = T)
  tbl <- as_tibble(tbl)
  tbl <- mutate(tbl, gid = 1:nrow(tbl))
  tbl <- gather(tbl, var, value, -gid, -x, -y)
  tbl <- separate(tbl, col = 'var', into = c('variable', 'date'), sep = '_')
  gds <- pull(tbl, gid) %>% unique()
  
  tbl %>% filter(date == '1990-01' & variable == 'prec') %>% dplyr::select(x, y, value) %>% rast(., type = 'xyz') %>% plot()
  
  options(future.globals.maxSize = 12000 * 1024^2)
  plan(cluster, workers = 20, gc = TRUE)
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
  
  crd <- distinct(tbl, gid, x, y)
  rsl <- bind_rows(rsl)
  rsl <- inner_join(rsl, crd, by = 'gid') %>% dplyr::select(x, y, gid, date, value) %>% mutate(variable = basename(dir), .after = gid) 
  write.csv(rsl, file = glue('../data/tbl/forecast/eaf/{basename(dir)}_frc_raw.csv'), row.names = F)
  cat('Done\n')
  
}

# Check the results 
fles <- as.character(dir_ls('../data/tbl/forecast/eaf', regexp = '.csv'))

# Summarise function
smmr.frcs <- function(x){
  d <- grep(x, fles, value = T) %>% 
    read_csv() %>% 
    mutate(month = str_sub(date, 6, nchar(date)), 
           month = as.numeric(month)) %>% 
  # d <- d %>% 
    group_by(gid, x, y, month) %>% 
    dplyr::summarise(value = mean(value, na.rm = T)) %>% 
    ungroup() %>% 
    dplyr::select(-gid) %>% 
    spread(month, value)  
  r <- rast(d, type = 'xyz')
  return(r)
}

prec <- smmr.frcs(x = 'prec')
tmax <- smmr.frcs(x = 'tmax')
tmin <- smmr.frcs(x = 'tmin')
tavg <- (tmax + tmin) / 2

tmax <- tmax - 273.15
tmin <- tmin - 273.15
tavg <- tavg - 273.15

plot(vect(zone))
plot(tmax[[1]], add = T)

terra::writeRaster(x = tmax, filename = '../data/tif/climate_monthly/forecast/eaf/tmax_frcs.tif', overwrite = T)
terra::writeRaster(x = tmin, filename = '../data/tif/climate_monthly/forecast/eaf/tmin_frcs.tif', overwrite = T)
terra::writeRaster(x = tavg, filename = '../data/tif/climate_monthly/forecast/eaf/tavg_frcs.tif', overwrite = T)
terra::writeRaster(x = prec, filename = '../data/tif/climate_monthly/forecast/eaf/prec_frcs.tif', overwrite = T)

tmin
tavg
tmax
prec

prec.bsln <- terra::rast('../data/tif/climate_monthly/baseline/eaf/prec/prec_1990.tif')
