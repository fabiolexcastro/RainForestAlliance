

## Extract by mask for the future data
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
dirs <- dir_ls('../data/tif/climate_monthly/future/waf', type = 'directory')
dirs <- as.character(dirs)
vars <- c('prec', 'tmax', 'tmea', 'tmin')
dirs <- grep(paste0(vars, collapse = '|'), dirs, value = T)