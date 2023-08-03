


# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================

dirs <- dir_ls('../data/tif/indices/eaf/coffee', type = 'directory') %>% as.character()
fles <- map(.x = dirs, .f = dir_ls, regexp = '.tif$') %>% map(., as.character) %>% unlist()
fles <- grep('bsl-ftr', fles, value = T)
fles <- grep('cls', fles, value = T)

# Read as raster 
rstr <- map(.x = fles, .f = function(i){
  rst <- rast(i)
  nme <- basename(i)
  ind <- str_split(nme, '_') %>% map_chr(1)
  names(rst) <- c(glue('{ind}_bsl'), glue('{ind}_ftr'))
  return(rst)
}) %>% 
  reduce(., c)

terra::writeRaster(x = rstr, filename = '../data/tif/indices/eaf/stck_indices_eaf-coffee.tif', overwrite = TRUE)

# To create the mask 
msks <- map(.x = 1:nlyr(rstr), .f = function(i){
  rst <- rstr[[i]]
  rst <- rst * 0 + 1
  return(rst)
}) %>% 
  reduce(., c)
msks <- app(msks, sum)
poly <- as.polygons(msks)
rstr <- terra::crop(rstr, poly) %>% terra::mask(., poly)
terra::writeRaster(x = rstr, filename = '../data/tif/indices/eaf/stck_indices_eaf-coffee.tif', overwrite = TRUE)



