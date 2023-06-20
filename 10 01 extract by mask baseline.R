
# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function ----------------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}

# To select the country ---------------------------------------------------
iso <- 'CIV'

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================
zone <- terra::vect('../data/gpkg/westafrica.gpkg')
limt <- zone[zone$GID_0 == iso,]

# Temperature -------------------------------------------------------------
dirs.tmpr <- '//catalogue/WFP_ClimateRiskPr1/1.Data/ERA5' %>% 
  dir_ls() %>% 
  grep('temperature', ., value = T) %>% 
  as.character()

year <- 1990:2021
 
map(.x = dirs.tmpr, .f = function(i){
  cat(green('To process: ', basename(i), '\n'))
  fles <- as.character(dir_ls(i, regexp = '.nc$')) 
  varb <- basename(i) %>% str_split(., '_') %>% map(., 4) %>% unlist() %>% str_sub(., 1, 3) %>% paste0('t', .)
  rstr <- map(.x = year, .f = function(y){
    cat(green('To process: ', y, '\t'))
    fls <- grep(paste0('_', y), fles, value = T)
    rst <- map(.x = 1:12, .f = function(m){
      cat('To process: ', month.abb[m], '\t')
      m <- ifelse(m < 10, paste0('0', m), as.character(m))
      r <- grep(glue('_{y}{m}'), fls, value = T)
      r <- rast(r)
      r <- crop(r, zone)
      r <- mask(r, zone)
      return(r)
    }) %>% 
      reduce(., c)
    dout <- glue('../data/tif/climate_daily/baseline/waf/{varb}')
    terra::writeRaster(x = rst, filename = glue('{dout}/{varb}_{y}.tif'), overwrite = TRUE)
  }) 
  cat('Done!\n')
})

# Precipitation -----------------------------------------------------------
dirs.prec <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps' 
year <- 1990:2021
fles.prec <- dir_ls(dirs.prec, regexp = '.tif$') %>% 
  as.character() %>% 
  grep(paste0(year, collapse = '|'), ., value = T)

map(.x = 1:length(year), .f = function(i){
  cat(green('To process: ', year[i], '\n'))
  fls <- grep(paste0('.', year[i], '.'), fles.prec, value = T)
  rst <- map(.x = 1:12, .f = function(m){
    cat(month.abb[m], '\t')
    m <- mnth(m)
    r <- grep(paste0(year[i], '.', m, '.'), fls, value = T)
    r <- rast(r)
    r <- terra::crop(r, limt)
    r <- terra::mask(r, limt)
    r <- ifel(r < 0, 0, r)
    return(r)
  }) %>% 
    reduce(., c)
  terra::writeRaster(x = rst, filename = glue('../data/tif/climate_daily/baseline/civ-gha-nga-cmr/prec/prec_{year[i]}.tif'), overwrite = TRUE)
})

rm(dirs.prec, dirs.tmpr, dout, fles, fles.prec, fls, varb, yr, r, rs, rst, zone, limt, rstr)

# =========================================================================
# Now to monthly  ---------------------------------------------------------
# =========================================================================

# Temperature 
dirs.tasm <- as.character(dir_ls('../data/tif/climate_daily/baseline/waf', type = 'directory'))
dirs.tasm <- grep('tm', dirs.tasm, value = T)

rstr.tasm <- map(.x = 1:length(dirs.tasm), .f = function(i){
  
  fles <- dirs.tasm[i] %>% dir_ls(., regexp = '.tif$') %>% as.character()
  varb <- dirname(fles) %>% basename %>% unique()
  
  trra <- purrr::map(.x = 1:length(fles), .f = function(y){
    
    rst <- terra::rast(fles[y])
    yea <- parse_number(basename(fles[y]))

    trr <- map(.x = 1:12, .f = function(m){
      
      cat('To process: ', m, '\n')
      m <- mnth(m)
      r <- rst[[grep(paste0(yea, '-', m, '-'), terra::time(rst), value = FALSE)]]
      r <- mean(r)
      names(r) <- glue('{varb}_{yea}-{m}')
      return(r)
      
    }) %>% 
      reduce(., c)
    
    dout <- glue('../data/tif/climate_monthly/baseline/waf/{varb}')
    dir_create(dout)
    terra::writeRaster(x = trr, filename =  glue('{dout}/{varb}_{yea}.tif'), overwrite = TRUE)
    return(trr)
    
  }) %>% 
    reduce(., c)
  
  return(trra)
  
})

# Precipitation
fles.prec <- dir_ls('../data/tif/climate_daily/baseline/prec') %>% 
  as.character() %>% 
  grep('.tif$', ., value = T)

map(.x = 1:length(year), .f = function(i){
  
  cat('To process: ', year[i], '\t')
  rst <- grep(paste0('_', year[i], '.tif'), fles.prec, value = T) %>% 
    terra::rast() 
  
  trr <- map(.x = 1:12, .f = function(m){
    mnt <- mnth(m)
    ps <- grep(paste0(year[i], '.', mnt, '.'), names(rst), value = F)
    rs <- rst[[ps]]
    rs <- app(rs, sum, na.rm = T)
    names(rs) <- glue('prec_{year[i]}-{mnt}')
    return(rs)
  }) %>% 
    reduce(., c)
  
  terra::writeRaster(x = trr, filename = glue('../data/tif/climate_monthly/baseline/prec/prec_{year[i]}.tif'), overwrite = T)
  cat('Done!\n')
    
})


# =========================================================================
# 2022 --------------------------------------------------------------------
# =========================================================================

# Precipitation -----------------------------------------------------------
library(chirps)
chrps <- chirps::get_chirps(object = zone, server = 'CHC', dates = c('2022-01-01', '2022-12-31'))

chrps <- terra::crop(chrps, zone)
chrps <- terra::mask(chrps, zone)
chrps <- terra::ifel(chrps < 0, 0, chrps)
plot(chrps[[1]])

# To write daily data 
terra::writeRaster(x = chrps, filename = '../data/tif/climate_daily/baseline/waf/prec/prec_2022.tif')

# To monthly data
chrps.mnth <- map(.x = 1:12, .f = function(i){
  m <- mnth(i)
  r <- chrps[[grep(paste0('2022.', m, '.'), names(chrps), value = F)]]
  r <- app(r, sum, na.rm = T)
  return(r)
}) %>% 
  reduce(., c)
mnth <- c(paste0('0', 1:9), 10:12)
names(chrps.mnth) <- glue('prec_2022-{mnth}')
terra::writeRaster(x = chrps.mnth, filename = '../data/tif/climate_monthly/baseline/waf/prec/prec_2022.tif')

# Temperature -------------------------------------------------------------

# Daily
dirs <- dir_ls('raw_era5') %>% as.character()
map(.x = dirs, .f = function(i){
  
  fls <- dir_ls(i) %>% as.character()
  var <- basename(fls) %>% str_split('-') %>% map(4) %>% unlist() %>% unique() %>% tolower %>% paste0('t', .)
  
  rst <- map(.x = 1:12, .f = function(m){
    
    cat('To process: ', month.abb[m], '\n')
    m <- mnth(m)
    f <- grep(paste0('2022', m), fls, value = T)
    r <- rast(f)
    r <- crop(r, zone)
    r <- mask(r, zone)
    return(r)
    
  }) %>% 
    reduce(., c)
  
  terra::writeRaster(x = rst, filename = glue('../data/tif/climate_daily/baseline/waf/{var}/{var}_2022.tif'), overwrite = TRUE)
  
})

# Monthly
dirs <- dir_ls('../data/tif/climate_daily/baseline/waf', type = 'directory', regexp = 'tm')
dirs <- as.character(dirs)

map(.x = 1:length(dirs), .f = function(i){
  fls <- dir_ls(dirs[i], regexp ='.tif$') %>% as.character() %>% grep('2022', ., value = T)
  rst <- rast(fls)
  
  trr <- map(.x = 1:12, .f = function(m){
    
    cat('To process: ', m, '\n')
    m <- mnth(m)
    r <- rst[[grep(paste0('2022-', m, '-'), terra::time(rst), value = F)]]
    r <- mean(r)
    return(r)
    
  }) %>% 
    reduce(., c)
  
  
  names(trr) <- glue('tmax_')
  
  
})

rast('../data/tif/climate_daily/baseline/waf/tmax/tmax_2022.tif')
