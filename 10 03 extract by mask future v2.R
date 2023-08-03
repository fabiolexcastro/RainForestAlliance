
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

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================
zone <- terra::vect('../data/gpkg/eastafrica.gpkg')

path <- '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/CMIP6/intermediate/monthly_annual'
fles <- dir_ls(path, regexp = '.tif$')
fles <- as.character(fles)
fles <- grep('ssp585', fles, value = T)

vars <- c('pr', 'tasmax', 'tasmin')
gcms <- basename(fles) %>% str_split(., '_') %>% map_chr(2) %>% unique()

# =========================================================================
# Extract by mask - Daily data --------------------------------------------
# =========================================================================
map(.x = vars, .f = function(v){
  
  fls <- grep(v, fles, value = T)
  
  map(.x = 1:length(fls), .f = function(f){
    
    gcm <- fls[f] %>% 
      basename %>% 
      str_split(., '_') %>% 
      map_chr(2)
    
    rst <- rast(fls[f])
    dte <- as.Date(names(rst), '%Y-%m-%d')
    time(rst) <- dte
    rst <- crop(rst, zone)
    rst <- mask(rst, zone)
    
    
    
  })
  
})








