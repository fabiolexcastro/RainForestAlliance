


## Extract by mask from the table - By each country
## By Fabio Castro-Llanos
## June 27th 2023

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
tble <- read_csv('../data/tbl/climate_index/climate-index_eaf-tea.csv')

# ISO
iso <- 'UGA'

# Filtering
limt <- zone[zone$sov_a3 == iso,]

# =========================================================================
# Extract by mask  --------------------------------------------------------
# =========================================================================
dfrm <- tble %>% 
  rast(., type = 'xyz') %>% 
  crop(., limt) %>%  
  mask(., limt) %>% 
  as.data.frame(., xy = T) %>% 
  as_tibble()
  
nrow(dfrm)
nrow(drop_na(dfrm))

rstr <- rast(dfrm, type = 'xyz')

write.csv(dfrm, file = glue('../data/tbl/climate_index/climate-index_{iso}-tea.csv'), row.names = FALSE)
