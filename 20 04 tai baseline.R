
# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function ----------------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.tai <- function(yr){
  
  # yr <- 1990
  cat('To process: ', yr, ' ', mn, '\n')
  tmx <- rast(as.character(grep(yr, fles.tmax, value = T)))
  tmn <- rast(as.character(grep(yr, fles.tmin, value = T)))
  ppt <- rast(as.character(grep(yr, fles.prec, value = T)))
  
  # To change the names
  names(tmx) <- paste0('TMAX_', c(paste0('0', 1:9), 10:12))
  names(tmn) <- paste0('TMIN_', c(paste0('0', 1:9), 10:12))
  names(ppt) <- paste0('PREC_', c(paste0('0', 1:9), 10:12))
  
  # Convert to celcius
  tmx <- tmx - 273.15
  tmn <- tmn - 273.15
  
  # To create average temperature 
  tav <- (tmx + tmn) / 2
  names(tav) <- gsub('TMAX', 'TMEAN', names(tav))
  
  # To calculate temperature range 
  rng <- abs(tmx - tmn)
  
  # Resample solar radiation 
  srd <- terra::resample(srad, rng, method = 'bilinear')
  
  # Assign precipitation names in envirem environment
  envirem::assignNames(solrad = 'SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
  
  # Thornthwaite's Aridity Index
  pet <- envirem::monthlyPET(stack(tav), stack(srd), stack(rng)) 
  pet <- raster::stack(pet)
  names(pet) <- paste0('PET_', c(paste0('0', 1:9), 10:12))
  pet <- raster::resample(pet, stack(ppt), method = 'bilinear')
  
  # Thornwaite aridity index 
  tai <- envirem::aridityIndexThornthwaite(stack(ppt), pet)
  names(tai) <- glue('tai_{yr}')
  return(tai)
  
}

# Load data ---------------------------------------------------------------

# Vector data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- terra::vect('../data/gpkg/westafrica.gpkg')
fles.prec <- dir_ls('../data/tif/climate_monthly/baseline/waf/prec')

# Climate data
fles.tmax <- dir_ls('../data/tif/climate_monthly/baseline/waf/tmax')
fles.tmin <- dir_ls('../data/tif/climate_monthly/baseline/waf/tmin')

# Solar radiation
srad <- dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad', regexp = 'et_solrad_')
srad <- as.character(srad)
srad <- mixedsort(srad)
srad <- rast(srad)
srad <- terra::crop(srad, zone)
srad <- terra::mask(srad, zone)
names(srad) <- glue('srad_{c(paste0("0", 1:9), 10:12)}')

# To create TAI index -----------------------------------------------------
tais <- map(1990:2022, calc.tai)
tais <- map(tais, rast)
tais <- reduce(tais, c)

# Average -----------------------------------------------------------------
tais.avrg <- app(tais, mean, na.rm = T)

# To write ----------------------------------------------------------------
dir.create('../data/tif/indices/waf/tai')
terra::writeRaster(x = tais, filename = '../data/tif/indices/waf/tai/tai_bsl_tsr.tif')
terra::writeRaster(x = tais.avrg, filename = '../data/tif/indices/waf/tai/tai_bsl_avg.tif')

# To reclassify -----------------------------------------------------------
pnts <- suppressMessages(read_csv('../data/tbl/presences/cocoa_westafrica_dupv.csv'))
vles <- terra::extract(tais.avrg, pnts[,1:2])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1), na.rm = T)
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)
clsf <- terra::classify(tais.avrg, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
pnts <- mutate(pnts, value = vles)
write.csv(mtrx, '../data/tbl/reclassify/tai_values.csv', row.names = F)

# =========================================================================
# To make the map ---------------------------------------------------------
# =========================================================================
tble <- terra::as.data.frame(clsf, xy = T) %>% 
  as_tibble() %>% 
  setNames(c('x', 'y', 'value')) %>% 
  mutate(value = factor(value))

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = value)) + 
  scale_fill_manual(values = brewer.pal(n = 3, name = 'YlOrRd'), labels = c('Low', 'Medium', 'High')) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey20') +
  labs(x = 'Lon', y = 'Lat', fill = 'Class') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) +
  ggtitle(label = 'Thornwaite Aridity Index - Westafrica', 
          subtitle = 'Baseline') +
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(face = 'bold', size = 8),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6)) +
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) 

ggsave(plot = gmap, filename = '../png/maps/indices/tai_bsl_rcl_waf.png', units = 'in', width = 9, height = 7, dpi = 300)

# ========================================================================
# Stat ecdf  -------------------------------------------------------------
# ========================================================================
g_stt <- ggplot(pnts) +
  stat_ecdf(geom = "step", pad = FALSE, aes(value)) + 
  labs(x = 'Value', y = 'Percentile') +
  ggtitle(label = 'Thornwaite Aridity Index - Baseline') + 
  theme_minimal() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
        axis.title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', hjust = 0.5, angle = 90),
        legend.position = 'bottom') 

g_stt

ggsave(plot = g_stt, filename = '../png/graphs/indices/tai_bsl_waf.png', units = 'in', width = 8, height = 6, dpi = 300)











