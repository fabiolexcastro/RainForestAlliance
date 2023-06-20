

# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, crayon, rgeos, RColorBrewer, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# =========================================================================
# Function to use ---------------------------------------------------------
# =========================================================================
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.ndd <- function(yr){
  
  # yr <- 1990
  cat('To process: ', yr, '\n')
  rst <- grep(paste0(yr, '.tif'), fles, value = T) %>% terra::rast()
  
  trr <- map(.x = 1:12, .f = function(m){
    
    cat('To process: ', m, '\n')
    m <- mnth(m)
    r <- rst[[grep(paste0(yr, '.', m,  '.'), names(rst), value = FALSE)]]
    f <- terra::app(x = r, fun = function(x){ ndd <- sum(x < 1, na.rm = T); return(ndd)})
    return(f)  
    
  }) %>% 
    reduce(., c)
  
  mnt <- c(glue('0{1:9}'), 10:12)
  names(trr) <- glue('ndd_{yr}-{mnt}')
  return(trr)
  
}

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================
path <- '../data/tif/climate_daily/baseline/waf/prec'
fles <- dir_ls(path) %>% as.character()

# Spatial data
zone <- terra::vect('../data/gpkg/westafrica.gpkg')
wrld <- ne_countries(returnclass = 'sf', scale = 50)

# To apply the function ---------------------------------------------------
ndd.tsr <- map(.x = 1990:2022, .f = calc.ndd)
ndd.tsr <- reduce(ndd.tsr, c)
dou <- glue('../data/tif/indices/waf/ndd')
dir_create(dou)
terra::writeRaster(x = ndd.tsr, filename = glue('{dou}/ndd_bsl_tsr.tif'), overwrite = TRUE)

# =========================================================================
# To calculate the multiaverage monthly -----------------------------------
# =========================================================================
ndd.mnt <- map(.x = 1:12, .f = function(m){
  
  cat('To process: ', month.abb[m], '\n')
  m <- mnth(m)
  r <- ndd.tsr[[grep(paste0('-', m, '$'), names(ndd.tsr), value = F)]]
  r <- app(r, mean, na.rm = T)
  names(r) <- glue('ndd_{m}')
  return(r)
  
}) %>% 
  reduce(., c)

ndd.mnt
terra::writeRaster(x = ndd.mnt, filename = glue('{dou}/ndd_bsl_mnt.tif'), overwrite = TRUE)

# =========================================================================
# To calculate the total average ------------------------------------------
# =========================================================================
ndd.mnt <- terra::rast(glue('{dou}/ndd_bsl_mnt.tif'))
ndd.avg <- terra::app(ndd.mnt, sum, na.rm = T)

terra::writeRaster(x = ndd.avg, filename = glue('{dou}/ndd_bsl_avg.tif'), overwrite = TRUE)

ndd.mnt <- terra::rast(glue('{dou}/ndd_bsl_mnt.tif'))
ndd.avg <- terra::rast(glue('{dou}/ndd_bsl_avg.tif'))

# =========================================================================
# Presences ---------------------------------------------------------------
# =========================================================================
pnts <- suppressMessages(read_csv('../data/tbl/presences/cocoa_westafrica.csv'))[,2:3]
dupv <- duplicated(terra::extract(ndd.avg, pnts[,1:2], cell = T)[,3])
pnts <- pnts[!dupv,]
write.csv(pnts, '../data/tbl/presences/cocoa_westafrica_dupv.csv', row.names = F)

# Extract the value for the presences -------------------------------------
vles <- terra::extract(ndd.avg, pnts[,1:2])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1))
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)
clsf <- terra::classify(ndd.avg, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
pnts <- mutate(pnts, value = vles)

dir.create('../data/tbl/reclassify')
write.csv(mtrx, '../data/tbl/reclassify/ndd_values.csv', row.names = F)

# =========================================================================
# To make the map ---------------------------------------------------------
# =========================================================================
tble <- terra::as.data.frame(clsf, xy = T) %>% 
  as_tibble() %>% 
  setNames(c('x', 'y', 'value')) %>% 
  mutate(value = factor(value))

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = value)) + 
  scale_fill_manual(values = brewer.pal(n = 3, name = 'BrBG'), labels = c('Low', 'Medium', 'High')) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey20') +
  labs(x = 'Lon', y = 'Lat', fill = 'NDD Class') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) +
  ggtitle(label = 'Number of dry days - Westafrica', 
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

gmap

ggsave(plot = gmap, filename = '../png/maps/indices/ndd_bsl_rcl_waf.png', units = 'in', width = 9, height = 7, dpi = 300)

# ========================================================================
# Stat ecdf  -------------------------------------------------------------
# ========================================================================
g_stt <- ggplot(pnts) +
  stat_ecdf(geom = "step", pad = FALSE, aes(value)) + 
  labs(x = 'Value', y = 'Percentile') +
  ggtitle(label = 'Number of dry days for Westafrica - Baseline') + 
  theme_minimal() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
        axis.title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', hjust = 0.5, angle = 90),
        legend.position = 'bottom') 

g_stt

ggsave(plot = g_stt, filename = '../png/graphs/indices/ndd_bsl_waf.png', units = 'in', width = 8, height = 6, dpi = 300)
