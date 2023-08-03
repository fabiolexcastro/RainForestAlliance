

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
  trr <- sum(trr)
  names(trr) <- glue('ndd_{yr}')
  return(trr)
  
}

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================
path <- '../data/tif/climate_daily/baseline/eaf/prec'
fles <- dir_ls(path) %>% as.character()

# Spatial data
wrld <- rnaturalearth::ne_countries(returnclas = 'sf', scale = 50)
zone <- terra::vect('../data/gpkg/westafrica.gpkg')
zone <- terra::vect('../data/gpkg/eastafrica.gpkg')
zone <- wrld[wrld$sov_a3 %in% c('KEN', 'UGA'),]

# To apply the function ---------------------------------------------------
plan(cluster, workers = 12, gc = TRUE)
# ndd.tsr <- map(.x = 1990:2022, .f = calc.ndd)
ndd.tsr <- future_map(.x = 1990:2022, .f = calc.ndd)
future:::ClusterRegistry('stop')
ndd.tsr <- reduce(ndd.tsr, c)
dou <- glue('../data/tif/indices/eaf/ndd')
dir_create(dou)
terra::writeRaster(x = ndd.tsr, filename = glue('{dou}/ndd_bsl_tsr.tif'), overwrite = TRUE)

# =========================================================================
# To calculate the total average ------------------------------------------
# =========================================================================

ndd.avg <- terra::app(ndd.tsr, mean, na.rm = T)
ndd.avg <- terra::crop(ndd.avg, zone) %>% terra::mask(., zone)
terra::writeRaster(x = ndd.avg, filename = glue('../data/tif/indices/eaf/ndd/ndd_bsl_avg.tif'), overwrite = TRUE)

ndd.bsl <- terra::rast(glue('../data/tif/indices/eaf/ndd/ndd_bsl_avg.tif'))
ndd.ftr <- terra::rast(glue('../data/tif/indices/eaf/ndd/ndd_ftr_avg.tif'))

# =========================================================================
# Presences ---------------------------------------------------------------
# =========================================================================

# Cocoa
pnts <- suppressMessages(read_csv('../data/tbl/presences/cocoa_westafrica.csv'))[,2:3]
pnts <- suppressMessages(read_csv('../data/tbl/presences/tea_east_africa.csv'))
pnts <- filter(pnts, iso %in% c('UGA', 'KEN'))
pnts <- dplyr::select(pnts, x, y)
# dupv <- duplicated(terra::extract(ndd.bsl, pnts[,1:2], cell = T)[,3])
# pnts <- pnts[!dupv,]

# Coffee
pnts <- read_csv('//alliancedfs.alliance.cgiar.org/CL9_Coffee_Cocoa2/_coffeeEastAfrica/tbl/occ/occ_swd')
pnts <- mutate(pnts, iso = terra::extract(vect(zone), pnts[,c('Lon', 'Lat')])$sov_a3, .before = Lon)
pnts <- filter(pnts, iso %in% c('UGA', 'KEN'))
pnts <- dplyr::select(pnts, Lon, Lat)

# Extract the value for the presences -------------------------------------
vles <- terra::extract(ndd.bsl, pnts[,1:2])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1))
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)

dir.create('../data/tbl/reclassify/eaf/coffee')
write.csv(mtrx, '../data/tbl/reclassify/eaf/coffee/coffee_ndd_values.csv', row.names = F)

# To reclassify -----------------------------------------------------------
names(ndd.bsl) <- 'ndd_bsl'
names(ndd.ftr) <- 'ndd_ftr'
ndd.ftr <- terra::resample(ndd.ftr, ndd.bsl)

stk <- c(ndd.bsl, ndd.ftr)

rcl <- terra::classify(stk, mtrx, include.lowest = T)
dir.create('../data/tif/indices/eaf/coffee/ndd')
terra::writeRaster(rcl, filename = glue('../data/tif/indices/eaf/coffee/ndd/ndd_bsl-ftr_cls.tif'))

# =========================================================================
# To make the map ---------------------------------------------------------
# =========================================================================
tble <- terra::as.data.frame(rcl, xy = T) %>% 
  as_tibble() %>% 
  setNames(c('x', 'y', 'Baseline', 'Future')) %>% 
  gather(var, value, -x, -y) %>%
  mutate(var = factor(var, levels = c('Baseline', 'Future')),
         value = factor(value))

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = value)) + 
  scale_fill_manual(values = rev(brewer.pal(n = 3, name = 'BrBG')), labels = c('Low', 'Medium', 'High')) +
  facet_wrap(.~ var) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey20') +
  labs(x = 'Lon', y = 'Lat', fill = 'NDD Class') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) +
  ggtitle(label = 'Number of dry days - Eastafrica', 
          subtitle = '') +
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        strip.text = element_text(face = 'bold', hjust = 0.5),
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

ggsave(plot = gmap, filename = '../png/maps/indices/ndd_bsl-ftr_rcl_eaf-coffee.png', units = 'in', width = 12, height = 5, dpi = 300)

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
