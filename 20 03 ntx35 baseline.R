
# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(RColorBrewer, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function ----------------------------------------------------------------
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.ntx <- function(yr, mn, thr = 35){
  
  # yr <- 2017; mn <- 10
  cat('To process: ', yr, ' ', mn, '\n')
  rst <- rast(grep(yr, fles, value = T))
  rst <- rst[[grep(glue('-{mnt}-'), terra::time(rst), value = F)]]
  rst <- rst - 273.15
  rsl <- terra::app(x = rst, fun = function(x){ntxval = sum(x >= thr, na.rm = T); return(ntxval)})
  rsl <- terra::crop(rsl, zone)
  rsl <- terra::mask(rsl, zone)
  return(rsl)
  
}

# Load data ---------------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- terra::vect('../data/gpkg/eastafrica.gpkg')
fles <- dir_ls('../data/tif/climate_daily/baseline/eaf/tmax', regexp = '.tif$')
fles <- as.character(fles)
dtes <- expand.grid(month = c(glue('0{1:9}'), 10:12), year = 1990:2022)

# To apply the function  --------------------------------------------------

# Time series
ntxr <- purrr::map(.x = 1:nrow(dtes), .f = function(i){
  try(expr = {
    cat("To process row: ", i, '\n')
    dts <- dtes[i,]
    mnt <- pull(dts, 1)
    mnt <- as.character(mnt)
    yea <- pull(dts, 2)
    ntx <- calc.ntx(yr = yea, mn = mnt, thr = 35)
    names(ntx) <- glue('ntx_{yea}-{mnt}')
    return(ntx)
  })
})
ntxr <- reduce(ntxr, c)
dir.create('../data/tif/indices/eaf/ntx35', recursive = TRUE)
terra::writeRaster(x = ntxr, filename = '../data/tif/indices/eaf/ntx35/ntx35_bsl_tsr.tif')

ntxr <- terra::rast('../data/tif/indices/eaf/ntx35/tea/ntx35_bsl_tsr.tif')

# Yearly average 
year <- 1990:2022
ntxr.yrs <- map(.x = 1:length(year), .f = function(y){
  
  cat('To process: ', year[y], '\n')
  r <- ntxr[[grep(glue('_{year[y]}-'), names(ntxr), value = F)]]
  a <- terra::app(r, sum, na.rm = T)
  names(a) <- glue('ntx_{year[y]}')
  return(a)
  
}) %>% 
  reduce(., c)

terra::writeRaster(x = ntxr.yrs, filename = '../data/tif/indices/eaf/ntx35/ntx35_bsl_tsr-avg.tif')

# Total average
ntxr.avg <- app(ntxr.yrs, mean, na.rm = T)
names(ntxr.avg) <- glue('ntxr_avg')
terra::writeRaster(x = ntxr.avg, filename = '../data/tif/indices/eaf/ntx35/ntx35_bsl_avg.tif')

# To read the future and make a stack with both periods -------------------
ntxr.bsl <- ntxr.avg
ntxr.bsl <- rast('../data/tif/indices/eaf/ntx35/ntx35_bsl_avg.tif')
ntxr.ftr <- terra::rast('../data/tif/indices/eaf/ntx35/ntx_ftr_avg.tif')

ntxr.bsl <- terra::resample(ntxr.bsl, ntxr.ftr, method = 'bilinear')
stck <- c(ntxr.bsl, ntxr.ftr)

# To reclassify -----------------------------------------------------------
pnts <- suppressMessages(read_csv('../data/tbl/presences/cocoa_westafrica_dupv.csv'))
pnts <- suppressMessages(read_csv('../data/tbl/presences/tea_east_africa.csv')) %>% filter(iso %in% c('KEN', 'UGA'))
vles <- terra::extract(ntxr.bsl, pnts[,1:2])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1), na.rm = T)
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)

# Coffee at east africa
pnts <- read_csv('//alliancedfs.alliance.cgiar.org/CL9_Coffee_Cocoa2/_coffeeEastAfrica/tbl/occ/occ_swd')
pnts <- mutate(pnts, iso = terra::extract(zone, pnts[,c('Lon', 'Lat')])$sov_a3, .before = Lon)
pnts <- filter(pnts, iso %in% c('UGA', 'KEN'))
pnts <- dplyr::select(pnts, Lon, Lat)

vles <- terra::extract(ntxr.bsl, pnts[,1:2])
qntl <- quantile(vles, c(0, 0.5, 0.75, 1), na.rm = T)
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)

# To make the reclassify
clsf <- terra::classify(stck, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
names(clsf) <- c('Baseline', 'Future')
pnts <- mutate(pnts, value = vles)

write.csv(mtrx, '../data/tbl/reclassify/ntx35_values-coffee.csv', row.names = F)

dir.create('../data/tif/indices/eaf/ntx35/coffee')
terra::writeRaster(x = clsf, filename = '../data/tif/indices/eaf/coffee/ntx35/ntx35_bsl-ftr_cls.tif')

# =========================================================================
# To make the map ---------------------------------------------------------
# =========================================================================
tble <- terra::as.data.frame(clsf, xy = T) %>% 
  as_tibble() %>% 
  gather(var, value, -x, -y) %>% 
  mutate(value = factor(value), var = factor(var, levels = c('Baseline', 'Future')))

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = factor(value))) + 
  scale_fill_manual(values = brewer.pal(n = 3, name = 'YlOrRd'), labels = c('Low', 'High'), na.translate = F) +
  facet_wrap(~var) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey20') +
  labs(x = 'Lon', y = 'Lat', fill = 'Class') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) +
  ggtitle(label = 'Heat stress for crops (days), Tmax > 35ºC - Eastafrica (Tea)', 
          subtitle = '') +
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

ggsave(plot = gmap, filename = '../png/maps/indices/ntx35_bsl-ftr_rcl_eaf-coffee.png', units = 'in', width = 9, height = 7, dpi = 300)

# ========================================================================
# Stat ecdf  -------------------------------------------------------------
# ========================================================================
g_stt <- ggplot(pnts) +
  stat_ecdf(geom = "step", pad = FALSE, aes(value)) + 
  labs(x = 'Value', y = 'Percentile') +
  ggtitle(label = 'Human Heat Stress - Baseline') + 
  theme_minimal() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
        axis.title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', hjust = 0.5, angle = 90),
        legend.position = 'bottom') 

g_stt

ggsave(plot = g_stt, filename = '../png/graphs/indices/ntx35_bsl_waf.png', units = 'in', width = 8, height = 6, dpi = 300)











