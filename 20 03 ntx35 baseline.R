
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
zone <- terra::vect('../data/gpkg/westafrica.gpkg')
fles <- dir_ls('../data/tif/climate_daily/baseline/waf/tmax', regexp = '.tif$')
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
dir.create('../data/tif/indices/waf/ntx35', recursive = TRUE)
terra::writeRaster(x = ntxr, filename = '../data/tif/indices/waf/ntx35/ntx35_bsl_tsr.tif')

# Monthly 
ntxr.mnth <- map(.x = 1:12, .f = function(m){
  
  cat('To process: ', month.abb[m], '\n')
  m <- mnth(m)
  r <- ntxr[[grep(glue('-{m}$'), names(ntxr), value = F)]]
  a <- terra::app(r, mean, na.rm = T)
  return(a)
  
}) %>% 
  reduce(., c)

terra::writeRaster(x = ntxr.mnth, filename = '../data/tif/indices/waf/ntx35/ntx35_bsl_mnt.tif')

# Total average
ntxr.summ <- app(ntxr.mnth, sum, na.rm = T)
terra::writeRaster(x = ntxr.summ, filename = '../data/tif/indices/waf/ntx35/ntx35_bsl_sum.tif')

# To reclassify -----------------------------------------------------------
pnts <- suppressMessages(read_csv('../data/tbl/presences/cocoa_westafrica_dupv.csv'))
vles <- terra::extract(ntxr.summ, pnts[,1:2])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1))
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)
clsf <- terra::classify(ntxr.summ, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
pnts <- mutate(pnts, value = vles)

write.csv(mtrx, '../data/tbl/reclassify/ntx35_values.csv', row.names = F)

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
  ggtitle(label = 'Heat stress for crops (days), Tmax > 35ºC - Westafrica', 
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

ggsave(plot = gmap, filename = '../png/maps/indices/ntx35_bsl_rcl_waf.png', units = 'in', width = 9, height = 7, dpi = 300)

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











