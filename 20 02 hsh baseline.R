
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
calc.hsh.bsl <- function(yr){
  
  # yr <- 1990
  
  cat("To process: ", yr, '\n')
  
  # Files 
  fl.hm <- grep(paste0('_', yr), fles.hmdt, value = T)
  fl.tx <- grep(paste0('_', yr), fles.tmax, value = T)
  fl.tn <- grep(paste0('_', yr), fles.tmin, value = T)
  
  # As raster fiels and extract by mask 
  rs.hm <- rast(fl.hm) %>% crop(., zone) %>% mask(., zone)
  rs.tx <- rast(fl.tx) %>% crop(., zone) %>% mask(., zone)
  rs.tn <- rast(fl.tn) %>% crop(., zone) %>% mask(., zone)
  rs.hm <- terra::resample(rs.hm, rs.tx, method = 'bilinear')
  
  rs.tx <- rs.tx - 273.15
  rs.tn <- rs.tn - 273.15
  
  # Constant
  c1 = -8.78469475556
  c2 =  1.61139411
  c3 =  2.33854883889
  c4 = -0.14611605
  c5 = -0.012308094
  c6 = -0.0164248277778
  c7 =  2.211732 * 10^(-3)
  c8 =  7.2546 * 10^(-4)
  c9 = -3.582 * 10^(-6)
  
  # To calculate the average temperature 
  rs.ta <- (rs.tx + rs.tn) / 2
  
  hs <- purrr::map(.x = 1:12, .f = function(i){
    heat_idx <- function(ta, rh){
      hi <- ifelse(ta >= 25, c1 + (c2*ta) + (c3*rh) + (c4*ta*rh) + (c5*ta^2) + (c6*rh^2) + (c7*ta^2*rh) + (c8*ta*rh^2) + (c9*ta^2*rh^2), ta)
      return(hi)
    }
    hm <- rs.hm[[i]]
    ta <- rs.ta[[i]]
    HI <- terra::lapp(terra::sds(ta, hm), fun = heat_idx)
    HI_avg <- mean(HI, na.rm = T) %>% terra::mask(., zone)
    HI_max <- max(HI, na.rm = T) %>% terra::mask(., zone)
    return(list(HI_avg, HI_max))
  })
  
  hs.av <- do.call('c', map(hs, 1))
  hs.mx <- do.call('c', map(hs, 2))
  names(hs.av) <- glue('hs.avg_{yr}_{1:12}')
  names(hs.mx) <- glue('hs.max_{yr}_{1:12}')
  
  # plot(hs.av)
  dout <- glue('../data/tif/indices/waf/hsh')
  # terra::writeRaster(x = hs.av, filename = glue('{dout}/hsh_{yr}.tif'), overwrite = T)
  cat('Done everything!\n')
  return(hs.av)
  
}

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================
path <- '../data/tif/climate_daily/baseline/eaf'
dirs <- dir_ls(path, type = 'directory') %>% grep('tm', ., value = T) %>% as.character()

# Humidity
path.hmdt <- '//catalogue/WFP_ClimateRiskPr1/1.Data/ERA5/2m_relative_humidity'
fles.hmdt <- dir_ls(path.hmdt, regexp = '.nc$') %>% as.character() # fles.hmdt <- dir_ls('./raw_era5/2m_relative_humidity', regexp = '.nc$') %>% as.character()

# Tmax 
path.tmax <- '../data/tif/climate_monthly/baseline/eaf/tmax'
fles.tmax <- dir_ls(path.tmax, regexp = '.tif$') %>% as.character()

# Tmin
path.tmin <- '../data/tif/climate_monthly/baseline/eaf/tmin'
fles.tmin <- dir_ls(path.tmin, regexp = '.tif$') %>% as.character()

# Spatial data
zone <- terra::vect('../data/gpkg/eastafrica.gpkg')
wrld <- ne_countries(returnclass = 'sf', scale = 50)

# =========================================================================
# Humidity to monthly -----------------------------------------------------
# =========================================================================
map(.x = 1990:2021, .f = function(i){
  
  cat('To process: ', i, '\n')
  fls <- grep(paste0('_', i), fles.hmdt, value = T)
  rst <- map(.x = 1:12, .f = function(m){
    cat('Process: ', m, '\t')
    m <- mnth(m)
    r <- grep(paste0('_', i, m), fles.hmdt, value = T)
    r <- rast(r)
    r <- terra::crop(r, zone)
    r <- terra::mask(r, zone)    
    r <- app(r, mean, na.rm = T)
    return(r)
  }) %>% 
    reduce(., c)
  
  dir <- glue('../data/tif/climate_monthly/baseline/eaf/hmdt')
  mns <- c(glue('0{1:9}'), 10:12)
  names(rst) <- glue('hmdt_{i}-{mns}')
  terra::writeRaster(x = rst, filename = glue('{dir}/hmdt_{i}.tif'), overwrite = T)
  cat('Done!\n')  
  
})

# -------------------------------------------------------------------------
# Continue ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Humidity 
path.hmdt <- '../data/tif/climate_monthly/baseline/eaf/hmdt'
fles.hmdt <- dir_ls(path.hmdt, regexp = '.tif$') %>% as.character()

# To apply the function  --------------------------------------------------

# Time series
hsh.tsr <- map(.x = 1990:2021, .f = calc.hsh.bsl)
hsh.tsr <- reduce(hsh.tsr, c)
terra::writeRaster(x = hsh.tsr, filename = '../data/tif/indices/eaf/hsh/hsh_bsl_tsr.tif', overwrite = TRUE)

# To calculate the average 
yrs <- 1990:2021

hsh.avg <- map(.x = 1:length(yrs), .f = function(i){
  cat('To process: ', yrs[i],  '\n')
  r <- hsh.tsr[[grep(yrs[i], names(hsh.tsr), value = F)]]
  a <- app(r, mean, na.rm = T)
  names(a) <- glue('hsh_{yrs[i]}')
  return(a)
}) %>% 
  reduce(., c)
terra::writeRaster(x = hsh.avg, filename = '../data/tif/indices/eaf/hsh/hsh_bsl_tsr-avg.tif', overwrite = T)

# Total average 
hsh.avg <- app(hsh.avg, mean, na.rm = T)
terra::writeRaster(x = hsh.avg, filename = '../data/tif/indices/eaf/hsh/hsh_bsl_avg.tif', overwrite = TRUE)

# Read future data --------------------------------------------------------
hsh.bsl <- rast('../data/tif/indices/eaf/hsh/hsh_bsl_avg.tif')
hsh.ftr <- rast('../data/tif/indices/eaf/hsh/hsh_ftr_avg.tif')

hsh.bsl <- terra::resample(hsh.bsl, hsh.ftr, method = 'bilinear')
names(hsh.bsl) <- 'hsh_bsl'
stk <- c(hsh.bsl, hsh.ftr)

# To reclassify -----------------------------------------------------------
mtrx <- matrix(c(-Inf, 25, 1, 25, 30, 2, 30, 35, 3, 35, Inf, 4), byrow = T, ncol = 3)

# To classify
clsf <- terra::classify(stk, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
names(clsf) <- c('Baseline', 'Future')

terra::writeRaster(clsf, filename = '../data/tif/indices/eaf/hsh/hsh_bsl-ftr_cls.tif')

# Mild or no stress: < 25
# Caution: 25–30
# Danger: 30-35
# Extreme danger: > 35

write.csv(mtrx, '../data/tbl/reclassify/hsh_values.csv', row.names = F)

# =========================================================================
# To make the map ---------------------------------------------------------
# =========================================================================
tble <- terra::as.data.frame(clsf, xy = T) %>% 
  as_tibble() %>% 
  gather(var, value, -x, -y) %>% 
  mutate(var = factor(var, levels = c('Baseline', 'Future')))

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = factor(value))) + 
  scale_fill_manual(values = brewer.pal(n = 4, name = 'YlOrRd'), labels = c('Mild or no stress', 'Caution', 'Danger', 'Extreme danger')) +
  facet_wrap(~var) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey20') +
  labs(x = 'Lon', y = 'Lat', fill = 'Class') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) +
  ggtitle(label = 'Human Heat Stress - Westafrica', 
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

ggsave(plot = gmap, filename = '../png/maps/indices/hsh_bsl-ftr_rcl_waf.png', units = 'in', width = 12, height = 5, dpi = 300)

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

ggsave(plot = g_stt, filename = '../png/graphs/indices/hsh_bsl_waf.png', units = 'in', width = 8, height = 6, dpi = 300)









