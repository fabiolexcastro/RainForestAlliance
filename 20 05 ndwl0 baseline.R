
# =========================================================================
# Load libraries ----------------------------------------------------------
# =========================================================================
require(pacman)
pacman::p_load(RColorBrewer, envirem, terra, fs, sf, tidyverse, crayon, rgeos, cmocean, ggspatial, gtools, raster, glue, lubridate, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Function ----------------------------------------------------------------
peest <- function(srad, tmin, tmean, tmax){
  
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  
  # Soil heat flux parameters
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Net radiation
  rn <- (1-albedo) * srad
  
  # Soil heat flux
  eslope <- a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
  
  # Estimate vpd
  esat_min <- 0.61120*exp((17.67*tmin)/(tmin+243.5))
  esat_max <- 0.61120*exp((17.67*tmax)/(tmax+243.5))
  vpd <- vpd_cte*(esat_max-esat_min) #kPa
  
  # Priestley-Taylor
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  pt_coef <- pt_fact*pt_const
  pt_coef <- 1 + (pt_coef-1) * vpd / vpd_ref
  
  #*10^6? To convert fluxes MJ to J
  #rlat_ht? Latent heat flux to water flux
  #100/rho_w? Kg/m^2 to cm
  et_max <- (pt_coef * rn * eslope/(eslope+psycho) * 10^6 / rlat_ht * 100/rho_w)*10 #in mm
  return(et_max)
}
eabyep_calc <<- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
  
  avail   <- min(avail, soilcp)
  
  # ERATIO
  percwt <- min(avail/soilcp*100, 100)
  percwt <- max(percwt, 1)
  eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
  
  demand  <- eratio * evap
  result  <- avail + rain - demand
  logging <- result - soilcp
  logging <- max(logging, 0)
  logging <- min(logging, soilsat)
  # runoff  <- result - logging + soilcp
  avail   <- min(soilcp, result)
  avail   <- max(avail, 0)
  # runoff  <- max(runoff, 0)
  
  out     <- list(Availability = c(AVAIL, avail),
                  # Demand       = demand,
                  # Eratio       = eratio,
                  Logging      = logging
  )
  return(out)
}
mnth <- function(x){ifelse(x < 10, paste0('0', x), as.character(x))}
calc.ndwl0 <- function(yr, mn){
  
  # yr <- 1990
  # mn <- 1
  
  cat('To process: ', yr, ' ', mn, '\n')
  mnt <- mnth(mn)
  
  # Solar radiation
  srd <- rast(as.character(grep(glue('_{yr}{mnt}'), fles.srad, value = T)))
  srd <- terra::crop(srd, zone)
  srd <- terra::mask(srd, zone)
  srd <- srd/1000000
  
  # Temperature
  tmx <- rast(grep(glue('_{yr}'), fles.tmax, value = T))
  tmx <- tmx[[grep(glue('{yr}-{mnt}'), time(tmx), value = F)]]
  
  tmn <- rast(grep(glue('_{yr}'), fles.tmin, value = T))
  tmn <- tmn[[grep(glue('{yr}-{mnt}'), time(tmn), value = F)]]
  
  tmx <- tmx - 273.15
  tmn <- tmn - 273.15
  
  # Precipitation
  ppt <- rast(grep(glue('_{yr}'), fles.prec, value = T))
  ppt <- ppt[[grep(glue('{yr}.{mnt}'), names(ppt), value = F)]]
  ref <- ppt[[1]] * 0
  names(ref) <- 'ref'
  
  # To make the resample 
  tmx <- terra::resample(tmx, ref, method = 'bilinear')
  tmn <- terra::resample(tmn, ref, method = 'bilinear')
  srd <- terra::resample(srd, ref, method = 'bilinear')
  
  # To calculate average temperature
  tav <- (tmx + tmn) / 2
  
  # To calculate ETP
  ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
  
  # Compute water balance model
  date <- paste0(yr,'-',mnt)
  if(date %in% c('1990-01','2026')){
    cat('Start\n')
    AVAIL <<- ref
    AVAIL[!is.na(AVAIL)] <- 0
  } else {
    cat('Continous\n')
    AVAIL <<- terra::rast(glue('../data/tif/soils/avail.tif'))
  }
  
  # Fill gaps SCP - SST
  scp <- terra::focal(scp, w = 9, fun = mean, na.policy = 'only', na.rm=TRUE) 
  sst <- terra::focal(sst, w = 9, fun = mean, na.policy = 'only', na.rm=TRUE) 
  scp <- terra::crop(scp, zone) %>% terra::mask(., zone)
  sst <- terra::crop(sst, zone) %>% terra::mask(., zone)
  
  scp <- terra::resample(scp, ref)
  sst <- terra::resample(sst, ref)
  ppt <- terra::resample(ppt, ref)
  
  watbal <- map(.x = 1:nlyr(ETMAX), .f = function(i){
    wtr <- eabyep_calc(soilcp  = scp,
                       soilsat = sst,
                       avail   = AVAIL[[terra::nlyr(AVAIL)]],
                       rain    = ppt[[i]],
                       evap    = ETMAX[[i]])
    AVAIL <<- wtr$Availability
    return(wtr)
  })
  
  LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
  NDWL0  <- sum(LOGGING > 0)  
  terra::writeRaster(x = AVAIL[[nlyr(AVAIL)]], glue('../data/tif/soils/avail.tif'), overwrite = TRUE)
  rm(ETMAX, LOGGING, ppt, tmx, tmn, tav, watbal)
  return(NDWL0)
  
}

# =========================================================================
# Load data ---------------------------------------------------------------
# =========================================================================

# Vector data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
zone <- terra::vect('../data/gpkg/westafrica.gpkg')

# Climate data
fles.tmax <- dir_ls('../data/tif/climate_daily/baseline/waf/tmax', regexp = '.tif$')
fles.tmin <- dir_ls('../data/tif/climate_daily/baseline/waf/tmin', regexp = '.tif$')
fles.prec <- dir_ls('../data/tif/climate_daily/baseline/waf/prec', regexp = '.tif$')

# Solar radiation 
fles.srad <- dir_ls('//catalogue/WFP_ClimateRiskPr1/1.Data/ERA5/solar_radiation_flux', regexp = '.nc$')

# Soils -------------------------------------------------------------------
scp <- rast('../data/tif/soils/africa_scp.tif') %>% terra::crop(., zone) %>% terra::mask(., zone)
sst <- rast('../data/tif/soils/africa_ssat.tif')%>% terra::crop(., zone) %>% terra::mask(., zone)

# =========================================================================
# To apply the function ---------------------------------------------------
# =========================================================================
years <- 1990:2021
month <- 1:12

ndwl <- map(.x = 1:length(years), .f = function(i){
  yea <- years[i]
  fnl <- map(.x = 1:length(month), .f = function(j){
    ndwl <- calc.ndwl0(yr = yea, mn = month[j])  
    return(ndwl)
  }) %>% 
    reduce(., c) %>% 
    app(., sum, na.rm = TRUE)
  names(fnl) <- glue('NDWL0_{yea}')
  return(fnl)
}) %>% 
  reduce(., c)

# To write time series
terra::writeRaster(x = ndwl, filename = '../data/tif/indices/waf/ndwl0/ndwl0_bsl_tsr.tif')

ndwl <- terra::rast('../data/tif/indices/eaf/ndwl0/ndwl0_bsl_tsr.tif')

# To calculate the average ------------------------------------------------
ndwl.avg <- app(ndwl, mean, na.rm = T)
terra::writeRaster(x = ndwl.avg, filename = '../data/tif/indices/waf/ndwl0/ndwl0_bsl_avg.tif')

# To read the future results ----------------------------------------------
ndwl.bsln <- ndwl.avg
ndwl.bsln <- rast('../data/tif/indices/eaf/ndwl0/ndwl0_bsl_avg.tif')
ndwl.bsln <- app(ndwl.bsln, sum)
ndwl.ftre <- rast('../data/tif/indices/eaf/ndwl0/ndwl0_ftr_avg.tif')

ndwl.bsln
ndwl.ftre

stck <- c(ndwl.bsln, ndwl.ftre)
names(stck) <- c('Baseline', 'Future')

# =========================================================================
# Presences ---------------------------------------------------------------
# =========================================================================
pnts <- suppressMessages(read_csv('../data/tbl/presences/cocoa_westafrica_dupv.csv'))
pnts <- suppressMessages(read_csv('../data/tbl/presences/tea_east_africa.csv')) %>% filter(iso %in% c('UGA', 'KEN'))

# Extract the value for the presences -------------------------------------
vles <- terra::extract(ndwl.bsln, pnts[,1:2])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1), na.rm = T)
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 367, 3), byrow = T, ncol = 3)

# To classify
clsf <- terra::classify(stck, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
pnts <- mutate(pnts, value = vles)

# Coffee ------------------------------------------------------------------
pnts <- read_csv('//alliancedfs.alliance.cgiar.org/CL9_Coffee_Cocoa2/_coffeeEastAfrica/tbl/occ/occ_swd')
pnts <- mutate(pnts, iso = terra::extract((zone), pnts[,c('Lon', 'Lat')])$sov_a3, .before = Lon)
pnts <- filter(pnts, iso %in% c('UGA', 'KEN'))
pnts <- dplyr::select(pnts, Lon, Lat)

# To extract the values
vles <- terra::extract(ndwl.bsln, pnts[,2:3])[,2]
qntl <- quantile(vles, c(0, 0.5, 0.75, 1), na.rm = T)
mtrx <- matrix(c(0, qntl[2], 1, qntl[2], qntl[3], 2, qntl[3], 100, 3), byrow = T, ncol = 3)
zone <- terra::vect('../data/gpkg/eastafrica.gpkg')

# To classify
clsf <- terra::classify(stck, mtrx, include.lowest = T)
clsf <- terra::crop(clsf, zone) %>% terra::mask(., zone)
dir.create('../data/tif/indices/eaf/coffee/ndwl0')
terra::writeRaster(x = clsf, filename = '../data/tif/indices/eaf/coffee/ndwl0/ndwl0_bsl-ftr_cls.tif', overwrite = TRUE)

# Raster to categorical
clsf <- terra::as.factor(clsf)
clsf.fcal <- terra::focal(clsf, w = 9, fun = modal, na.rm = T)
clsf.fcal <- round(clsf.fcal, 0)
mask <- clsf.fcal[[1]] * 0 + 1
poly <- as.polygons(mask)

factor(clsf, levels = c(1, 2, 3))
clsf.fcal <- terra::crop(clsf.fcal, zone) %>% terra::mask(., zone)

dir.create('../data/tbl/reclassify')
write.csv(mtrx, '../data/tbl/reclassify/ndwl0_values_eaf-coffee.csv', row.names = F)

# To write the raster
terra::writeRaster(x = clsf.fcal, filename = '../data/tif/indices/eaf/ndwl0/ndwl0_bsl-ftr_eaf-cls.tif', overwrite = T)

# =========================================================================
# To make the map ---------------------------------------------------------
# =========================================================================
clsf <- clsf.fcal
tble <- terra::as.data.frame(clsf, xy = T) %>% 
  as_tibble() %>% 
  gather(var, value, -x, -y) %>% 
  mutate(value = factor(value), var = factor(var, levels = c('Baseline', 'Future')))

wrld <- ne_countries(returnclass = 'sf', scale = 50)
gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = value)) + 
  facet_wrap(.~var) +
  scale_fill_manual(values = brewer.pal(n = 3, name = 'BrBG'), labels = c('Low', 'Medium', 'High')) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  geom_sf(data = st_as_sf(zone), fill = NA, col = 'grey20') +
  labs(x = 'Lon', y = 'Lat', fill = 'NDWL0 Class') +
  coord_sf(xlim = ext(zone)[1:2], ylim = ext(zone)[3:4]) +
  ggtitle(label = 'Number of water logging days - Eastafrica', 
          subtitle = '') +
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(size = 6),
        strip.text = element_text(face = 'bold', hjust = 0.5),
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

ggsave(plot = gmap, filename = '../png/maps/indices/ndwl0_bsl-ftr_rcl_eaf-coffee.png', units = 'in', width = 12, height = 6, dpi = 300)

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


