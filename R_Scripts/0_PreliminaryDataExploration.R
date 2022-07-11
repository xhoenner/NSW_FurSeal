rm(list=ls());
library(Hmisc); library(mapdata); library(DAAG); library(dplyr); library(foieGras); library(rgdal); library(raster); library(scales); library(ncdf4); library(fields); library(sf); library(lwgeom); library(spatstat); library(oce); library(gsw);
# source('~/Work/DoEE_ComplianceProject/code/vmsfunction/fn.vms.dist.to.land.R', chdir = TRUE)
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"; mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


## Read location and CTD datasets, along with metadata, bathymetry, and coastline
	## Primary datasets
	ct91.loc <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct91.mdb', tables = 'diag'); ct91.ctd <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct91.mdb', tables = 'ctd')
	ct101.loc <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct101.mdb', tables = 'diag'); ct101.ctd <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct101.mdb', tables = 'ctd')
	ct110.loc <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct110.mdb', tables = 'diag'); ct110.ctd <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct110.mdb', tables = 'ctd')
	metadata <- read.csv('~/Work/WIP/NSW_FurSeal/Data/NSW_FurSeal_TagMetadata.csv')
	metadata$mass <- ifelse(is.na(metadata$estimated_mass), metadata$actual_mass, metadata$estimated_mass)
	
	## Ancillary datasets
	bathymetry <- nc_open("~/Work/Resources/GIS/Bathymetry/GEBCO_2020/gebco_2020_n25.0_s-80.0_w70.0_e179.9.nc");	
	coast <- sf::st_read('~/Work/Resources/GIS/OSM_LandPolygonsComplete/land_polygons.shp'); coast_wgs84 <- st_transform(coast, crs = wgs.84); 
	mapBase <- maps::map("worldHires", fill = T, plot = F)
	
	## Primary datasets preprocessing
	loc <- rbind(ct91.loc, ct101.loc, ct110.loc); ctd <- rbind(ct91.ctd, ct101.ctd, ct110.ctd); # ctd.ssm <- rbind(ct91.ctd.ssm, ct101.ctd.ssm, ct110.ctd.ssm)
	loc <- loc[which(loc$LON > 139 & loc$LON < 163 & loc$LAT > -43 & loc$LAT < -30),]; ctd <- ctd[which(ctd$LON > 139 & ctd$LON < 163 & ctd$LAT > -43 & ctd$LAT < -30),]
	xr <- range(ctd$LON, na.rm = T); yr <- range(ctd$LAT, na.rm = T); expansion <- c(-2, 2);
	loc$ref <- as.character(loc$ref); ctd$ref <- as.character(ctd$ref); # ctd.ssm$ref <- as.character(ctd.ssm$ref)
	loc$D.DATE <- strptime(as.character(loc$D.DATE), '(%m/%d/%y %H:%M:%S)', tz = 'UTC'); ctd$END.DATE <- strptime(as.character(ctd$END.DATE), '(%m/%d/%y %H:%M:%S)', tz = 'UTC'); 
	loc$sattag_program <- sapply(strsplit(loc$ref, '-'), '[[', 1); ctd$sattag_program <- sapply(strsplit(ctd$ref, '-'), '[[', 1); # ctd.ssm$sattag_program <- sapply(strsplit(ctd.ssm$ref, '-'), '[[', 1); 
	sattag_programs <- unique(loc$sattag_program)
	ctd.all <- ctd[which(ctd$TEMP.DBAR != '' & !is.na(ctd$END.DATE)),]; id <- unique(ctd.all$ref); 
	

## Process Argos tracks through foieGras
	## 6 hourly
	for (i in 1:length(id)){
		dat <- loc[which(loc$ref == id[i]), c('ref', 'D.DATE', 'LQ', 'LON', 'LAT')] ## Select an individual
		colnames(dat) <- c('id', 'date', 'lc', 'lon', 'lat')
		dat$date <- as.POSIXct(dat$date)
		dat$lc <- as.character(dat$lc); dat$lc[which(dat$lc == '-1')] <- 'A'; dat$lc[which(dat$lc == '-2')] <- 'B'; dat$lc[which(dat$lc == '-9')] <- 'Z'
		dat$lon <- as.numeric(dat$lon); dat$lat <- as.numeric(dat$lat); dat$id <- as.character(dat$id);
		dat <- dat[which(!is.na(dat$date)),];
		
		dat_ssm_rw <- fit_ssm(dat, model = 'crw', time.step = 6)
		# plot(dat_ssm_rw$ssm[[1]])
		coords <- st_coordinates(st_transform(dat_ssm_rw$ssm[[1]]$predicted, 4326));
		coords.df <- data.frame(dat_ssm_rw$ssm[[1]]$predicted$date, coords); colnames(coords.df)[1] <- 'date'
	
		SSMproc <- data.frame(cbind(dat_ssm_rw$ssm[[1]]$predicted[, c('id', 'date')], coords))
		if (i == 1) {SSMproc.all <- SSMproc} else {SSMproc.all <- rbind(SSMproc.all, SSMproc)}
	}
	SSMproc.all$sattag_program <-  sapply(strsplit(as.character(SSMproc.all$id), '-'), '[[', 1); 

	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsFoieGras2022locations_6hrs.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
		for (i in 1:length(sattag_programs)){
			loc.sel <- loc[which(loc$sattag_program == sattag_programs[i]),]; SSMproc.all.sel <- SSMproc.all[which(SSMproc.all$sattag_program == sattag_programs[i]),]
	
			## Left
				maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
				for (j in 1:length(unique(loc.sel$ref))){
					lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
				}
			## Right			
				maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - FoieGras', sep = ''))
				for (j in 1:length(unique(SSMproc.all.sel$id))){
					lines(SSMproc.all.sel[which(SSMproc.all.sel$id == unique(SSMproc.all.sel$id)[j]), c('X', 'Y')], col = j)
				}
			# pause()
		}
	dev.off()

	## At time of CTD casts
	for (i in 1:length(id)){
		dat <- loc[which(loc$ref == id[i]), c('ref', 'D.DATE', 'LQ', 'LON', 'LAT')]; dat.ctd <- ctd.all[which(ctd.all$ref == id[i]), c('ref', 'END.DATE', 'MAX.DBAR', 'LON', 'LAT')] ## Select an individual
		colnames(dat) <- c('id', 'date', 'lc', 'lon', 'lat'); colnames(dat.ctd) <- c('id', 'date', 'max.dbar', 'lon', 'lat'); dat.ctd <- dat.ctd[order(dat.ctd$date),]
		dat$date <- as.POSIXct(dat$date)
		dat$lc <- as.character(dat$lc); dat$lc[which(dat$lc == '-1')] <- 'A'; dat$lc[which(dat$lc == '-2')] <- 'B'; dat$lc[which(dat$lc == '-9')] <- 'Z'
		dat$lon <- as.numeric(dat$lon); dat$lat <- as.numeric(dat$lat); dat$id <- as.character(dat$id);
		dat <- dat[which(!is.na(dat$date)),];
		
		dat_ssm_rw <- fit_ssm(dat, model = 'crw', time.step = dat.ctd[, 1:2])
		coords <- st_coordinates(st_transform(dat_ssm_rw$ssm[[1]]$predicted, 4326)); coords.df <- data.frame(dat_ssm_rw$ssm[[1]]$predicted$date, coords); colnames(coords.df)[1] <- 'date'
	
		SSMproc <- data.frame(cbind(dat_ssm_rw$ssm[[1]]$predicted[, c('id', 'date')], coords))
		if (i == 1) {SSMproc.all <- SSMproc} else {SSMproc.all <- rbind(SSMproc.all, SSMproc)}
	}
	SSMproc.all$sattag_program <-  sapply(strsplit(as.character(SSMproc.all$id), '-'), '[[', 1); 

	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsFoieGras2022locations_CTDdates.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
		for (i in 1:length(sattag_programs)){
			loc.sel <- loc[which(loc$sattag_program == sattag_programs[i]),]; SSMproc.all.sel <- SSMproc.all[which(SSMproc.all$sattag_program == sattag_programs[i]),]
	
			## Left
				maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
				for (j in 1:length(unique(loc.sel$ref))){
					lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
				}
			## Right			
				maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - FoieGras', sep = ''))
				for (j in 1:length(unique(SSMproc.all.sel$id))){
					lines(SSMproc.all.sel[which(SSMproc.all.sel$id == unique(SSMproc.all.sel$id)[j]), c('X', 'Y')], col = j)
				}
			# pause()
		}
	dev.off()

















	


# ## Before and After foieGras SSM fit
	# tmp <- fn.vms.dist.to.land(dataIN= loc, LonCol=6, LatCol=5);
	# loc$DISTkmLAND <- tmp$DISTkmLAND; loc$PointOnLand <- tmp$PointOnLand
	# loc <- loc[-which(loc$PointOnLand & loc$DISTkmLAND > 10),] ## Remove all points on land, farther than 10 km from the coastline

	
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsFoieGraslocations_Remove10km.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
		for (i in 1:length(sattag_programs)){
			loc.sel <- loc[which(loc$sattag_program == sattag_programs[i]),]; SSMproc.all.sel <- SSMproc.all[which(SSMproc.all$sattag_program == sattag_programs[i]),]
	
			## Left
				map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
				for (j in 1:length(unique(loc.sel$ref))){
					lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
				}
			## Right			
				map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - FoieGras', sep = ''))
				for (j in 1:length(unique(SSMproc.all.sel$id))){
					lines(SSMproc.all.sel[which(SSMproc.all.sel$id == unique(SSMproc.all.sel$id)[j]), c('X', 'Y')], col = j)
				}
			# pause()
		}
	dev.off()
	
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsFoieGraslocations.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
		for (i in 1:length(sattag_programs)){
			loc.sel <- loc[which(loc$sattag_program == sattag_programs[i]),]; ctd.ssm.sel <- ctd.ssm[which(ctd.ssm$sattag_program == sattag_programs[i]),]
	
			## Left
				map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
				for (j in 1:length(unique(loc.sel$ref))){
					lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
				}
			## Right			
				map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - FoieGras CTD', sep = ''))
				for (j in 1:length(unique(ctd.ssm.sel$ref))){
					lines(ctd.ssm.sel[which(ctd.ssm.sel$ref == unique(ctd.ssm.sel$ref)[j]), c('ssm_lon', 'ssm_lat')], col = j)
				}
			# pause()
		}
	dev.off()








































































	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	## FoieGras processed campaigns (see Ian's email dating 17 July 2020)
	ct91.ctd.ssm <- read.csv('~/Work/WIP/NSW_FurSeal/Data/FoieGrasProcessedCampaigns/ct91_hist/ctd_ct91_hist.csv')
	ct101.ctd.ssm <- read.csv('~/Work/WIP/NSW_FurSeal/Data/FoieGrasProcessedCampaigns/ct101_hist/ctd_ct101_hist.csv')
	ct110.ctd.ssm <- read.csv('~/Work/WIP/NSW_FurSeal/Data/FoieGrasProcessedCampaigns/ct110_hist/ctd_ct110_hist.csv')

## Preprocess datasets for subsequent use

	
	## Bathymetry dataset preprocessing
	# lon<- ncvar_get(bathymetry,'x_range')[1]; lat<- ncvar_get(bathymetry,'y_range')[2]; depth<- ncvar_get(bathymetry,'z')
	# lengthlon <- (ncvar_get(bathymetry,'x_range')[2]-ncvar_get(bathymetry,'x_range')[1])*60
	# lengthlat <- (ncvar_get(bathymetry,'y_range')[2]-ncvar_get(bathymetry,'y_range')[1])*60
	# for (i in 2:(lengthlon+1)){lon<-c(lon,lon[1]+(i/60))}
	# for (i in 2:(lengthlat+1)){lat<-c(lat,lat[1]-(i/60))}
	# xmin <- which.min(abs(lon-(round(xr[1],2)-2))); xmax <- which.min(abs(lon-(round(xr[2],2)+2)));	ymin <- which.min(abs(lat-(round(yr[1],2)-2))); ymax <- which.min(abs(lat-(round(yr[2],2)+2)))
	# long <- lon[xmin:xmax]; lati <- lat[ymax:ymin]
	# depthmat<-matrix(depth,ncol=length(lat),nrow=length(lon)); bathy <- depthmat[xmin:xmax,ymax:ymin]
	# dat1 <- list( ); dat1$x <- c( long); dat1$y <- c( lati); dat1$z <- bathy
	# raster_bathy <- raster( dat1$z, xmn = range( dat1[[2]])[1], xmx = range( dat1[[2]])[2], ymn = range( dat1[[1]])[1], ymx = range( dat1[[1]])[2])
	# raster_bathy <- t(raster_bathy);
	# raster_f <- rasterToPoints(raster_bathy); raster_f <- data.frame(raster_f); colnames(raster_f)[3] <- 'Depth'; raster_f$Depth[which(raster_f$Depth > 0)] <- NA

	raster_bathy <- raster("~/Work/Resources/GIS/Bathymetry/GEBCO_2020/gebco_2020_n25.0_s-80.0_w70.0_e179.9.nc")
	raster_bathy <- crop(raster_bathy, extent(xr[1] - 2, xr[2] + 2, yr[1] - 2, yr[2] + 2))
	raster_f <- rasterToPoints(raster_bathy); raster_f <- data.frame(raster_f); colnames(raster_f)[3] <- 'Depth'; raster_f$Depth[which(raster_f$Depth > 0)] <- NA
	
	
	## Map dataset processing
	mapBase <- st_as_sf(mapBase); # coerce it to an "sf" object
	mapBase <- st_make_valid(mapBase); cropMap <- st_crop(mapBase, xmin = min(raster_f$x), xmax = max(raster_f$x), ymin = min(raster_f$y), ymax = max(raster_f$y)) ## Crop
	coast <- crop(shp, extent(min(raster_f$x), max(raster_f$x), min(raster_f$y), max(raster_f$y))); 
	coast_p <- as(coast, 'SpatialPolygons'); coast <- as(coast, 'SpatialLines');
	
	
## Pre-visualise datasets
png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsCTDlocations.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
	for (i in 1:length(sattag_programs)){
		loc.sel <- loc[which(loc$sattag_program == sattag_programs[i]),]; ctd.sel <- ctd[which(ctd$sattag_program == sattag_programs[i]),]

		## Left
			map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
			for (j in 1:length(unique(loc.sel$ref))){
				lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
			}
		## Right			
			map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - CTD', sep = ''))
			for (j in 1:length(unique(ctd.sel$ref))){
				lines(ctd.sel[which(ctd.sel$ref == unique(ctd.sel$ref)[j]), c('LON', 'LAT')], col = j)
			}
		# pause()
	}
dev.off()
	
	
## Extract temperature, salinity and depth values
	ctd.all <- ctd[which(ctd$TEMP.DBAR != '' & !is.na(ctd$END.DATE)),]; id <- unique(ctd.all$ref); ctd.all.backup <- ctd.all
	for (j in 1:length(id)){
		ctd <- ctd.all[which(ctd.all$ref == id[j]),]
		x <- as.numeric(ctd$LON) # Longitude
		y <- as.numeric(ctd$LAT) # Latitude
		date <- ctd$END.DATE
		temp <- psal <- dbar <- NA
		for (s in 1:length(x)){
			lon <- rep(x[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			lat <- rep(y[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			datetime <- rep(date[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			temp <- as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]]) # SST	
			if (ctd$SAL.CORRECTED.VALS[s] == '') {psal <- rep(NA, length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))} else {
				psal <- as.numeric(strsplit(as.character(ctd$SAL.CORRECTED.VALS[s]), ',')[[1]])} # PSAL
			dbar <- as.numeric(strsplit(as.character(ctd$TEMP.DBAR[s]), ',')[[1]]) # SST	
			if(length(temp) != length(psal)) {print(s); break}
			if(s == 1) {lon.all <- lon; lat.all <- lat; date.all <- datetime; temp.all <- temp; psal.all = psal; dbar.all = dbar} else 
				if (s > 1)  {lon.all <- c(lon.all, lon); lat.all <- c(lat.all, lat); date.all <- c(date.all, datetime); temp.all <- c(temp.all, temp); psal.all = c(psal.all, psal); dbar.all = c(dbar.all, dbar)}
		}
		CTD <- data.frame(ref = rep(id[j], length(temp.all)), LON = lon.all, LAT = lat.all, DATE = date.all, TEMP = temp.all, PSAL = psal.all, DBAR = dbar.all)
		if (j == 1) {CTD.ALL <- CTD} else {CTD.ALL <- rbind(CTD.ALL, CTD)}
	}
	CTD.ALL$sattag_program <- sapply(strsplit(as.character(CTD.ALL$ref), '-'), '[[', 1);  CTD.ALL$DayOfYear <- as.numeric(format(CTD.ALL$DATE, '%j'))
	
	
## Extract bathymetry and distance from land for each individual. Potentially boxplot as above
	CTD.ALL$DEPTH_M <- extract(raster_bathy, cbind(CTD.ALL$LON, CTD.ALL$LAT))
	tmp <- fn.vms.dist.to.land(dataIN= CTD.ALL, LonCol=2, LatCol=3);
	CTD.ALL$DISTkmLAND <- tmp$DISTkmLAND
	
## Summarise CTD profile information
	id.summary <- CTD.ALL %>% 
				group_by(ref) %>% 
				summarise(MIN.DATE = as.Date(min(DATE)), MAX.DATE = as.Date(max(DATE)), DEPL.DUR = round(as.numeric(difftime(max(DATE), min(DATE), units = 'days')), 1),
						MIN.TEMP = min(TEMP), MAX.TEMP = max(TEMP), 
						MIN.PSAL = min(PSAL, na.rm = T), MAX.PSAL = max(PSAL, na.rm = T), 
						MEDIAN.DBAR = median(DBAR), MEAN.DBAR = mean(DBAR), MAX.DBAR = max(DBAR),
						NB.PROF = length(unique(DATE)), PROF.FREQ = round(length(unique(DATE))/as.numeric(difftime(max(DATE), min(DATE), units = 'days')), 1))
						
	id.summary <- merge(id.summary, metadata[, c('device_id', 'common_name', 'age_class', 'sex', 'length', 'mass')], by.x = 'ref', by.y = 'device_id', all.y = T)
	id.summary <- id.summary[, c('ref', 'common_name', 'age_class', 'sex', 'length', 'mass', 'MIN.DATE', 'MAX.DATE', 'DEPL.DUR', 'MAX.DBAR', 'NB.PROF', 'PROF.FREQ')]
	id.summary <- id.summary[order(id.summary$ref), ]
	print(id.summary, n = 'Inf')
	

	## Histograms and boxplots
	hist(id.summary$DEPL.DUR, br= 20, xlab = 'Deployment duration (days)', main = ''); abline(v= median(id.summary$DEPL.DUR), col = 'red')
	hist(id.summary$MAX.DBAR, br= 20, xlab = 'Max depth', main = ''); abline(v= median(id.summary$MAX.DBAR), col = 'red')
	
	boxplot.dbar <- CTD.ALL %>% group_by(ref, DATE, sattag_program) %>% summarise(MAX.DBAR = max(DBAR))
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/MaxCTDProfileDepthBoxplot.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
		ggplot(boxplot.dbar, aes(ref, MAX.DBAR, colour = sattag_program)) + geom_boxplot() + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.title = element_blank()) + labs(y= "Maximum CTD profile depth", x = "Device ID")
	dev.off()
	
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/CTDProfileTempBoxplot.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
		ggplot(CTD.ALL, aes(ref, TEMP, colour = sattag_program)) + geom_boxplot() + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.title = element_blank()) + labs(y= "Temperature (deg C)", x = "Device ID")
	dev.off()
	
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/CTDProfileSalinityBoxplot.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
		ggplot(CTD.ALL, aes(ref, PSAL, colour = sattag_program)) + geom_boxplot() + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.title = element_blank()) + labs(y= "Salinity", x = "Device ID")
	dev.off()
	
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/CTDProfileBathymetryBoxplot.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
		ggplot(CTD.ALL, aes(ref, DEPTH_M, colour = sattag_program)) + geom_boxplot() + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.title = element_blank()) + labs(y= "Bathymetry (m)", x = "Device ID")
	dev.off()

	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/CTDProfileDistLandBoxplot.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
		ggplot(CTD.ALL, aes(ref, DISTkmLAND, colour = sattag_program)) + geom_boxplot() + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90), legend.title = element_blank()) + labs(y= "Distance from land (km)", x = "Device ID")
	dev.off()
	

	## CTD profile transmission frequency per day of year
	tag.uplink.summary <- CTD.ALL %>%
					group_by(DayOfYear, ref, sattag_program) %>%
					summarise(NbTransmissions = length(unique(DATE)))

	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/NbCTDProfilesPerDay.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
		ggplot(tag.uplink.summary, aes(DayOfYear, NbTransmissions, color= ref)) + geom_line() + facet_grid(sattag_program ~.)
	dev.off()


## Map of all deployments with bathymetry
	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/StudyAreaTracksMap.png', width = 1920, height = 1080, units = "px", res=92, bg = "white");
	ggplot(cropMap) + geom_raster(data = raster_f, aes(x = x, y = y, fill = Depth)) + labs(y= "Latitude", x = "Longitude") + coord_fixed() + 
		geom_sf() + 
		geom_path(data = ctd.all, aes(x = LON, y = LAT, group = ref, colour = ref)) + 
		theme_bw() + 
		theme(panel.grid = element_line(colour = "transparent"),
	    axis.title = element_text(size = 18),
	    axis.text = element_text(size = 16),
	    legend.text = element_text(size = 14),
	    legend.title = element_text(size = 14))
	dev.off()




## TS plots ODV for animals having ventured far offshore through ODV? Or through my R script developed for Antarctica?
for (i in 1:length(id)){
	sel <- CTD.ALL[which(CTD.ALL$ref == id[i]),]; sel$profile_id <- as.factor(as.character(sel$DATE))
	profile_id <- sel$profile_id
	time <- sel$DATE
	lon <- sel$LON
	lat <- sel$LAT
	pres_adj <- sel$DBAR
	temp_adj <- sel$TEMP
	psal_adj <- sel$PSAL
	sigma0_adj <- gsw_sigma0(psal_adj, temp_adj); # compute sigma0
	ts_hr1 <- as.section(salinity = c(psal_adj), temperature = c(temp_adj), pressure = c(pres_adj), longitude = lon, latitude = lat, station = profile_id);

	# ## Plot CTD profiles
	# plotTS(ts_hr1, eos = 'gsw', type = 'n', levels = seq(round(min(sigma0_adj, na.rm=T),1) - .5, round(max(sigma0_adj, na.rm=T),1) + .5, 0.1),
		# Slim = range(psal_adj) + c(-.1, .1), Tlim = range(temp_adj) + c(-.1, .1))
	# for (j in unique(profile_id)){
		# lines(psal_adj[which(profile_id == j)], temp_adj[which(profile_id == j)], col = 'blue', lwd = 1.6)
	# }
	
	## Plot ocean transects
	GSg <- sectionGrid(ts_hr1)
	# plot(GSg, map.ylim=range(lat) + c(-.5, .5))
	
	# s <- plot(GSg, which="temperature"); ss <- plot(GSg, which="temperature")
	distance <- s[["distance", "byStation"]]
	depth <- s[["station", 1]][["depth"]]
	temperature <- matrix(s[["temperature"]], byrow=TRUE, nrow=length(s[["station"]]))
	salinity <- matrix(ss[["salinity"]], byrow=TRUE, nrow=length(s[["station"]]))
	
	temperature <- temperature[order(distance),]; salinity <- salinity[order(distance),];
	distance <- distance[order(distance)];
	
	png(paste('~/Work/Antarctica/Antarctica_DDU2020/SMRU_SatTagDepl_DDU2020/Outcomes/TSsections_', ind[i], '_Distance_', gsub('-','',Sys.Date()), '.png', sep = ''), width = 4000, height = 2400, units = 'px', res = 300)
	split.screen(c(1,2))
	screen(1)
	plot(GSg, which = 'temperature', xtype = 'distance', ztype = 'image', xlab = 'Distance from deployment site (km)', ylim = c(10,420))
	contour(distance, depth, temperature, add =T) ### Error in contour.default(distance, depth, temperature, add = T) : increasing 'x' and 'y' values expected
	screen(2)
	plot(GSg, which = 'salinity', xtype = 'distance', ztype = 'image', ylab = '', xlab = 'Distance from deployment site (km)', ylim = c(10,420))
	contour(distance, depth, salinity, add =T)
	close.screen(all= T)
	dev.off()
