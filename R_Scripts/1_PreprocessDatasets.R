## Objective: 
			  # * Preprocess and format datasets
			  # * Fit foieGras SSM to Argos (6 hourly) and CTD cast locations.
			  # * Calculate distance from land
			  # * Calculate distance from SMRU-assigned to SSM-predicted locations
			  # * Extract seafloor depth at each SSM-predicted CTD cast locations 


rm(list=ls());

####################################
## Load up config file and libraries
	Packages <- c('Hmisc', 'sf', 'parallel', 'doSNOW', 'spatstat', 'mapdata', 'foieGras', 'raster')
	pacman::p_load(Packages, character.only = TRUE);
	wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"; mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
	source('~/Work/WIP/NSW_FurSeal/R_Scripts/FindDistance.R', chdir = TRUE)

##################################################
##### Load up and prepare various datasets - START
	## Fur seal datasets
		ct91.loc <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct91.mdb', tables = 'diag'); ct91.ctd <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct91.mdb', tables = 'ctd') ## 9 individuals with CTD profiles
		ct101.loc <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct101.mdb', tables = 'diag'); ct101.ctd <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct101.mdb', tables = 'ctd') ## 14 individuals with CTD profiles
		ct110.loc <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct110.mdb', tables = 'diag'); ct110.ctd <- mdb.get('~/Work/WIP/NSW_FurSeal/Data/ct110.mdb', tables = 'ctd') ## 12 individuals with CTD profiles
	## Coastline dataset
		coast <- sf::st_read('~/Work/Resources/GIS/OSM_LandPolygonsComplete/land_polygons.shp'); coast_wgs84 <- st_transform(coast, crs = wgs.84); 
	## Bathymetry dataset
		raster_bathy <- raster("~/Work/Resources/GIS/Bathymetry/GEBCO_2020/gebco_2020_n25.0_s-80.0_w70.0_e179.9.nc")
##### Load up and prepare various datasets - END
################################################


##########################################
##### Preprocess fur seal datasets - START	
	loc <- rbind(ct91.loc, ct101.loc, ct110.loc); ctd <- rbind(ct91.ctd, ct101.ctd, ct110.ctd); 
	loc <- loc[which(loc$LON > 139 & loc$LON < 163 & loc$LAT > -43 & loc$LAT < -30),]; ctd <- ctd[which(ctd$LON > 139 & ctd$LON < 163 & ctd$LAT > -43 & ctd$LAT < -30),] ## 11253 CTD profiles across 32 individuals
	xr <- range(ctd$LON, na.rm = T); yr <- range(ctd$LAT, na.rm = T); expansion <- c(-2, 2);
	loc$ref <- as.character(loc$ref); ctd$ref <- as.character(ctd$ref); 
	loc$D.DATE <- strptime(as.character(loc$D.DATE), '(%m/%d/%y %H:%M:%S)', tz = 'UTC'); ctd$END.DATE <- strptime(as.character(ctd$END.DATE), '(%m/%d/%y %H:%M:%S)', tz = 'UTC'); 
	loc$sattag_program <- sapply(strsplit(loc$ref, '-'), '[[', 1); ctd$sattag_program <- sapply(strsplit(ctd$ref, '-'), '[[', 1); 
	sattag_programs <- unique(loc$sattag_program)
	ctd.all <- ctd[which(ctd$TEMP.DBAR != '' & ctd$SAL.DBAR != '' & !is.na(ctd$END.DATE)),]; id <- unique(ctd.all$ref); 
	loc <- loc[order(loc$ref, loc$D.DATE),]; ctd.all <- ctd.all[order(ctd.all$ref, ctd.all$END.DATE),]; ## 5656 CTD profiles across 32 individuals (= 50.2 % of all original CTD casts)
##### Preprocess fur seal datasets - END
########################################
	
	
##########################################
##### Calculate distance from land - START
	loc_sf <- st_as_sf(loc, coords = c('LON', 'LAT'), crs = wgs.84);
	bb <- st_bbox(loc_sf); bb[1:2] <- trunc(bb[1:2]) - 2; bb[3:4] <- ceiling(bb[3:4]) + 2; ## Crop coast based on data bounding box
	sf::sf_use_s2(FALSE); coast_wgs84.cropped <- st_crop(st_make_valid(coast_wgs84), bb)
	
	## Transforming to metre projections, multilinestrings, and ppp/psp objects
	coast_mollweide.cropped <- st_transform(coast_wgs84.cropped, crs = mollweide); 
	coast_mollweide.cropped.line <- coast_mollweide.cropped[grep('POLYGON', st_geometry_type(coast_mollweide.cropped)),] %>% st_cast('MULTILINESTRING'); 
	loc_sf_mollweide <- st_transform(loc_sf, crs = mollweide)	
	coast_mollweide.cropped.line.psp <- as.psp(coast_mollweide.cropped.line); 
	
	## Parallelising for distance to coast calculation and point on land
	ids <- sort(unique(loc_sf_mollweide$ref)); pb <- txtProgressBar(max = length(ids), style = 3); progress <- function(n) setTxtProgressBar(pb, n); opts <- list(progress = progress)	
	no_cores <- detectCores(logical = FALSE) - 1; cl <- makeCluster(no_cores); registerDoSNOW(cl); gc(reset = TRUE)  
	tmp <- foreach(ii= 1:length(ids), .combine = rbind, .multicombine = T, .options.snow = opts, .packages = c('spatstat', 'sf')) %dopar% {
		tmp <- loc_sf_mollweide[which(loc_sf_mollweide$ref == ids[ii]),];
		tmp2 <- st_join(tmp, coast_mollweide.cropped, join = st_intersects); tmp2 <- tmp2[which(regexpr("\\.[^\\.]*$", rownames(tmp2)) == -1),]
			return(data.frame(DISTkmLAND = nncross(as.ppp(tmp), coast_mollweide.cropped.line.psp, what = 'dist')/1000, PointOnLand = !is.na(tmp2$FID)));
		}
	close(pb); stopCluster(cl);		
	loc$DISTkmLAND <- tmp$DISTkmLAND; loc$DISTkmLAND[which(tmp$PointOnLand)] <- (-1) * loc$DISTkmLAND[which(tmp$PointOnLand)]; rm(tmp);	
	# plot(loc[, c('LON', 'LAT')], asp = 1, pch = 3, col = alpha('black', .1)); maps::map('worldHires', add = T, fill = T, col = 'grey'); points(loc[which(loc$DISTkmLAND < 0), c('LON', 'LAT')], pch = 3, col = alpha('red', .25))
##### Calculate distance from land - END
########################################


########################################
##### Fit SSM to Argos locations - START
	loc.new <- loc[which(loc$DISTkmLAND > -.5),] ## Only retaining points farther than XX km from the coastline - If using 3 km then losing 1/2 of all Argos locations. Chat to Clive and Mark
	for (i in 1:length(id)){
		dat <- loc.new[which(loc.new$ref == id[i]), c('ref', 'D.DATE', 'LQ', 'LON', 'LAT')] ## Select an individual
		colnames(dat) <- c('id', 'date', 'lc', 'lon', 'lat'); dat$date <- as.POSIXct(dat$date)
		dat$lc <- as.character(dat$lc); dat$lc[which(dat$lc == '-1')] <- 'A'; dat$lc[which(dat$lc == '-2')] <- 'B'; dat$lc[which(dat$lc == '-9')] <- 'Z'
		dat$lon <- as.numeric(dat$lon); dat$lat <- as.numeric(dat$lat); dat$id <- as.character(dat$id);
		dat <- dat[which(!is.na(dat$date)),];
		
		dat_ssm_rw <- fit_ssm(dat, model = 'crw', time.step = 6) ## 6 hourly
		coords <- st_coordinates(st_transform(dat_ssm_rw$ssm[[1]]$predicted, 4326)); coords.df <- data.frame(dat_ssm_rw$ssm[[1]]$predicted$date, coords); colnames(coords.df)[1] <- 'date'
		tmp <- data.frame(cbind(dat_ssm_rw$ssm[[1]]$predicted[, c('id', 'date')], coords))
		if (i == 1) {loc.ssm <- tmp} else {loc.ssm <- rbind(loc.ssm, tmp)}
	}
	loc.ssm$sattag_program <-  sapply(strsplit(as.character(loc.ssm$id), '-'), '[[', 1); 

	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsFoieGras2022locations_6hrsFiltered.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
		for (i in 1:length(sattag_programs)){
			loc.sel <- loc.new[which(loc.new$sattag_program == sattag_programs[i]),]; loc.ssm.sel <- loc.ssm[which(loc.ssm$sattag_program == sattag_programs[i]),]
			## Left
			maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
			for (j in 1:length(unique(loc.sel$ref))){
				lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
			}
			## Right			
			maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - FoieGras', sep = ''))
			for (j in 1:length(unique(loc.ssm.sel$id))){
				lines(loc.ssm.sel[which(loc.ssm.sel$id == unique(loc.ssm.sel$id)[j]), c('X', 'Y')], col = j)
			}
		}
	dev.off()	
##### Fit SSM to Argos locations - END
######################################

	
###############################################
##### Fit SSM to CTD profile timestamps - START
	for (i in 1:length(id)){
		dat <- loc[which(loc$ref == id[i]), c('ref', 'D.DATE', 'LQ', 'LON', 'LAT')]; dat.ctd <- ctd.all[which(ctd.all$ref == id[i]), c('ref', 'END.DATE', 'MAX.DBAR', 'LON', 'LAT')] ## Select an individual
		colnames(dat) <- c('id', 'date', 'lc', 'lon', 'lat'); colnames(dat.ctd) <- c('id', 'date', 'max.dbar', 'lon', 'lat'); dat.ctd <- dat.ctd[order(dat.ctd$date),]
		dat$date <- as.POSIXct(dat$date)
		dat$lc <- as.character(dat$lc); dat$lc[which(dat$lc == '-1')] <- 'A'; dat$lc[which(dat$lc == '-2')] <- 'B'; dat$lc[which(dat$lc == '-9')] <- 'Z'
		dat$lon <- as.numeric(dat$lon); dat$lat <- as.numeric(dat$lat); dat$id <- as.character(dat$id);
		dat <- dat[which(!is.na(dat$date)),];
		
		dat_ssm_rw <- fit_ssm(dat, model = 'crw', time.step = dat.ctd[, 1:2])
		coords <- st_coordinates(st_transform(dat_ssm_rw$ssm[[1]]$predicted, 4326)); coords.df <- data.frame(dat_ssm_rw$ssm[[1]]$predicted$date, coords); colnames(coords.df)[1] <- 'date'
		tmp <- data.frame(cbind(dat_ssm_rw$ssm[[1]]$predicted[, c('id', 'date')], coords))
		if (i == 1) {loc.ssm.ctd <- tmp} else {loc.ssm.ctd <- rbind(loc.ssm.ctd, tmp)}
	}
	loc.ssm.ctd$sattag_program <-  sapply(strsplit(as.character(loc.ssm.ctd$id), '-'), '[[', 1); 

	png(file = '~/Work/WIP/NSW_FurSeal/Outcomes/PreliminaryMaps_ArgosVsFoieGras2022locations_CTDdatesFiltered.png', width = 1920, height = 1920, units = "px", res=92, bg = "white");
		par(mfrow = c(3, 2), mar = c(4, 4, 1.5, .5))
		for (i in 1:length(sattag_programs)){
			loc.sel <- loc[which(loc$sattag_program == sattag_programs[i]),]; loc.ssm.ctd.sel <- loc.ssm.ctd[which(loc.ssm.ctd$sattag_program == sattag_programs[i]),]	
			## Left
			maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - Argos', sep = ''))
			for (j in 1:length(unique(loc.sel$ref))){
				lines(loc.sel[which(loc.sel$ref == unique(loc.sel$ref)[j]), c('LON', 'LAT')], col = j)
			}
			## Right			
			maps::map('worldHires', fill = T, col = 'grey', xlim = xr + expansion, ylim = yr + expansion); map.axes(); title(main = paste(sattag_programs[i], ' - FoieGras', sep = ''))
			for (j in 1:length(unique(loc.ssm.ctd.sel$id))){
				lines(loc.ssm.ctd.sel[which(loc.ssm.ctd.sel$id == unique(loc.ssm.ctd.sel$id)[j]), c('X', 'Y')], col = j)
			}
		}
	dev.off()
	
	loc.ssm.ctd <- loc.ssm.ctd[order(loc.ssm.ctd$id, loc.ssm.ctd$date),]
	ctd.all <- cbind(ctd.all, loc.ssm.ctd[, c('X', 'Y')])
	ctd.all$DISTkm <- FindDistance(ctd.all$LAT, ctd.all$LON, ctd.all$Y, ctd.all$X); length(which(ctd.all$DISTkm > 10)) ## Calculate distance, in km, between CTD and SSM locations
##### Fit SSM to CTD profile timestamps - END
#############################################


#######################################################
##### Calculate distance from land for CTD data - START
	ctd.all_sf <- st_as_sf(ctd.all, coords = c('X', 'Y'), crs = wgs.84);
	ctd.all_sf_mollweide <- st_transform(ctd.all_sf, crs = mollweide);
	
	## Parallelising for distance to coast calculation and point on land
	ids <- sort(unique(ctd.all_sf_mollweide$ref)); pb <- txtProgressBar(max = length(ids), style = 3); progress <- function(n) setTxtProgressBar(pb, n); opts <- list(progress = progress)	
	no_cores <- detectCores(logical = FALSE) - 1; cl <- makeCluster(no_cores); registerDoSNOW(cl); gc(reset = TRUE)  
	tmp <- foreach(ii= 1:length(ids), .combine = rbind, .multicombine = T, .options.snow = opts, .packages = c('spatstat', 'sf')) %dopar% {
		tmp <- ctd.all_sf_mollweide[which(ctd.all_sf_mollweide$ref == ids[ii]),];
		tmp2 <- st_join(tmp, coast_mollweide.cropped, join = st_intersects); tmp2 <- tmp2[which(regexpr("\\.[^\\.]*$", rownames(tmp2)) == -1),]
			return(data.frame(DISTkmLAND = nncross(as.ppp(tmp), coast_mollweide.cropped.line.psp, what = 'dist')/1000, PointOnLand = !is.na(tmp2$FID)));
		}
	close(pb); stopCluster(cl);		
	ctd.all$DISTkmLAND <- tmp$DISTkmLAND; ctd.all$DISTkmLAND[which(tmp$PointOnLand)] <- (-1) * ctd.all$DISTkmLAND[which(tmp$PointOnLand)]; rm(tmp);	
	# plot(ctd.all[, c('X', 'Y')], asp = 1, pch = 3, col = alpha('black', .1)); maps::map('worldHires', add = T, fill = T, col = 'grey'); points(ctd.all[which(ctd.all$DISTkmLAND < 0), c('X', 'Y')], pch = 3, col = alpha('red', .25))
##### Calculate distance from land for CTD data - END
#####################################################


###############################################
## Extract bathymetry for each CTD cast - START
	raster_bathy <- crop(raster_bathy, extent(xr[1] - 2, xr[2] + 2, yr[1] - 2, yr[2] + 2))
	raster_f <- rasterToPoints(raster_bathy); raster_f <- data.frame(raster_f); colnames(raster_f)[3] <- 'Depth'; raster_f$Depth[which(raster_f$Depth > 0)] <- NA
	ctd.all$DEPTH_M <- extract(raster_bathy, cbind(ctd.all$X, ctd.all$Y))
## Extract bathymetry for each CTD cast - END
#############################################


#####################
##### Export datasets
	save(loc.ssm, ctd.all, file= "~/Work/WIP/NSW_FurSeal/Data/1_ArgosCTD_SSM.Rdata");