## Objective: Fit foieGras SSM to Argos and CTD cast locations 


rm(list=ls());

####################################
## Load up config file and libraries
	Packages <- c('gsw', 'oce', 'rLakeAnalyzer')
	pacman::p_load(Packages, character.only = TRUE);

	## Read location and CTD datasets, along with deployment metadata
	load("~/Work/WIP/NSW_FurSeal/Data/1_ArgosCTD_SSM.Rdata")
	metadata <- read.csv('~/Work/WIP/NSW_FurSeal/Data/NSW_FurSeal_TagMetadata.csv')
	metadata$mass <- ifelse(is.na(metadata$estimated_mass), metadata$actual_mass, metadata$estimated_mass)
	

#########################################################
## Extract temperature, salinity and depth values - START
	id <- unique(ctd.all$ref); # ctd.all.backup <- ctd.all
	for (j in 1:length(id)){
		ctd <- ctd.all[which(ctd.all$ref == id[j]),]
		x <- as.numeric(ctd$X) # Longitude
		y <- as.numeric(ctd$Y) # Latitude
		date <- ctd$END.DATE
		DISTkm <- ctd$DISTkm; DISTkmLAND <- ctd$DISTkmLAND; DEPTH_M <- ctd$DEPTH_M
		temp <- psal <- dbar <- NA
		for (s in 1:length(x)){
			prof_id <- rep(s, length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			uprof_id <- rep(rownames(ctd)[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			lon <- rep(x[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			lat <- rep(y[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			datetime <- rep(date[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			distkm <- rep(DISTkm[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			distkmland <- rep(DISTkmLAND[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			depthm <- rep(DEPTH_M[s], length = length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))
			temp <- as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]]) # SST	
			if (ctd$SAL.CORRECTED.VALS[s] == '') {psal <- rep(NA, length(as.numeric(strsplit(as.character(ctd$TEMP.VALS[s]), ',')[[1]])))} else {
				psal <- as.numeric(strsplit(as.character(ctd$SAL.CORRECTED.VALS[s]), ',')[[1]])} # PSAL
			dbar <- as.numeric(strsplit(as.character(ctd$TEMP.DBAR[s]), ',')[[1]]) # SST
			if(length(temp) != length(psal)) {print(s); break}
			if(s == 1) {prof_id.all <- prof_id; uprof_id.all <- uprof_id; lon.all <- lon; lat.all <- lat; date.all <- datetime; distkm.all <- distkm; distkmland.all <- distkmland; depthm.all <- depthm; temp.all <- temp; psal.all = psal; dbar.all = dbar} else 
				if (s > 1)  {prof_id.all <- c(prof_id.all, prof_id); uprof_id.all <- c(uprof_id.all, uprof_id); lon.all <- c(lon.all, lon); lat.all <- c(lat.all, lat); date.all <- c(date.all, datetime); 
								distkm.all <- c(distkm.all, distkm); distkmland.all <- c(distkmland.all, distkmland); depthm.all <- c(depthm.all, depthm); 
								temp.all <- c(temp.all, temp); psal.all = c(psal.all, psal); dbar.all = c(dbar.all, dbar)}
		}
		CTD <- data.frame(ref = rep(id[j], length(temp.all)), unique_profile_id = uprof_id.all, profile_id = prof_id.all, LON = lon.all, LAT = lat.all, DATE = date.all, TEMP = temp.all, PSAL = psal.all, DBAR = dbar.all, 
							DISTkm = distkm.all, DISTkmLAND = distkmland.all, DEPTH_M = depthm.all)
		if (j == 1) {CTD.ALL <- CTD} else {CTD.ALL <- rbind(CTD.ALL, CTD)}
	}
	CTD.ALL$sattag_program <- sapply(strsplit(as.character(CTD.ALL$ref), '-'), '[[', 1);  CTD.ALL$DayOfYear <- as.numeric(format(CTD.ALL$DATE, '%j'))
## Extract temperature, salinity and depth values - END
#######################################################


##################################################
## Remove dubious profiles, ad-hoc process - START
	id2 <- 'ct110-255-12'
	test <- subset(CTD.ALL, ref == id2); 
		test[which(test$PSAL < 35.2),]; # test[which(test$TEMP < 10),]
		test[which(test$unique_profile_id %in% c(9948)),]
	ProfToFilterOut <- c(651, 652, 1029, 1103, 1613, 874, 3999, 4123, 4366, 4421, 4962, 5011, 8771, 8927, 9225, 9351, 9360, 5216, 6130, 6642, 6644, 6711, 6865, 6867, 10784, 10853, 10859, 10874, 10648, 10077, 10220, 9448, 11165)
	CTD.ALL <- CTD.ALL[-which(CTD.ALL$unique_profile_id %in% ProfToFilterOut),] ## Removing 33 profiles manually
## Remove dubious profiles, ad-hoc process - END
################################################

	
###################
## Produce TS plots
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
	
	################
	## Plot profiles
	png(paste0('~/Work/WIP/NSW_FurSeal/Outcomes/CTD/', id[i], '_TempPSALprofiles.png'), width = 4000, height = 2400, units = 'px', res = 300)
		split.screen(c(2,2))
		## Temperature
		screen(1)		
		plot(temp_adj, pres_adj, type = 'n', ylab = 'Pressure (dbar)', xlab = 'Temperature (deg C)', ylim = rev(range(pres_adj)))
		for (j in unique(profile_id)){
			lines(temp_adj[which(profile_id == j)], pres_adj[which(profile_id == j)], col = 'blue', lwd = 1.6)
		}
		## Salinity
		screen(2)
		plot(psal_adj, -pres_adj, type = 'n', ylab = 'Pressure (dbar)', xlab = 'Salinity (psu)', ylim = rev(range(pres_adj)))
		for (j in unique(profile_id)){
			lines(psal_adj[which(profile_id == j)], pres_adj[which(profile_id == j)], col = 'blue', lwd = 1.6)
		}
		## Density
		screen(3)
		plot(sigma0_adj, -pres_adj, type = 'n', ylab = 'Pressure (dbar)', xlab = 'Potential density anomaly (kg.m-3)', ylim = rev(range(pres_adj)))
		for (j in unique(profile_id)){
			lines(sigma0_adj[which(profile_id == j)], pres_adj[which(profile_id == j)], col = 'blue', lwd = 1.6)
		}	
		## TS 
		screen(4)
		plotTS(ts_hr1, eos = 'gsw', type = 'n', levels = seq(round(min(sigma0_adj, na.rm=T),1) - .5, round(max(sigma0_adj, na.rm=T),1) + .5, 0.1),
			Slim = range(psal_adj, na.rm = T) + c(-.1, .1), Tlim = range(temp_adj, na.rm = T) + c(-.1, .1))
		for (j in unique(profile_id)){
			lines(psal_adj[which(profile_id == j)], temp_adj[which(profile_id == j)], col = 'blue', lwd = 1.6)
		}
		close.screen(all = T)
	dev.off()
	
	#######################
	## Plot ocean transects
	GSg <- sectionGrid(ts_hr1)
	# plot(GSg, map.ylim=range(lat) + c(-.5, .5))
	
	s <- plot(GSg, which="temperature"); ss <- plot(GSg, which="temperature"); dev.off()
	distance <- s[["distance", "byStation"]]
	depth <- s[["station", 1]][["depth"]]
	temperature <- matrix(s[["temperature"]], byrow=TRUE, nrow=length(s[["station"]]))
	salinity <- matrix(ss[["salinity"]], byrow=TRUE, nrow=length(s[["station"]]))
	
	temperature <- temperature[order(distance),]; salinity <- salinity[order(distance),];
	distance <- distance[order(distance)];
	
	png(paste0('~/Work/WIP/NSW_FurSeal/Outcomes/CTD/', id[i], '_TSsectionsLongitude.png'), width = 4000, height = 2400, units = 'px', res = 300)
		split.screen(c(2,2))
		screen(1)
		plot(GSg, which = 'temperature', xtype = 'longitude', ztype = 'image', xlab = 'Longitude')
		# contour(distance, depth, temperature, add =T)
		screen(2)
		plot(GSg, which = 'salinity', xtype = 'longitude', ztype = 'image', ylab = '', xlab = 'Longitude')
		# contour(distance, depth, salinity, add =T)
		screen(3)
		plot(GSg, which = 'density', xtype = 'longitude', ztype = 'image', ylab = '', xlab = 'Longitude')
		# contour(distance, depth, salinity, add =T)	
		screen(4)
		plot(lon, lat, pch = 3, asp = 1, col = alpha('black', .1), xlab = 'Longitude', ylab = 'Latitude')
		maps::map('worldHires', add = T, col = 'grey', fill = T)
		lines(lon, lat, cex = .5, lty = 'dashed')
		close.screen(all= T)
	dev.off()
}