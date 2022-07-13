rm(list = ls())
library(ncdf4); library(mapdata);

source('~/Work/WIP/NSW_FurSeal/R_Scripts/FindDistance.R', chdir = TRUE)
load('~/Work/WIP/NSW_FurSeal/Data/1_ArgosCTD_SSM.Rdata')
setwd('~/Work/WIP/NSW_FurSeal/Data/IMOS_-_Moorings_-_Hourly_time-series_product_source_files')
files <- list.files()

sumnc <- data.frame(matrix(nrow = length(files), ncol = 8))
for (i in 1:length(files)){
	nc <- nc_open(files[i])
	time <- ncvar_get(nc, varid = 'TIME'); time <- as.POSIXct(time * 3600 * 24,origin='1950-01-01 00:00')
	lon <- ncvar_get(nc, varid = 'LONGITUDE')
	lat <- ncvar_get(nc, varid = 'LATITUDE')
	pres <- ncvar_get(nc, varid = 'DEPTH')
	PSAL <- ifelse(length(which(names(nc$var) == 'PSAL')) > 0, TRUE, FALSE)
	TEMP <- ifelse(length(which(names(nc$var) == 'TEMP')) > 0, TRUE, FALSE)
	if(PSAL) {p <- which(!is.na(ncvar_get(nc, varid = 'PSAL')))}
	sumnc[i, 1] <- files[i]
	sumnc[i, 2] <- ifelse(PSAL, as.character(min(time[p])), NA); sumnc[i, 3] <- ifelse(PSAL, as.character(max(time[p])), NA)
	sumnc[i, 4] <- min(lon); sumnc[i, 5] <- min(lat)
	sumnc[i, 6] <- min(pres, na.rm = T); sumnc[i, 7] <- max(pres, na.rm = T)
	sumnc[i, 8] <- ifelse(PSAL, length(p), NA)
}

## Plot all mooring locations
	plot(sumnc[,4:5], pch = 3, col = 'red', asp = 1, type = 'n', xlab = 'Longitude', ylab = 'Latitude')
	map('worldHires', add = T, col = 'grey', fill = T)
	points(sumnc[3:nrow(sumnc),4:5], pch = 3, col = 'red')
	points(sumnc[1:2,4:5], pch = 3, col = 'dark green') ## Only two with salinity data during the time period of interest

## Plot moorings with salinity and seals' CTD casts
plot(sumnc[1:2,4:5], pch = 3, col = 'dark green', asp = 1, xlab = 'Longitude', ylab = 'Latitude', xlim = range(sumnc[1:2,4]) + c(-1, 1), ylim = range(sumnc[1:2,5]) + c(-1, 1))
	map('worldHires', add = T, col = 'grey', fill = T)
	points(ctd.all[which(format(ctd.all$END.DATE, '%Y') == '2012'), c('X', 'Y')], pch = 20, cex = .5, col = 'orange')
	points(ctd.all[which(format(ctd.all$END.DATE, '%Y') == '2013'), c('X', 'Y')], pch = 20, cex = .5, col = 'red')
	points(ctd.all[which(format(ctd.all$END.DATE, '%Y') == '2014'), c('X', 'Y')], pch = 20, cex = .5, col = 'purple')	

ctd.all$DISTkmMooring1 <- FindDistance(ctd.all$Y, ctd.all$X, sumnc[1,5], sumnc[1,4]); length(which(ctd.all$DISTkmMooring1 < 5)) ## Calculate distance, in km, between CTD and mooring locations
ctd.all$DISTkmMooring2 <- FindDistance(ctd.all$Y, ctd.all$X, sumnc[2,5], sumnc[2,4]); length(which(ctd.all$DISTkmMooring2 < 5)) ## Calculate distance, in km, between CTD and mooring locations

ctd.all[which(ctd.all$DISTkmMooring1 < 5),]


ctd.all$DISTkmMooring3 <- FindDistance(ctd.all$Y, ctd.all$X, -34.1192, 151.2267); length(which(ctd.all$DISTkmMooring3 < 5)) ## Calculate distance, in km, between CTD and mooring locations