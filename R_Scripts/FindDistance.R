FindDistance <- function(LatsStart, LongsStart,LatsEnd,LongsEnd){
  #this function calculates the great circle distance between two points or vectors of
  #points.  Note that it assumes the earth is spherical, so has an error of ~0.3%
  #formulae are from http://www.movable-type.co.uk/scripts/latlong.html
  
  EarthRadius <- 6371 #in km, so answer is in km
  
  #convert lat, long and bearings from degrees to radians using rad = pi*deg/180
  LatsStart <- (22/7)*LatsStart/180
  LongsStart <- (22/7)*LongsStart/180
  LatsEnd <- (22/7)*LatsEnd/180
  LongsEnd <- (22/7)*LongsEnd/180
  
  DeltaLat <-  LatsEnd - LatsStart
  DeltaLong <-  LongsEnd - LongsStart
  a <- (sin(DeltaLat/2))^2 + (cos(LatsStart)*cos(LatsEnd)*((sin(DeltaLong/2))^2))
  c <- 2*atan2(a^0.5,(1-a)^0.5)
  Distance <- EarthRadius * c
  return(Distance)
}