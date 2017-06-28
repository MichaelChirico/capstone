
download.file('https://data.seattle.gov/api/views/kzjm-xkqj/rows.csv?accessType=DOWNLOAD',
              destfile = 'data/seattle_911.csv')

library(data.table)
library(sp)
library(rgdal)
library(rgeos)

# skip -- malformed second row
calls = fread('data/seattle_911.csv')
NN = nrow(calls)

type_table = fread('data/type_category_match.csv')

calls[ , Type := tolower(gsub('^\\s+|\\s+$', '', gsub('\\s{2,}', ' ', gsub('[[:punct:]]', ' ', Type))))]
calls[type_table, category := i.category, on = 'Type']

calls = calls[!is.na(category)]
cat('Rows eliminated for unmatched category:', NN - nrow(calls), 'of', NN)
NN = nrow(calls)

#eliminate calls missing geospatial tags
calls = calls[!is.na(Longitude)]
cat('\nRows eliminated for missing geospatial information:', NN - nrow(calls), 'of', NN)
NN = nrow(calls)

calls[ , date := as.POSIXct(Datetime, format = '%m/%d/%Y %I:%M:%S %p +0000')]

setnames(calls, c('Latitude', 'Longitude'), c('y_lat', 'x_lon'))
keepcol = c('date', 'y_lat', 'x_lon', 'category')

#a relatively small (<1000) number of rows have mal-formatted Datetime ala
#  YYYY-MM-DDT24:MM:SS+0000 (MM >= 0) probably means after midnight?
calls = calls[!is.na(date), ..keepcol]
cat('\nRows eliminated for malformed timestamp:', NN - nrow(calls), 'of', NN)
NN = nrow(calls)

calls[ , yr := year(date)]
calls[ , wk := week(date)]
#assign pd_no to skip grouping by .(yr, wk) all the time
setorder(calls, yr, wk)
calls[ , pd_no := .GRP, by = .(yr, wk)]
setkey(calls, pd_no)

#eliminate events not within Seattle
callsSP = SpatialPointsDataFrame(
  calls[ , cbind(x_lon, y_lat)], data = calls,
  proj4string = CRS('+init=epsg:4326')
)
#neighborhoods shapefile from
#  https://data.seattle.gov/dataset/Neighborhoods/2mbt-aqqx
seattle = readOGR('data', 'Neighborhoods', verbose = FALSE)

#re-cast coordinates for spatial join
callsSP = spTransform(callsSP, proj4string(seattle))
calls[ , c('x_lon', 'y_lat') := as.data.table(coordinates(callsSP))]

#spatial join, assign neighborhood to each crime
callsSP$nbhd = (callsSP %over% seattle)$OBJECTID

#re-record to data.table
calls[ , nbhd := callsSP$nbhd]
#eliminate outlier crimes not in Seattle boundary
calls = calls[!is.na(nbhd)]
cat('\nRows eliminated for lying outside Seattle:', NN - nrow(calls), 'of', NN)
NN = nrow(calls)

#create data.tables fire & medic
attach(split(calls, by = 'category'))
cat('\nFire Rows:', nrow(calls[category == 'fire']))
cat('\nMedical Rows:', nrow(calls[category == 'medic']))

fwrite(fire, 'fire.csv')
cat('\nFire file written;', round(file.size('fire.csv')/1e6), 'MB')
fwrite(medic, 'medical.csv')
cat('\nMedical file written;', round(file.size('medical.csv')/1e6), 'MB')

library(rgeos)

seattle = gUnaryUnion(gBuffer(
  #buffer to eliminate water features
  #  (2e18 by trial and elimination)
  seattle, width = 2e18*.Machine$double.eps
))

#number of polygon constituents
#  (including holes, islands)
subP = seattle@polygons[[1L]]@Polygons
npoly = length(subP)
out = which.max(sapply(subP, function(p) p@area * (!p@hole)))
#extract the actual boundary polygon to a
#  new object (que feo there's gotta be a better way??)
seattle.boundary = 
  SpatialPolygons(list(Polygons(subP[out], ID = 'boundary')),
                  proj4string = CRS(proj4string(seattle)))

#only really need the coordinates
fwrite(as.data.table(
  seattle.boundary@polygons[[1L]]@Polygons[[1L]]@coords),
  'data/seattle_coords.csv')

plot(seattle, main = 'Spatial Concentration of Fire Emergencies')
invisible(fire[ , smoothScatter(x_lon, y_lat, add = TRUE,
                                colramp = colorRampPalette(c('white', 'red')))])
plot(seattle, add = TRUE)

plot(seattle, main = 'Spatial Concentration of Medical Emergencies')
invisible(medic[ , smoothScatter(x_lon, y_lat, add = TRUE,
                                 colramp = colorRampPalette(c('white', 'blue')))])
plot(seattle, add = TRUE)
