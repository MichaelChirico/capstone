library(spatstat, quietly = TRUE)
library(splancs, quietly = TRUE)
library(rgeos)
library(rgdal)
library(data.table, warn.conflicts = FALSE, quietly = TRUE)
library(maptools)

#from random.org
set.seed(19775046)

delx = 607.908596915976; dely = 620.742535363564; 
alpha = 0.263157894736842; eta = 1.08594374790609; 
lt = 24.9125476663891; theta = 1.75891299661207;
features = 39; l1 = 0; l2 = 0; kde.bw = 496.257052984183; 
kde.lags = 3; kde.win = 13; call.type = 'fire'

incl_mos = c(10L, 11L, 12L, 1L, 2L, 3L)

aa = delx*dely #forecasted area
lx = eta*250
ly = eta*250

#these files created with get_data notebook
calls = fread(paste0(call.type, '.csv'))
calls[ , date := as.IDate(date)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
# ROTATION ----
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=

#rotation formula, relative to a point (x_0, y_0) that's not origin:
#  [x_0, y_0] + R * [x - x_0, y - y_0]
#  (i.e., rotate the distances from (x_0, y_0) about that point,
#   then offset again by (x_0, y_0))
#  Equivalently (implemented below):
#  (I - R)[x_0, y_0] + R[x, y]
rotate = function(x, y, theta, origin)
  matrix(origin, nrow = length(x), 
         ncol = 2L, byrow = TRUE) %*% (diag(2L) - RT(theta)) + 
  cbind(x, y) %*% RT(theta)
#use the transpose of the rotation matrix to multiply against
#  column vectors of coordinates
RT = function(theta) matrix(c(cos(theta), -sin(theta), 
                              sin(theta), cos(theta)), 
                            nrow = 2L, ncol = 2L)

#use the lower-left corner of data as the origin
#  through which to rotate
#  (side-effect -- the same theta on different
#   call types/horizons will result in different
#   rotated grids, since point0 will likely differ)
point0 = calls[ , c(min(x_lon), min(y_lat))]
calls[ , c('x_lon', 'y_lat') :=
          as.data.table(rotate(x_lon, y_lat, theta, point0))]

#boundary coordinates of seattle,
#  as a matrix (rotated immediately)
seattle = 
  with(fread('data/seattle_coords.csv'),
       rotate(x, y, theta, point0))

#record range here (after rotation), so that
#  we have the same range 
#  after we subset below
#use full boundary range to be sure
#  we eventually cover the output polygon
xrng = range(seattle[ , 1L])
yrng = range(seattle[ , 2L])
    
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
# CREATE GRID ----
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=

getGTindices <- function(gt) {
  # Obtain indices to rearange data from image (eg. result frim pixellate)
  # so that it conforms with data from GridTopology objects (eg. results
  # from using spkernel2d).
  # Input: gt is a grid topology.
  # Returns an index.
  dimx <- gt@cells.dim[1L]
  dimy <- gt@cells.dim[2L]
  c(matrix(seq_len(dimx*dimy), ncol = dimy, byrow = TRUE)[ , dimy:1L])
}

# from create GridTopology corresponding to pixel image used for crime counts
#   (so that we're sure the output from spatstat & sp overlap properly)
grdtop <- as(as.SpatialGridDataFrame.im(
  pixellate(ppp(xrange=xrng, yrange=yrng), eps=c(delx, dely))), "GridTopology")

# index to rearrange rows in pixellate objects
idx.new <- getGTindices(grdtop)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
# CREATE DATA TABLE OF CALLS ----
# 1) aggregate at lag.window level for included periods
# 2) aggregate at forecasting horizon level for included periods
# 3) merge previous results
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=

#Before subsetting, get indices of ever-crime cells
## Per here, these are always sorted by x,y:
##   https://github.com/spatstat/spatstat/issues/37
## NB: this must be done within-loop
##   since it depends on delx & dely
incl_ids = 
  with(calls, as.data.table(pixellate(ppp(
    #pixellate counts dots over each cell,
    #  and appears to do so pretty quickly
    x = x_lon, y = y_lat,
    #ppp complains when given crimes with 
    #  identical coordinates (which is fine
    #  for our purposes), so set check = F
    xrange = xrng, yrange = yrng, check = FALSE),
    #idx.new here converts from spatstat
    #  indexing to GridTopology index order
    #  of the grid cells
    eps = c(delx, dely)))[idx.new]
    #find cells that ever have a crime (value > 0)
  )[value > 0, which = TRUE]

## Associating each event with an interval
##   of width equal to `horizon`, 
##   cascading backwards from the 
##   competition forecasting horizon

# how long is one period for this horizon?
#   NB: these are not great approximations
#       of the horizon lengths, but what is
#       crucial is to line up with the
#       ultimate forecasting horizons,
#       e.g. 1m is March 1 - 31
pd_length = 7L
# how many periods are there in one year for this horizon?
one_year = 52L
# *** TO DO ***
# how many total periods are there in the data?
#   2013/14/15/16/(17), though 17 not used here
n_pds = one_year

#actually easier/quicker to deal with
#  integer representation of the dates
calls[ , date_int := unclass(date)]
#all dates on which an event occurred
unq_crimes = calls[ , unique(date_int)]

march117 = unclass(as.IDate('2017-03-01'))
#all possible period start dates
start = march117 - (seq_len(n_pds) - 1L) * pd_length
#eliminate irrelevant (summer) data
#  and non-testable data after testing dates in 2016
start = start[month(as.IDate(start, origin = '1970-01-01')) %in% incl_mos & 
                start <= march117]
#all period end dates (nonoverlapping with starts)
end = start + pd_length - 1L
#for feeding to foverlaps
windows = data.table(start, end, key = 'start,end')

call_start_map = data.table(date_int = unq_crimes)
call_start_map[ , start_date := 
                   foverlaps(data.table(start = date_int, 
                                        end = date_int),
                             windows)$start]

#finally, add the interval start date associated with
#  each date of occurrence
calls[call_start_map, start_date := i.start_date, on = 'date_int']

#Create the LHS variable -- counts of crimes
#  in each cell, in each interval (indexed by start_date)
X = calls[!is.na(start_date), as.data.table(pixellate(ppp(
  x = x_lon, y = y_lat, xrange = xrng, yrange = yrng, check = FALSE),
  #reorder using GridTopology - im mapping
  eps = c(x = delx, dely)))[idx.new],
  #subset to eliminate never-crime cells
  keyby = start_date][ , I := rowid(start_date)][I %in% incl_ids]

X[ , train_17 := start_date < march117]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
# Add lagged KDE covariates
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=

# create sp object of crimes
prj = CRS("+init=epsg:2926")
to.spdf = function(dt) {
  SpatialPointsDataFrame(
    coords = dt[ , cbind(x_lon, y_lat)],
    data = dt[ , -c('x_lon', 'y_lat')],
    proj4string = prj)
}
calls.sp = to.spdf(calls)
#convert to DT to speed up repeated
#  subsetting by sorting the data
#  to leverage good DT internals
calls.sp@data = setDT(calls.sp@data)
setkey(calls.sp@data, date_int)

#Given a SpatialPointsDF, a start_date,
#  and a lag number, calculate the
#  KDE corresponding to the kde.win
#  days spanning that horizon
compute.kde <- function(pts, start, lag.no) {
  #subset using data.table for speed;
  #  first get the indices via [.data.table
  #  of which crimes are included in that lag number
  idx = pts@data[date_int %between% 
                   #between start - kde.win*lag.no
                   #  and start - (lag.no - 1)*kde.win - 1
                   (start - kde.win*lag.no + c(0, kde.win - 1L)), 
                 which = TRUE]
  #if no events found, just return 0
  #  (spkernel2d handles this with varying degrees of success)
  if (!length(idx)) return(rep(0, length(incl_ids)))
  kde = spkernel2d(pts = pts[idx, ],
                   #immediately exclude never-crime cells
                   poly = seattle, h0 = kde.bw, grd = grdtop)[incl_ids]
}

#start_lag facilitates using within-group lapply
#  in the next command, see below
start_lag = CJ(start = start, lag = seq_len(kde.lags))

#for up through kde.lags total lags to include,
#  compute the associated KDE and add it as a
#  column associated with the grid cell;
#  do this by start_date to complete
#  the lagged-KDE specification on the RHS of the model
RHS = start_lag[ , 
    c(I = list(incl_ids),
      lapply(setNames(lag, paste0('lag', lag)), compute.kde,
             pts = calls.sp, start = .BY$start)), by = start]
    
#join to our main data.table
X = X[RHS, on = c(start_date = 'start', 'I')]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Add Random Fourier Features ----
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#project -- these are the omega * xs
proj = X[ , cbind(x, y, start_date)] %*% 
  (matrix(rt(3L*features, df = 2.5), nrow = 3L)/c(lx, ly, lt))

#now convert to data.table to use fwrite,
#  also moving into VW-parseable format
incl.kde = grep("^lag", setNames(nm = names(X)), value = TRUE)

# check for NAs in the features -- still happens
#   for long, thin cell sizes since we didn't buffer
#   the boundary coordinates enough, and these examples
#   end up returning some NA KDE values because of
#   how spkernel2d decides which cells are within/not
#   the designated boundary coordinates
stopifnot(all(!is.na(X$start_date)))

phi.dt =
  X[ , {
    #some nonsense about how get works in j --
    #  if we define coln_to_vw in global environment,
    #  lapply(incl.kde, coln_to_vw) fails because
    #  get doesn't find the variables.
    #  Probably some workaround, but w/e
    coln_to_vw = function(vn) { 
      V = get(vn)
      #scale up to minimize wasted 0s
      #  (eschewing z-scores somewhat arbitrarily)
      val = V * 10^(abs(round(mean(log10(V[V>0])))))
      #another check on whether we created any
      #  missing KDE values -- you might notice
      #  this happened a lot, much to our chagrin
      if (any(is.nan(val)))
        stop('NaNs detected! Current parameters:',
             paste(args, collapse = '/'))
      sprintf("%s:%.5f", vn, val)
    }
    c(list(v = value, 
           #|kdes and |rff define namespaces,
           #  which we didn't ultimately use
           #  (though could have -- more of
           #   a time constraint than anything)
           l = paste0(I, "_", start_date, "|kdes")), 
      lapply(incl.kde, coln_to_vw),
      list(rff_namespace = '|rff'))
  }]

if (features == 0) phi.dt[ , rff_namespace := NULL]
#about to assign a lot of columns using set --
#  will fail if we don't warn data.table that
#  phi.dt is going to grow substantially
if (features > 500L) invisible(alloc.col(phi.dt, 3L*features))
#create the features
#  previously explored alternative:
#  assign cos/sin projection as matrix:
#  phi = cbind(cos(proj), sin(proj))/sqrt(features)
#  then assign to phi.dt column-wise,
#  but this _appears_ to be slower than implicitly
#  creating this as below by taking sin/cos 
#  simultaneously with assigning to phi.dt.
fkt = 1/sqrt(features)
for (jj in seq_len(features)) {
  pj = proj[ , jj]
  set(phi.dt, j = paste0(c("cos", "sin"), jj), 
      #these are the paired random fourier features;  
      #  truncating at 5 digits to rein in file size
      value = list(sprintf("cos%i:%.5f", jj, fkt*cos(pj)),
                   sprintf("sin%i:%.5f", jj, fkt*sin(pj))))
}
rm(proj)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# WRITE VW FILES ----
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
tdir = "tmpf"

#when we're at the minimum forecast area, we must round up
#  to be sure we don't undershoot; when at the max,
#  we must round down; otherwise, just round
which.round = function(x)
  if (x > 0) {if (x < 1) round else floor} else ceiling


test_idx = !X$train_17

#temporary files -- include lots of
#  potential debugging info in the file name
filename = paste('arws', call.type, delx,
                 dely, eta, lt, theta, features, kde.bw,
                 kde.lags, kde.win, sep = '_')
#write the training data to:
train.vw = paste0(tdir, '/train_', filename)
#write the test data to:
test.vw = paste0(tdir,'/test_', filename)
#simply append .cache suffix to make it easier
#  to track association when debugging
cache = paste0(train.vw, '.cache')
#write the predictions to:
pred.vw = paste0(tdir, '/pred_', filename)

fwrite(phi.dt[!test_idx], train.vw,
       sep = " ", quote = FALSE, col.names = FALSE)
fwrite(phi.dt[test_idx], test.vw,
       sep = " ", quote = FALSE, col.names = FALSE)
rm(phi.dt)

X = X[test_idx]
NN = X[ , sum(value)]

model = tempfile(tmpdir = tdir, pattern = "model")
#train with VW
call.vw = paste('vw --loss_function poisson --l1', l1, 
                '--l2', l2, train.vw, '--cache_file', cache, 
                '--passes 200 -f', model)
system(call.vw, ignore.stderr = TRUE)
system(paste('vw -t -i', model, '-p', pred.vw,
             test.vw, '--loss_function poisson'),
       ignore.stderr = TRUE)

preds =
  fread(pred.vw, sep = " ", header = FALSE, col.names = c("pred", "I_start"))
invisible(file.remove(pred.vw, model, cache, test.vw))
#wrote 2-variable label with _ to fit VW guidelines;
#  now split back to constituents so we can join
preds[ , c("I", "start_date", "I_start") :=
         c(lapply(tstrsplit(I_start, split = "_"), as.integer),
           list(NULL))]

#VW predictions are in logs for poisson regression,
#  so exponentiate to get counts
X[preds, pred.count := exp(i.pred), on = c("I", "start_date")]
rm(preds)

#Calculate all cell rankings, then
#  loop through alpha to include different numbers of them

#Here are the rankings of the cells according to our predictions
X[X[ , .(tot.pred = sum(pred.count)), by = I
    ][order(-tot.pred), .(I, rank = .I)],
  rank := i.rank, on = 'I']

#Here are the TRUE rankings of the cells
X[X[ , .(tot.crimes = sum(value)), by = I
    ][order(-tot.crimes), .(I, true_rank = .I)],
  true_rank := i.true_rank, on = 'I']

n.cells = as.integer(which.round(alpha)(6969600*(1+2*alpha)/aa))
nn = X[rank <= n.cells, sum(value)]
#PEI denominator -- total crimes in the TRUE top n.cells
N_star = X[true_rank <= n.cells, sum(value)]

grdSPDF = SpatialPolygonsDataFrame(
  as.SpatialPolygons.GridTopology(grdtop, proj4string = prj),
  data = data.frame(I = seq_len(prod(grdtop@cells.dim)),
                    row.names = sprintf('g%d', seq_len(prod(grdtop@cells.dim)))), 
  match.ID = FALSE
)

grdSPDF = elide(grdSPDF, rotate = 180/pi * theta,
                center = point0)

grdSPDF@data = setDT(grdSPDF@data)
grdSPDF$pred.count = 0
grdSPDF@data[X, pred.count := i.pred.count, on = 'I']

grdSPDF$hotspot_pred = grdSPDF$I %in% X[rank <= n.cells, I]
grdSPDF$hotspot_actu = grdSPDF$I %in% X[true_rank <= n.cells, I]

seattle = gUnaryUnion(
  readOGR('data', 'Neighborhoods', verbose = FALSE)
)

grdSPDF@data = setDF(grdSPDF@data)
grdSPDF = 
  SpatialPolygonsDataFrame(
    gIntersection(grdSPDF, seattle, byid = TRUE),
    data = grdSPDF@data[gIntersects(grdSPDF, seattle, byid = TRUE), ],
    match.ID = FALSE
  )
proj4string(grdSPDF) = prj

writeOGR(grdSPDF, dsn = '.', layer = 'rff_pred_plot', 
         driver = 'ESRI Shapefile')

spplot(grdSPDF, zcol = 'pred.count', col.regions = viridis(64))