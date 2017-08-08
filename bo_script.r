#!/usr/bin/env Rscript
library(spatstat, quietly = TRUE)
library(splancs, quietly = TRUE)
library(rgeos)
library(data.table, warn.conflicts = FALSE, quietly = TRUE)
library(maptools)
library(rBayesianOptimization)

#character(0) if run in live RStudio session
args = commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  # integer version of Mar. 1, 2017
  forecast_start = 17226L
  start_str = '20170301'
  bo_file = 'scores/20170301/bo_runs.csv'
  cat("**********************\n",
      "* TEST PARAMETERS ON *\n",
      "**********************\n")
} else {
  forecast_start = 
    unclass(as.IDate(commandArgs(trailingOnly = TRUE)[1L],
                     format = '%Y%m%d'))
  bo_file = file.path("scores", start_str, "bo_runs.csv")
}

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

  # create sp object of crimes
to.spdf = function(DT) {
  SpatialPointsDataFrame(
    coords = DT[ , cbind(x_lon, y_lat)],
    data = DT[ , -c('x_lon', 'y_lat')],
    proj4string = CRS("+init=epsg:4326"))
}

#when we're at the minimum forecast area, we must round up
#  to be sure we don't undershoot; when at the max,
#  we must round down; otherwise, just round
which.round = function(x)
  if (x > 0) {if (x < 1) round else floor} else ceiling

evaluate_pei = function(delx, dely, eta, lt, theta, k,
                        kde.bw, kde.lags, kde.win) {
  #from random.org
  set.seed(19775046)

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
    spkernel2d(pts = pts[idx, ],
               #immediately exclude never-crime cells
               poly = seattle, h0 = kde.bw, grd = grdtop)[incl_ids]
  }
  
  aa = delx*dely #forecasted area
  lx = eta*250
  ly = eta*250
  
  #these files created with get_data notebook
  calls = fread('fire.csv')
  calls[ , date := as.IDate(date)]
  
  # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
  # ROTATION ----
  # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
    
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
  pd_length = 7L
  # how many periods are there in one year for this horizon?
  one_year = 52L
  half_year = 26L
  yr_length = one_year * pd_length
  
  # count the number of training periods
  n_train = uniqueN(year(calls[date <= forecast_start, date])) - 1L
  
  #actually easier/quicker to deal with
  #  integer representation of the dates
  calls[ , date_int := unclass(date)]
  #all dates on which an event occurred
  unq_crimes = calls[ , unique(date_int)]
  
  #all possible period start dates --
  #  the half_year periods preceding the yr_start for each of
  #  the n_train periods in the _training_ data;
  #  yr_start is built from the input parameter start_str
  start = c(sapply(forecast_start - seq_len(n_train) * yr_length,
                   function(yr_start) yr_start - 
                     (seq_len(half_year) - 1L) * pd_length))
  #windows for feeding to foverlaps to assign a start_date
  #  to each unique event occurrence date
  windows = 
    #all period end dates (nonoverlapping with starts)
    data.table(start, end = start + pd_length - 1L,
               key = 'start,end')
  
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
  
  #Use multiple holdout periods -- one for each possible
  #  year -- to stabilize jumpy prediction validity
  for (ii in seq_len(n_train)) {
    test_start = forecast_start - ii * one_year*pd_length
    X[start_date <= test_start, 
      sprintf('train_%02d', ii) := start_date < test_start]
  }
  
  # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
  # Add lagged KDE covariates
  # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>=
  
  calls.sp = to.spdf(calls)
  #convert to DT to speed up repeated
  #  subsetting by sorting the data
  #  to leverage good DT internals
  calls.sp@data = setDT(calls.sp@data)
  setkey(calls.sp@data, date_int)
    
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
    (matrix(rt(3L*k, df = 2.5), nrow = 3L)/c(lx, ly, lt))
  
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
  
  if (k == 0) phi.dt[ , rff_namespace := NULL]
  #about to assign a lot of columns using set --
  #  will fail if we don't warn data.table that
  #  phi.dt is going to grow substantially
  if (k > 500L) invisible(alloc.col(phi.dt, 3L*k))
  #create the features
  #  previously explored alternative:
  #  assign cos/sin projection as matrix:
  #  phi = cbind(cos(proj), sin(proj))/sqrt(features)
  #  then assign to phi.dt column-wise,
  #  but this _appears_ to be slower than implicitly
  #  creating this as below by taking sin/cos 
  #  simultaneously with assigning to phi.dt.
  fkt = 1/sqrt(k)
  for (jj in seq_len(k)) {
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
  
  #run over these values of alpha (very fast to do within-loop)
  alpha_variations = seq(0, 1, length.out = 20)
  
  #these are all the indicators for training holdouts
  train_variations = grep('^train', names(X), value = TRUE)
  
  #also relatively fast to recycle through VW
  #  runs once we have a data set, so
  #  leverage the fixed cost of data
  #  creation to do so
  l1s = l2s = c(0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3)
  vw_variations = CJ(l1 = l1s, l2 = l2s)
    
  #for storing the output, taking care to
  #  store things in the correct order
  #  relative to our work flow
  scores =
    CJ(train_set = train_variations,
       delx = delx, dely = dely,
       alpha = alpha_variations,
       eta = eta, lt = lt, theta = theta,
       k = k, l1 = l1s, l2 = l2s,
       kde.bw = kde.bw, kde.lags = kde.lags,
       kde.win = kde.win, pei = 0, PEI = NA)
  
  #for easy/clean updating syntax
  setkey(scores, train_set, alpha, l1, l2)
  
  #loop over using each year's holdout test set to calculate PEI/PAI
  for (train in train_variations) {
    #these are the test data
    test_idx = !X[[train]]
    
    #temporary files -- include lots of
    #  potential debugging info in the file name
    filename = paste('arws', train, delx, dely,
                     eta, lt, theta, k, kde.bw,
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
    
    # used to delete phi.dt & subset X here to cut down on
    #   RAM hit. optionally, could split this into two
    #   for loops -- one to write out the files (could then
    #   delete phi.dt) and then another to run VW
    
    # looping over calls to VW
    for (ii in seq_len(nrow(vw_variations))) {
      # print(ii)
      L1 = vw_variations[ii, l1]
      L2 = vw_variations[ii, l2]
      #need a new model file for each run
      model = tempfile(tmpdir = tdir, pattern = "model")
      #train with VW
      call.vw = paste('vw --loss_function poisson --l1', L1, 
                      '--l2', L2, train.vw, '--cache_file', cache, 
                      '--passes 200 -f', model)
      #ingore.stderr since this is where system outputs,
      #  and having all of the VW output (save for debugging)
      #  was just cluttering up the output files
      
      
      system(call.vw, ignore.stderr = TRUE)
      #training data now stored in cache format,
      #  so can delete original (don't need to, but this is a useful
      #  check to force an error if s.t. wrong with cache)
      if (file.exists(train.vw)) invisible(file.remove(train.vw))
      #test with VW
      system(paste('vw -t -i', model, '-p', pred.vw,
                   test.vw, '--loss_function poisson'),
             ignore.stderr = TRUE)
      invisible(file.remove(model))
      
      preds =
        fread(pred.vw, sep = " ", header = FALSE, col.names = c("pred", "I_start"))
      invisible(file.remove(pred.vw))
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
      X[X[test_idx, .(tot.pred = sum(pred.count)), by = I
          ][order(-tot.pred), .(I, rank = .I)],
        rank := i.rank, on = 'I']
      
      #Here are the TRUE rankings of the cells
      X[X[test_idx, .(tot.crimes = sum(value)), by = I
          ][order(-tot.crimes), .(I, true_rank = .I)],
        true_rank := i.true_rank, on = 'I']
      
      # loop over alpha values
      ##this loop can probably be vectorized better...
      for (AA in alpha_variations) {
        #6969600 ft^2 = .25 mi^2 (minimum forecast area);
        #triple this is maximum forecast area
        n.cells = as.integer(which.round(AA)(6969600*(1+2*AA)/aa))
        #n.cells = as.integer(ceiling(6969600/aa))
        
        #how well did we do? lower-case n in the PEI/PAI calculation
        nn = X[rank <= n.cells & test_idx, sum(value)]
        #PEI denominator -- total crimes in the TRUE top n.cells
        N_star = X[true_rank <= n.cells & test_idx, sum(value)]
        
        scores[.(train, AA, L1, L2), pei := nn/N_star]
      }
      #reset ranking for next run
      X[ , rank := NULL]
    }
    invisible(file.remove(cache, test.vw))
  }
  
  #exponentially discounted weighted average mutes
  #  the influence of earlier years in favor of
  #  more recent training sets to allow for trend evolution;
  #  not normalizing since the argmax is unaffected
  best_pei_wtd = scores[ , .(pei = sum(.95^seq_len(.N) * pei)), 
                         keyby = .(alpha, l1, l2)][ , max(pei)]
  fwrite(scores, bo_file, append = file.exists(bo_file))
  return(list(Score = best_pei_wtd, Pred = 0))
}

prior_wk = 
  file.path('scores',
            format(as.IDate(forecast_start - 7L,
                            origin = '1970-01-01'), '%Y%m%d'),
            'bo_runs.csv')

all_params = c('delx', 'dely', 'eta', 'lt', 'theta', 'k',
               'kde.bw', 'kde.lags', 'kde.win')
if (file.exists(bo_file)) {
  current = fread(bo_file)
  current = current[ , .(Value = max(pei)), keyby = all_params]
} else {
  current = data.table()
}

if (file.exists(prior_wk)) {
  old = fread(prior_wk)
  old = old[ , .(Value = max(pei)), keyby = all_params
             #take top 20 parameter sets from prior week
             ][order(-Value)[1:20]]
  #blank the score for these parameters to force
  #  BO to re-calculate scores for these parameters
  old[ , Value := NA_real_]
} else {
  old = data.table()
}

# if current PEI values have been calculated for the old
#   best already, no need to re-calculate them
initial_scores = unique(rbind(current, old), by = all_params)

bounds = list(delx = c(125, 800), dely = c(125, 800),
              eta = c(.5, 1.5), lt = c(3, 30),
              theta = c(0, pi), k = c(0L, 150L),
              kde.bw = c(125, 800), kde.lags = c(1L, 6L),
              kde.win = c(3L, 14L))

fire_opt = BayesianOptimization(evaluate_pei, bounds,
                                initial_scores, 
                                n_iter = 50, verbose = TRUE)
save(fire_opt, paste0('BO_output_', start_str))
