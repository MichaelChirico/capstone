library(data.table)

NN = 200
params = data.table(delx = runif(NN, 125, 800),
                    dely = runif(NN, 125, 800),
                    eta = runif(NN, .75, 1.25),
                    lt = sample(28, NN, replace = TRUE),
                    theta = runif(NN, 0, pi),
                    features = sample(c(50, 100, 150, 200), NN, replace = TRUE),
                    kde.bw = runif(NN, 125, 800),
                    kde.lags = sample(6, NN, replace = TRUE),
                    kde.win = sample(3:14, NN, replace = TRUE))
params = params[(delx*dely) %between% c(250*250,600*600)]

outf = 'random_search_runs.sh'
cat('#! /bin/bash\n', file = outf)
params[ , {
  cat(paste('echo "', .I, 'of', .N, 
            '"\ntime Rscript predict_with_params.r', 
            delx, dely, eta, lt, theta, features, kde.bw,
            kde.lags, kde.win, 'medical',
            '\nrm tmpf/*'), sep = '\n',
      file = outf, append = TRUE)
}]

# randomly search near optima from random search

rscores = fread('scores/fire_random_search.csv')

outf = 'random_search_runs_local.sh'
cat('#! /bin/bash\n', file = outf)
rscores[ , .(pei = mean(pei)),
         by = .(delx, dely, alpha, eta, lt, theta, k,
                l1, l2, kde.bw, kde.lags, kde.win)
         ][order(-pei)[1:5], 
           .(delx = delx + rnorm(10, sd = 5),
             dely = dely + rnorm(10, sd = 5),
             eta = eta + rnorm(10, sd = .01),
             lt = lt + rnorm(10, sd = 2),
             theta = theta + rnorm(10, sd = pi/64),
             k = pmax(k + sample(-50:50, 10, replace = TRUE), 0),
             kde.bw = kde.bw + rnorm(10, sd = 5),
             kde.lags = pmax(kde.lags + sample(-5:5, 10, replace = TRUE), 0),
             kde.win = pmax(kde.win + sample(-7:7, 10, replace = TRUE), 0)),
           by = seq_len(5L)
          ][ , {
              cat(paste('echo "', .I, 'of 50', 
                        '"\ntime Rscript predict_with_params.r', 
                        delx, dely, eta, lt, theta, k, kde.bw,
                        kde.lags, kde.win, 'fire'), sep = '\n',
                  file = outf, append = TRUE)
             }]
