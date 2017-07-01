library(data.table)

NN = 150
params = data.table(delx = runif(NN, 125, 800),
                    dely = runif(NN, 125, 800),
                    theta = runif(NN, 0, pi),
                    kde.bw = runif(NN, 125, 800))
params = params[(delx*dely) %between% c(250*250,600*600)]

outf = 'random_search_runs.sh'
cat('#! /bin/bash\n', file = outf)
params[ , {
  cat(paste('echo "', .I, 'of', .N, 
            '"\ntime Rscript predict_with_params.r', 
            delx, dely, '1 1', theta, '0', kde.bw,
            '1 7 fire'), sep = '\n',
      file = outf, append = TRUE)
}]

# randomly search near optima from random search

rscores = fread('scores/fire_kde_only.csv')

outf = 'random_search_runs_local.sh'
cat('#! /bin/bash\n', file = outf)
rscores[ , .(pei = mean(pei), pai = mean(pai)),
         by = .(delx, dely, alpha, theta, l1, l2, kde.bw)
         ][ , .(pei = max(pei), pai = max(pai)),
            by = c('delx', 'dely', 'theta', 'kde.bw')
            ][ , {
              .SD[order(-pei)[1:5], 
                  .(delx = delx + rnorm(10, sd = 10),
                    dely = dely + rnorm(10, sd = 10),
                    theta = theta + rnorm(10, sd = pi/32),
                    kde.bw = kde.bw + rnorm(10, sd = 10)),
                  by = seq_len(5L)
                  ][ , {
                    cat(paste('echo "', .I, 'of 100', 
                              '"\ntime Rscript predict_with_params.r', 
                              delx, dely, '1 1', theta, '0', kde.bw,
                              '1 7 fire'), sep = '\n',
                        file = outf, append = TRUE)
                  }]
              .SD[order(-pai)[1:5], 
                  .(delx = delx + rnorm(10, sd = 10),
                    dely = dely + rnorm(10, sd = 10),
                    theta = theta + rnorm(10, sd = pi/32),
                    kde.bw = kde.bw + rnorm(10, sd = 10)),
                  by = seq_len(5L)
                  ][ , {
                    cat(paste('echo "', .I + 50, 'of 100', 
                              '"\ntime Rscript predict_with_params.r', 
                              delx, dely, '1 1', theta, '0', kde.bw,
                              '1 7 fire'), sep = '\n',
                        file = outf, append = TRUE)
                  }]
              }]
