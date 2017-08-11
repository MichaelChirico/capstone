library(data.table)

ff = 'output_call.sh'
file.remove(ff)
cat('#Create plots for Shiny\n', file = ff)

for (wk in sprintf('201703%02d', c(1, 8, 15, 22))) {
  
  scores = fread(file.path('scores', wk, 'bo_runs.csv'))
  scores[ , sum(.95^seq_len(.N) * pei),
          keyby = .(delx, dely, alpha, eta, lt, theta,
                    k, l1, l2, kde.bw, kde.lags, kde.win)
          ][which.max(V1), cat(
            do.call(paste, c(list('Rscript show_prediction.r'), .SD)),
            paste0(' ', wk), '\n', sep = '', file = ff, append = TRUE
          ), .SDcols = !'V1']
}
