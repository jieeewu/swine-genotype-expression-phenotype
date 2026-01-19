
######################################################################
###### 
###### Calculate empirical hyper priors for 
###### X -> M, X -> Y, and M -> Y
######
######################################################################
source("calc_effect_size.R")

calc_trio_effect_sizes <- function(y, m, X, Z = NULL, Z_y = NULL, Z_m = NULL,
                                    w = NULL, w_y = NULL, w_m = NULL,
                                    make_intercept = TRUE, align_data = TRUE) {
  processed_data <- bmediatR:::process_data(y = y, M = m, X = X,
                                             Z_y = Z_y, Z_M = Z_m,
                                             w_y = w_y, w_M = w_m, 
                                             align_data = align_data,
                                             verbose = FALSE)
  
  y <- processed_data$y
  m <- processed_data$M[,1]
  X <- processed_data$X
  Z_y <- processed_data$Z_y
  Z_m <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_m <- processed_data$w_M
  
  results <- c(
    calc_effect_size(Z = Z_m, y = m, X = X, w = w_m, make_intercept = make_intercept),
    calc_effect_size(Z = Z_y, y = y, X = m, w = w_y, make_intercept = make_intercept),
    calc_effect_size(Z = Z_y, y = y, X = X, w = w_y, make_intercept = make_intercept)
  )
  
  setNames(results, c("x_m", "m_y", "x_y"))
}
