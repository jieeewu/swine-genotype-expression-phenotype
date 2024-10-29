######################################################################
###### 
###### Calculate the proportion of variance of y explained by X
######
######################################################################

calc_effect_size <- function(y, X, Z = NULL, w = NULL, make_intercept = TRUE) {
  if (is.matrix(y)) {
    y <- y[,1]
  }
  
  int_mat <- if (make_intercept) {
    matrix(rep(1, length(y)), ncol = 1, dimnames = list(names(y), NULL))
  } else {
    NULL
  }
  
  fit_alt <- lm.wfit(x = cbind(int_mat, Z, X), y = y, w = w)
  fit_null <- lm.wfit(x = cbind(int_mat, Z), y = y, w = w)
  
  1 - sum((y - fit_alt$fitted.values)^2) / sum((y - fit_null$fitted.values)^2)
}
