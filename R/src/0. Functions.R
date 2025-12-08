
#function to calcualte standard error
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

#define function for calculating cv
cv <- function(x) {
  sd(x) / mean(x)
}

#Proportional Variability (PV)
p_cv <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  z_values <- outer(x, x, function(a, b) 1 - pmin(a, b) / pmax(a, b))
  P_CV <- 2 * sum(z_values[upper.tri(z_values)], na.rm = TRUE) / (n * (n - 1))
  return(P_CV)
}

#Kvalseths CV
k_cv <- function(x) {
  if (length(x) < 2) return(NA)
  CV <- sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  k_CV <- sqrt( (CV^2) / (1 + (CV^2)) )
  return(k_CV)
  
}