invgamma_u <- function(u, xm, se) {
  a <- (xm^2) / (se^2)
  b <- (se^2) / xm
  qgamma(u, shape = a, rate = b)
}