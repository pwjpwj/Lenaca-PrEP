invbeta_u <- function(u, xm = 0, se = 0) {
  S <- xm * (1 - xm) / (se^2) - 1
  a <- xm * S
  b <- (1 - xm) * S
  qbeta(u, shape1 = a, shape2 = b)
}