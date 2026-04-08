invlognorm_u <- function(u, xm = 0, se = 0) {
  a <- exp(xm + (se^2) / 2)
  b <- ((exp(se^2)) - 1) * exp((2 * xm + se^2))
  qlnorm(u, meanlog = a, sdlog = b)
}