rdirichlet_u <- function(u_vec, alpha) {
  stopifnot(length(u_vec) == length(alpha))
  g <- qgamma(u_vec, shape = alpha, rate = 1)
  g / sum(g)
}