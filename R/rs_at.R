rs_at <- function(QQ, idx, t) {
  if (length(idx) == 1L) QQ[, idx, t]
  else rowSums(QQ[, idx, t, drop = FALSE])
}