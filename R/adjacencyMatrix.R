adjacencyMatrix <- function (g) 
{
  d <- length(g$V$names)
  v.all <- cbind(do.call(c, lapply(1:(d - 1), function(i) seq.int(i))), 
                 do.call(c, lapply(1:(d - 1), function(i) rep(i + 1, i))))
  vals <- apply(v.all, 1, set_weight, E = g$E)
  M <- matrix(0, d, d)
  M[upper.tri(M)] <- vals
  M <- M + t(M)
  diag(M) <- Inf
  M
}