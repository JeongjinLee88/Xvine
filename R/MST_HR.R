MST_HR <- function(Dat_Pareto,quan)
{
  data_k=data2mpareto(Dat_Pareto,p = 1-quan) # Pareto scale
  Emp_Var=-emp_vario(data_k) # find the whole variogram matrix
  all.pairs <- combn(1:ncol(data_k), 2)  # take a pair of variables
  edge.ws <- apply(all.pairs, 2, function(ind) Emp_Var[ind[1],ind[2]])
  rel.nobs <- apply(all.pairs, 2, function(ind) mean(!is.na(data_k[,ind[1]] + data_k[, ind[2]])))
  W <- diag(ncol(data_k))
  W[lower.tri(W)] <- edge.ws
  W <- t(W)
  colnames(W) <- rownames(W) <- colnames(data_k)
  gg=graphFromWeightMatrix(W)
  MST <- findMST(gg, mode = "RVine")
  return(MST)
}
