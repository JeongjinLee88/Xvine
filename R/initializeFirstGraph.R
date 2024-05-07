InitializeFirstGraph <- function (data, treecrit, weights) 
{
  all.pairs <- combn(1:ncol(data), 2)
  edge.ws <- apply(all.pairs, 2, function(ind) treecrit(data[, 
                                                             ind[1]], data[, ind[2]], weights))
  rel.nobs <- apply(all.pairs, 2, function(ind) mean(!is.na(data[, 
                                                                 ind[1]] + data[, ind[2]])))
  edge.ws <- edge.ws
  W <- diag(ncol(data))
  W[lower.tri(W)] <- edge.ws
  W <- t(W)
  colnames(W) <- rownames(W) <- colnames(data)
  graphFromWeightMatrix(W)
}