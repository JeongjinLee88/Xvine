# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
graphFromWeightMatrix <- function (W) 
{
  d <- ncol(W)
  nms <- colnames(W)
  if (is.null(nms)) 
    nms <- paste0("V", 1:d)
  E <- cbind(do.call(c, sapply(1:(d - 1), function(i) seq.int(i))), 
             do.call(c, sapply(1:(d - 1), function(i) rep(i + 1, i))))
  E.names <- apply(E, 1, function(x) paste(nms[x[1]], nms[x[2]], 
                                           sep = ","))
  w <- W[upper.tri(W)]
  list(V = list(names = nms, conditionedSet = NULL, conditioningSet = NULL), 
       E = list(nums = E, names = E.names, weights = w, conditionedSet = lapply(1:nrow(E), 
                                                                                function(i) E[i, ]), conditioningSet = NULL))
}