# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
makeFullGraph <- function (d) 
{
  E <- cbind(do.call(c, lapply(1:(d - 1), function(i) rep(i, 
                                                          d - i))), do.call(c, lapply(1:(d - 1), function(i) (i + 
                                                                                                                1):d)))
  E <- matrix(E, ncol = 2)
  list(V = list(names = NULL, conditionedSet = NULL, conditioningSet = NULL), 
       E = list(nums = E, names = NULL, weights = NULL, conditionedSet = E, 
                conditioningSet = NULL))
}
