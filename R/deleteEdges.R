# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
deleteEdges <- function (g) 
{
  keep <- which(!g$E$todel)
  E <- list(nums = matrix(g$E$nums[keep, ], ncol = 2), names = g$E$names[keep], 
            weights = g$E$weights[keep], conditionedSet = g$E$conditionedSet[keep], 
            conditioningSet = g$E$conditioningSet[keep])
  list(V = g$V, E = E)
}