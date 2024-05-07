deleteEdges <- function (g) 
{
  keep <- which(!g$E$todel)
  E <- list(nums = matrix(g$E$nums[keep, ], ncol = 2), names = g$E$names[keep], 
            weights = g$E$weights[keep], conditionedSet = g$E$conditionedSet[keep], 
            conditioningSet = g$E$conditioningSet[keep])
  list(V = g$V, E = E)
}