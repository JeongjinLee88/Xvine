BuildNextGraph <- function (oldVineGraph, treecrit, weights = NA, truncated = FALSE) 
{
  d <- nrow(oldVineGraph$E$nums)
  g <- makeFullGraph(d)
  g$V$names <- oldVineGraph$E$names
  g$V$conditionedSet <- oldVineGraph$E$conditionedSet
  g$V$conditioningSet <- oldVineGraph$E$conditioningSet
  
  out <- lapply(seq_len(nrow(g$E$nums)), getEdgeInfo, g = g, 
                oldVineGraph = oldVineGraph, treecrit = treecrit, weights = weights, 
                truncated = truncated)
  g$E$weights <- sapply(out, function(x) x$w)
  g$E$names <- sapply(out, function(x) x$name)
  g$E$conditionedSet <- lapply(out, function(x) x$nedSet)
  g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
  g$E$todel <- sapply(out, function(x) x$todel)
  deleteEdges(g)
}