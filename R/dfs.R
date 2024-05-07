dfs <- function (adjacencyList, v, e = NULL, dfsorder = list()) 
{
  dfsorder$V <- c(dfsorder$V, v)
  dfsorder$E <- rbind(dfsorder$E, e)
  for (u in adjacencyList[[v]]) {
    if (!is.element(u, dfsorder$V)) {
      dfsorder <- dfs(adjacencyList, u, c(u, v), dfsorder)
    }
  }
  return(dfsorder)
}
