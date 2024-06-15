# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
findMST <- function (g, mode = "RVine", truncated = FALSE) 
{
  if (truncated == FALSE) {
    A <- adjacencyMatrix(g)
    d <- ncol(A)
    if (mode == "RVine") {
      tree <- NULL
      edges <- matrix(NA, d - 1, 2)
      w <- numeric(d - 1)
      i <- 1
      for (k in 1:(d - 1)) {
        tree <- c(tree, i)
        m <- apply(as.matrix(A[, tree]), 2, min)
        a <- apply(as.matrix(A[, tree]), 2, function(x) order(rank(x)))[1, 
        ]
        b <- order(rank(m))[1]
        j <- tree[b]
        i <- a[b]
        edges[k, ] <- c(j, i)
        w[k] <- A[i, j]
        for (t in tree) A[i, t] <- A[t, i] <- Inf
      }
      edges <- t(apply(edges, 1, function(x) sort(x)))
      edges <- edges[order(edges[, 2], edges[, 1]), ]
      E <- g$E$nums
      in.tree <- apply(matrix(edges, ncol = 2), 1, function(x) which((x[1] == 
                                                                        E[, 1]) & (x[2] == E[, 2])))
      MST <- g
      g$E$todel <- rep(TRUE, nrow(E))
      if (any(g$E$todel)) {
        g$E$todel[in.tree] <- FALSE
        MST <- deleteEdges(g)
      }
    }
    else if (mode == "CVine") {
      A <- adjacencyMatrix(g)
      diag(A) <- 0
      sumtaus <- rowSums(A)
      root <- which.min(sumtaus)
      g$E$todel <- !((g$E$nums[, 2] == root) | (g$E$nums[, 
                                                         1] == root))
      MST <- g
      if (any(g$E$todel)) 
        MST <- deleteEdges(g)
    }
    else {
      stop("vine not implemented")
    }
  }
  else {
    MST <- g
    edgesList <- g$E$nums
    uid <- sort(unique(as.vector(g$E$nums)))
    luid <- length(uid)
    if (luid > 2) {
      adjacencyList <- lapply(uid, function(u) c(edgesList[edgesList[, 
                                                                     1] == u, 2], edgesList[edgesList[, 2] == u, 1]))
      dfsorder <- dfs(adjacencyList, 1)
      newEdgesList <- t(apply(dfsorder$E, 1, sort))
      MST$E$todel <- !duplicated(rbind(newEdgesList, edgesList))[-seq(1, 
                                                                      luid - 1)]
    }
    if (any(MST$E$todel)) 
      MST <- deleteEdges(MST)
  }
  MST
}