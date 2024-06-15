# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
getEdgeInfo <- function (i, g, oldVineGraph, treecrit, weights, truncated = FALSE) 
{
  con <- g$E$nums[i, ]
  temp <- oldVineGraph$E$nums[con, ]
  ok <- FALSE
  if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 
                                                        1])) {
    ok <- TRUE
    same <- temp[2, 1]
  }
  else {
    if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 
                                                          2])) {
      ok <- TRUE
      same <- temp[2, 2]
    }
  }
  w <- nedSet <- ningSet <- name <- NA
  todel <- TRUE
  if (ok) {
    l1 <- c(g$V$conditionedSet[[con[1]]], g$V$conditioningSet[[con[1]]])
    l2 <- c(g$V$conditionedSet[[con[2]]], g$V$conditioningSet[[con[2]]])
    nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
    ningSet <- intersect(l1, l2)
    todel <- FALSE
    if (truncated == FALSE) {
      if (temp[1, 1] == same) {
        zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
      }
      else {
        zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
      }
      if (temp[2, 1] == same) {
        zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
      }
      else {
        zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
      }
      else {
        zr1a <- zr1
        zr2a <- zr2
      }
      keine_nas <- !(is.na(zr1a) | is.na(zr2a))
      w <- treecrit(zr1a[keine_nas], zr2a[keine_nas], weights)
      name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
      name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]
      nmdiff <- c(setdiff(name.node1, name.node2), setdiff(name.node2, 
                                                           name.node1))
      nmsect <- intersect(name.node1, name.node2)
      name <- paste(paste(nmdiff, collapse = ","), paste(nmsect, 
                                                         collapse = ","), sep = " ; ")
    }
    else {
      w <- 1
    }
  }
  list(w = w, nedSet = nedSet, ningSet = ningSet, name = name, 
       todel = todel)
}
