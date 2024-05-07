#' Redefine pseudo-observations for the conditional pair-copula above the second level tree
#'
#' @description 
#' `Condfit.SubTree()` defines pseudo-observations for the conditional pair-copulas with arguments properly conditioned.
#' For each edge, the conditioning set`condSet` is passed to the function `pseudo.ite()` where we condition all variables in the conditioning set being less than 1.
#' Note that the previously defined pseudo-observations are conditioned on one less variable being less than 1, which are used to estimate pair-copulas in \eqn{T_{i-1}} but not in \eqn{T_i}.
#' 
#' @param tree Numeric; a tree level.
#' @param MST A list of maximum spanning trees.
#' @param data A \eqn{n\times d} inverted-Pareto sample data set.
#' @param VineTree A list of selected vine structures and fitted X-vine models.
#'
#' @return A list of redefined pseudo-observations for the conditional pair-copulas with all variables in the conditioning set being less than 1.
#'  
#' @export
#'
Condfit.SubTree <- function(tree, MST, data, VineTree){
  #when j equals tree=3
  Vine_final=list()
  d=nrow(MST[[tree]]$E$nums)
  for(i in 1:d){ #i=2
    con <- MST[[tree]]$E$nums[i, ]  #MST[[3]]
    PairEdge <- MST[[tree-1]]$E$nums[con, ]
    common=ifelse(test = any(PairEdge[1,]==PairEdge[2,1]), PairEdge[2,1], PairEdge[2,2])
    if(con[1]){  # C4|1;2(lam4|2,lam1|2)
      condSet <- union(MST[[tree-1]]$E$conditionedSet[[con[1]]],MST[[tree-1]]$E$conditioningSet[[con[1]]]) #tree>2  
      #condSet <- union(MST[[tree-1]]$V$conditionedSet[[con[1]]],MST[[tree-1]]$V$conditionedSet[[con[2]]]) 
      Vine_temp <- pseudo.ite(condSet = condSet,treeMinus = tree-1,data = data,MST = MST,VineTree = VineTree)
      #when j=tree, take Vine_temp above whenever the condSet is given.
      if (PairEdge[1, 1] == common) {
        zr1 <- Vine_temp[[tree-1]]$CondData.2[con[1]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[1]]
      }else {
        zr1 <- Vine_temp[[tree-1]]$CondData.1[con[1]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[1]]
      }
      if (PairEdge[2, 1] == common) {
        zr2 <- Vine_temp[[tree-1]]$CondData.2[con[2]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[2]]
      }else {
        zr2 <- Vine_temp[[tree-1]]$CondData.1[con[2]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1 <- as.vector(zr1[[1]])
        zr2 <- as.vector(zr2[[1]])
        n1 <- as.vector(n1[[1]])
        n2 <- as.vector(n2[[1]])
      }else {
        zr1 <- zr1
        zr2 <- zr2
        n1 <- n1
        n2 <- n2
      }
      CoPam  <- do.call(rbind,MST[[tree]]$E$Copula.param)[,1]
      CoFam  <- MST[[tree]]$E$Copula.type
      #CondOn.2 <- CondCop(u1 = zr2,u2 = zr1,par = CoPam[i],family = CoFam[i])  
      CondOn.2 <- BiCopHfunc2(zr2, zr1, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
      #VineTree[[tree]]$E$Copula.CondData.2[i] <- list(CondOn.2)
      #Vine_temp[[tree]]$CondData.2[i] <- list(CondOn.2)
      Vine_final$CondData.2[i] <- list(CondOn.2)
    }
    if(con[2]){  # C1|4;2(lam1|2,lam4|2)
      condSet <- union(MST[[tree-1]]$E$conditionedSet[[con[2]]],MST[[tree-1]]$E$conditioningSet[[con[2]]]) #tree>2  
      #condSet <- union(MST[[tree-1]]$V$conditionedSet[[con[1]]],MST[[tree-1]]$V$conditionedSet[[con[2]]]) 
      Vine_temp <- pseudo.ite(condSet = condSet,treeMinus = tree-1,data = data,MST = MST,VineTree = VineTree)
      if (PairEdge[1, 1] == common) {
        zr1 <- Vine_temp[[tree-1]]$CondData.2[con[1]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[1]]
      }else {
        zr1 <- Vine_temp[[tree-1]]$CondData.1[con[1]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[1]]
      }
      if (PairEdge[2, 1] == common) {
        zr2 <- Vine_temp[[tree-1]]$CondData.2[con[2]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[2]]
      }else {
        zr2 <- Vine_temp[[tree-1]]$CondData.1[con[2]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1 <- as.vector(zr1[[1]])
        zr2 <- as.vector(zr2[[1]])
        n1 <- as.vector(n1[[1]])
        n2 <- as.vector(n2[[1]])
      }else {
        zr1 <- zr1
        zr2 <- zr2
        n1 <- n1
        n2 <- n2
      }
      CoPam  <- do.call(rbind,MST[[tree]]$E$Copula.param)[,1]
      CoFam  <- MST[[tree]]$E$Copula.type
      #CondOn.1 <- CondCop(u1 = zr1,u2 = zr2,par = CoPam[i],family = CoFam[i])
      CondOn.1 <- BiCopHfunc2(zr1, zr2, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
      #VineTree[[tree]]$E$Copula.CondData.1[i] <- list(CondOn.1)
      #Vine_temp[[tree]]$CondData.1[i] <- list(CondOn.1)
      Vine_final$CondData.1[i] <- list(CondOn.1)
    }
  }
  Vine_final
}