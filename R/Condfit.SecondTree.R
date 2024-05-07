#' Define pseudo-observations of the conditional pair-copulas for the second tree
#'
#' @description 
#' `Condfit.SecondTree()` defines pseudo-observations for the conditional pair-copulas with arguments properly conditioned.
#' For each edge, the conditioning set`condSet` is passed to the function `Condfit.FirstTree()` where we condition two variables in the conditioning set being less than 1.
#' Note that the previously defined pseudo-observations are conditioned on one variable being less than 1, which are used to estimate pair-copulas in \eqn{T_2} but not in \eqn{T_3}.
#' 
#' 
#' @param tree Numeric; a tree level.
#' @param data A \eqn{n\times d} inverted-Pareto sample data set.
#' @param MST A list of maximum spanning trees.
#' @param VineTree A list of selected vine structures and fitted X-vine models.
#'
#' @return A list of pseudo-observations such that two variables are less than 1 where two variables are nodes in \eqn{T_2}.
#'  These pseudo-data will be used for estimating pair-copulas in the next tree \eqn{T_3}.
#' @export
#'
Condfit.SecondTree <- function(tree, data, MST, VineTree)
{
  # Goal: MST[[1]] and MST[[2]] with VineTree[[1]] -> define VineTree[[2]]
  VineTemp=list()
  #Vine_temp=lapply(1:tree,function(i)NULL)
  #Vine_temp[[1]] all bivariate tail copula densities in T1
  #Vine_temp[[2]] all pair-copulas in T2 (single conditioning variable)
  #Vine_temp[[3]] all pair-copulas in T3 (two conditioning variables)
  d=nrow(MST[[tree]]$E$nums)
  for(i in 1:d){ #i=123
    con <- MST[[tree]]$E$nums[i, ]  #MST[[2]]
    PairEdge <- MST[[tree-1]]$E$nums[con, ] #VineTree[[1]]
    #PairEdge <- VineTree[[tree-1]]$E$nums[con, ] #VineTree[[1]]
    common=ifelse(test = any(PairEdge[1,]==PairEdge[2,1]), yes = PairEdge[2,1], no = PairEdge[2,2])  
    # Given T1, define T2
    if(con[1]){ 
      ##  Define the conditioning set
      condSet <- MST[[tree]]$V$conditionedSet[[con[1]]]
      ##  Only condition on X1 and X2 not X2 and X4
      Output=Condfit.FirstTree(condSet = condSet,data = data,VineTree = VineTree[[tree-1]])
      CondData11=Output$Pseudo1
      CondData12=Output$Pseudo2
      if (PairEdge[1, 1] == common) {
        zr1 <- CondData12[[con[1]]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[1]]
      }else {
        zr1 <- CondData11[[con[1]]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[1]]
      }
      if (PairEdge[2, 1] == common) {
        zr2 <- CondData12[[con[2]]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[2]]
      }else {
        zr2 <- CondData11[[con[2]]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
        n1a <- as.vector(n1[[1]])
        n2a <- as.vector(n2[[1]])
      }else {
        zr1a <- zr1
        zr2a <- zr2
        n1a <- n1
        n2a <- n2
      }
      CoPam  <- do.call(rbind,MST[[tree]]$E$Copula.param)[,1] # need a VineTree[[2]]
      CoFam  <- MST[[tree]]$E$Copula.type
      #CoPam  <- do.call(rbind,VineTree[[tree]]$E$Copula.param)[,1] # need a VineTree[[2]]
      #CoFam  <- VineTree[[tree]]$E$Copula.type
      #CondOn.2 <- CondCop(u1 = zr2a,u2 = zr1a,par = CoPam[i],family = CoFam[i])  
      CondOn.2 <- BiCopHfunc2(zr2a, zr1a, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
      VineTemp$CondData.2[i] <- list(CondOn.2)
      #VineTree[[tree]]$E$Copula.CondData.2[i] <- list(CondOn.2)
    }
    if(con[2]){
      condSet <- MST[[tree]]$V$conditionedSet[[con[2]]]
      ##  Only condition on X2 and X4 not X1 and X2
      Output2=Condfit.FirstTree(condSet = condSet,data = data,VineTree = VineTree[[tree-1]])
      CondData11=Output2$Pseudo1
      CondData12=Output2$Pseudo2
      if (PairEdge[1, 1] == common) {
        zr1 <- CondData12[[con[1]]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[1]]
      }else {
        zr1 <- CondData11[[con[1]]]
        n1 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[1]]
      }
      if (PairEdge[2, 1] == common) {
        zr2 <- CondData12[[con[2]]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[2]]
      }else {
        zr2 <- CondData11[[con[2]]]
        n2 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
        n1a <- as.vector(n1[[1]])
        n2a <- as.vector(n2[[1]])
      }else {
        zr1a <- zr1
        zr2a <- zr2
        n1a <- n1
        n2a <- n2
      }
      CoPam  <- do.call(rbind,MST[[tree]]$E$Copula.param)[,1] # need a VineTree[[2]]
      CoFam  <- MST[[tree]]$E$Copula.type
      #CondOn.1 <- CondCop(u1 = zr1a,u2 = zr2a,par = CoPam[i],family = CoFam[i])
      CondOn.1 <- BiCopHfunc2(zr1a, zr2a, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
      VineTemp$CondData.1[i] <- list(CondOn.1)
      #VineTree[[tree]]$E$Copula.CondData.1[i] <- list(CondOn.1)
    }
  }  
  VineTemp
}




