#' Fitting subsequent maximum spanning trees to a multivariate inverted-Pareto data set 
#'
#' @description 
#' Fits all previous maximum spanning trees to a multivariate inverted-Pareto data set. 
#' The function uses the recursive relationship for conditional distributions.
#'  The properly selected arguments for conditional pair-copulas are passed to the function 'copSelect'.
#' After selecting the class of pair-copula models with the lowest AIC, stores the output of fitted models in the list of MST, including parameter estimates, model classesl (numeric) in the list of MST.
#' After that, the function stores properly defined pseudo-copula observations through functions 'repseudo.SecTree' for the second tree and 'pseudo.Subtree' for subsequent trees in the list of MST.
#' 
#' @param data A \eqn{n\times d} inverted-Pareto sample data set.
#' @param MST A list of maximum spanning trees.
#' @param VineTree A list of fitted vine tree structures and fitted models.
#' @param copfamset Numeric; indicates the class of pair-copula models users consider.
#' @param tree Numeric value of tree levels.
#' @param selectioncrit Character string indicating the model selection criteria, e.g. 'AIC'.
#' @param progress Logical; whether the vine tree selection is reported.
#' @param weights Logical; whether the weights are applied to missing values.
#' @param se Logical; whether standard errors for ML estimators are reported.
#' @param cores Numeric; indicates the number of cores for parallel computing (optional).
#' @param si Numeric; a tuning parameter for mBIC (\eqn{\si_0=0.9}; default).
#' @param effsampsize Numeric; the specified effective sample size for the independence copula (\eqn{n_{D_e}}<10; default).
#' @param tau_threshold Numeric; the specified Kendall's tau value for the independence copula (\eqn{\hat{\tau}_e}<0.05; default).
#' 
#' @return A nested list of MST including the maximum spanning tree and fitted vine models.
#' @export
#'
fit.SubTree <- function(data, MST, VineTree, copfamset, tree, si, selectioncrit, progress, effsampsize, tau_threshold, weights = NA, se = FALSE, cores = 1) 
{
  #VineTree[[tree]] <- MST[[tree]]
  d <- nrow(MST[[tree]]$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)
  for (i in 1:d) { #i=1
    con <- MST[[tree]]$E$nums[i, ]
    PairEdge <- MST[[tree-1]]$E$nums[con, ]
    #PairEdge <- VineTree[[tree-1]]$E$nums[con, ]
    common=ifelse(test = any(PairEdge[1,]==PairEdge[2,1]), yes = PairEdge[2,1], no = PairEdge[2,2])  
    if (PairEdge[1, 1] == common) {
      zr1 <- VineTree[[tree-1]]$E$Copula.CondData.2[con[1]]
      n1 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[1]]
    }else {
      zr1 <- VineTree[[tree-1]]$E$Copula.CondData.1[con[1]]
      n1 <- VineTree[[tree-1]]$E$Copula.CondName.1[con[1]]
    }
    if (PairEdge[2, 1] == common) {
      zr2 <- VineTree[[tree-1]]$E$Copula.CondData.2[con[2]]
      n2 <- VineTree[[tree-1]]$E$Copula.CondName.2[con[2]]
    }else {
      zr2 <- VineTree[[tree-1]]$E$Copula.CondData.1[con[2]]
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
    if (progress == TRUE) 
      message(n1a, " + ", n2a, " --> ", MST[[tree]]$E$names[i])
    pc.data[[i]]$zr1 <- zr1a
    pc.data[[i]]$zr2 <- zr2a
    MST[[tree]]$E$Copula.CondName.1[i] <- n1a
    MST[[tree]]$E$Copula.CondName.2[i] <- n2a
  }

  if (cores > 1){
    cl <- makeCluster()
    lapply <- function(...) parallel::parLapply(cl, 
                                                ...)
  } 
  pc.fits <- lapply(pc.data, copSelect, familyset = copfamset, tree=tree, si=si,
                    selectioncrit = selectioncrit, se = se, effsampsize=effsampsize, tau_threshold=tau_threshold)
  
  for (i in 1:d) {
    MST[[tree]]$E$Copula.param[[i]] <- c(pc.fits[[i]]$par, pc.fits[[i]]$par2)
    MST[[tree]]$E$Copula.type[i] <- pc.fits[[i]]$family
    MST[[tree]]$E$EffectSize[i] <- pc.fits[[i]]$Eff_k
    MST[[tree]]$E$KendallTau[i] <- pc.fits[[i]]$Tau
    MST[[tree]]$E$AIC[i] <- pc.fits[[i]]$AIC
    MST[[tree]]$E$BIC[i] <- pc.fits[[i]]$BIC
    MST[[tree]]$E$mBIC[i] <- pc.fits[[i]]$mBIC
    MST[[tree]]$E$logLik[i] <- pc.fits[[i]]$logLik
    #MST[[tree]]$E$fits[[i]] <- pc.fits[[i]]
  }
  
  ##  Redefine pseudo-observations
  if(tree==2){
    Out <- Condfit.SecondTree(tree = tree, MST = MST, data = data, VineTree = VineTree)
    MST[[tree]]$E$Copula.CondData.1 <- Out$CondData.1
    MST[[tree]]$E$Copula.CondData.2 <- Out$CondData.2
  }
  if(tree > 2){
    Out <- Condfit.SubTree(tree = tree, MST = MST, data = data, VineTree = VineTree)
    MST[[tree]]$E$Copula.CondData.1 <- Out$CondData.1
    MST[[tree]]$E$Copula.CondData.2 <- Out$CondData.2
  }
  MST[[tree]]
}

copSelect <- function (CopulaParams, ...) 
{
  return(pcFamSel(CopulaParams$zr1, CopulaParams$zr2, 
                   ...))
}


