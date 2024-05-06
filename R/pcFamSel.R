#' Selecting the family of pair-copula models for each edge and fit its model.
#'
#' @description
#' Given the family of pair-copula models users consider, `pcFamSel()` selects
#'  the family of pair-copula models with the lowest AIC for a specific edge, and returns
#'   parameter estimates, likelihood value, AIC, and pseud-observations used for the next tree.
#'   
#' @param u1 A numeric vector, indicating the first argument for the conditional pair-copulas for each edge in each tree.
#' @param u2 A numeric vector, indicating the second argument for the conditional pair-copulas for each edge in each tree.
#' @param familyset A numeric vector, indicating the list of bivariate pair-copulas.
#' @param tree A numeric value indicating the tree level.
#' @param selectioncrit Character string, indicating the selection criteria for bivariate pair copulas.
#' @param se Logical; whether standard errors for ML estimators are reported.
#'
#' @return A list object consisting of the class type, parameter estimates, log-likelihood value, AIC, and pseudo-observations for selected bivariate pair-copulas.
#' @export
#'
pcFamSel <- function(u1, u2, familyset, si=0.9, tree, selectioncrit="AIC", effsampsize=10, tau_threshold=0.05, se=FALSE){
  #u1=pc.data[[1]]$zr1 #if the effective size is too small, Kendall's tau could be -1 or 1
  #u2=pc.data[[1]]$zr2
  optiout <- list()
  lls <- rep(Inf, length(familyset))
  AICs <- rep(Inf, length(familyset))
  BICs <- rep(Inf, length(familyset))
  mBICs <- rep(Inf, length(familyset))
  obj <- list()
  
  Z=cbind(u1,u2)
  Eff_k=nrow(Z) #effective sample sample size = k
  KendallTau=cor(x = u1,y = u2,method = 'kendall')
  ##  If the effective sample size < 10 or tau < 0.05, then the pair-copula is set to the independence copula
  if((Eff_k < effsampsize) | (abs(KendallTau) < tau_threshold)){
    obj$family <- 0
    obj$par <- 0
    obj$par2 <- 0
    obj$logLik <- 0 # for independence copula, ll or AIC = 0 since log(1)=0
    obj$AIC <- 0
    obj$BIC <- 0
    #obj$mBIC <- -2*log(1-si^tree)
    obj$mBIC <- -2*log(1-si^(tree-1))
    obj$Eff_k <- Eff_k
    obj$Tau <- KendallTau
    obj$CondOn.1 <- BiCopHfunc2(u1, u2, family = 0, par = 0, par2 = 0, check.pars = FALSE)
    obj$CondOn.2 <- BiCopHfunc2(u2, u1, family = 0, par = 0, par2 = 0, check.pars = FALSE)
  }else{ ## select the family of pair copula models via AIC
    for(j in 1:length(familyset)){
      MLest=BiCopEst(u1,u2,family = familyset[j],method = "mle",se = se)$par
      optiout[[j]]=list(family=familyset[j], par=MLest, par2=0)
      lls[j] <- BiCopEst(u1,u2,family = familyset[j],method = "mle",se = se)$logLik
      AICs[j] <- BiCopEst(u1,u2,family = familyset[j],method = "mle",se = se)$AIC
      BICs[j] <- BiCopEst(u1,u2,family = familyset[j],method = "mle",se = se)$BIC
      mBICs[j] <- BiCopEst(u1,u2,family = familyset[j],method = "mle",se = se)$BIC-2*(tree-1)*log(si)
      #mBICs[j] <- BiCopEst(u1,u2,family = familyset[j],method = "mle",se = se)$BIC-2*tree*log(si)
    }
    sel <- switch(selectioncrit, logLik = which.max(lls), 
                  AIC = which.min(AICs), BIC = which.min(BICs))
    obj$family <- optiout[[sel]]$family
    obj$par <- optiout[[sel]]$par
    obj$par2 <- optiout[[sel]]$par2
    obj$logLik <- lls[sel] # for independence copula, ll or AIC = 0 since log(1)=0
    obj$AIC <- AICs[sel]
    obj$BIC <- BICs[sel]
    obj$mBIC <- mBICs[sel]
    obj$Eff_k <- Eff_k
    obj$Tau <- KendallTau
    obj$CondOn.1 <- BiCopHfunc2(u1, u2, family = optiout[[sel]]$family, par = optiout[[sel]]$par, par2 = 0, check.pars = FALSE)
    obj$CondOn.2 <- BiCopHfunc2(u2, u1, family = optiout[[sel]]$family, par = optiout[[sel]]$par, par2 = 0, check.pars = FALSE)
  }
  #obj$CondOn.1 <- CondCop(u1 = u1,u2 = u2,par = optiout[[sel]]$par,family = optiout[[sel]]$family)
  #obj$CondOn.2 <- CondCop(u1 = u2,u2 = u1,par = optiout[[sel]]$par,family = optiout[[sel]]$family)  
  obj
}