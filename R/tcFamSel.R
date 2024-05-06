#' Select the family of bivariate tail copula models for each edge and fit the model
#'
#' @description 
#' Given the family of bivariate tail copula models users consider, `tcFamSel()` selects
#'  the family of bivariate tail copula models with the lowest AIC for a specific edge
#'  , and returns parameter estimates, likelihood value, AIC, and pseud-observations used for the next tree.
#' 
#' @param x1 A numeric vector; indicates the first variable corresponding to the first node of the edge \eqn{e=\{a,b\}} in the first tree.
#' @param x2 A numeric vector; indicates the second variable corresponding to the second node of the edge \eqn{e=\{a,b\}} in the first tree. 
#' @param familyset Numeric vector; indicates the list of bivariate tail copula models.
#' @param selectioncrit Character string; indicates the selection criterion for the family of bivariate tail copula models.
#'
#' @return A list object consisting of the class type, parameter estimates, log-likelihood value, AIC, pseudo-observations
#'  for selected bivariate tail copula models.
#' @export
#'
tcFamSel <- function(x1, x2, si=0.9, familyset, selectioncrit="AIC"){
  
  optiout <- list()
  lls <- rep(Inf, length(familyset))
  AICs <- rep(Inf, length(familyset))
  BICs <- rep(Inf, length(familyset))
  mBICs <- rep(Inf, length(familyset))
  obj <- list()
  
  Z=cbind(x1,x2)
  Zstar=Z[apply(Z<=1,1,any),]
  Eff_k=nrow(Zstar)
  Z1.2=Z[Z[,2]<=1,] #lambda1.2 given X2 < 1
  Z2.1=Z[Z[,1]<=1,c(2,1)] #lambda2.1 given X1 < 1
  
  for(j in 1:length(familyset)){
    ml1.2=mleBiTC(ft = LL.BiTC, family=familyset[j], data = Z1.2, range = ParRangeTC(familyset[j]))
    ml2.1=mleBiTC(ft = LL.BiTC, family=familyset[j], data = Z2.1, range = ParRangeTC(familyset[j]))
    MLest=(ml1.2+ml2.1)/2 # averaged MLE
    
    optiout[[j]]=list(family=j, par=MLest, par2=0)
    
    #lls[j] <- LL.BiTC(par = MLest,x = Zstar,family = familyset[j]) #need to check
    lls[j] <- -0.5*(LL.BiTC(par = MLest,x = Z1.2,family = familyset[j])+LL.BiTC(par = MLest,x = Z2.1,family = familyset[j]))
    npars <- 1
    AICs[j] <- -2 * lls[j] + 2 * npars # 2*(#par) - 2log-likelihood
    BICs[j] <- -2 * lls[j] + log(Eff_k) * npars # log(k)*(#par) - 2log-likelihood
    #mBICs[j] <- -2 * lls[j] + log(Eff_k) * npars -2*log(si) # log(k)*(#par) - 2log-likelihood
    mBICs[j] <- -2 * lls[j] + log(Eff_k) * npars # log(k)*(#par) - 2log-likelihood
  }
  sel <- switch(selectioncrit, logLik = which.max(lls), 
                AIC = which.min(AICs), BIC = which.min(BICs), mBIC = which.min(mBICs))
  obj$family <- optiout[[sel]]$family
  obj$par <- optiout[[sel]]$par
  obj$par2 <- optiout[[sel]]$par2
  obj$logLik <- lls[sel]
  obj$AIC <- AICs[sel]
  obj$BIC <- BICs[sel]
  obj$mBIC <- mBICs[sel]
  obj$Eff_k <- Eff_k
  
  obj$CondOn.1 <- CondTC(x1 = Z1.2[,1],x2 = Z1.2[,2],par = optiout[[sel]]$par,family = optiout[[sel]]$family)
  obj$CondOn.2 <- CondTC(x1 = Z2.1[,1],x2 = Z2.1[,2],par = optiout[[sel]]$par,family = optiout[[sel]]$family)
  obj
}