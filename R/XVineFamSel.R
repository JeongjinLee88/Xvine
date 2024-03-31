#' Model selection for bivariate parametric families along a vine structure
#'
#' @description
#' Assuming that the vine structure is given, `XVineFamSel` chooses bivariate parametric families for each edge in each tree.
#' Specifically, in \eqn{T_1}, we consider the list of four candidate bivariate tail copula densities: the H\eqn{\"{u}}sler-Reiss exponent measure, negative logistic model, logistic model, and Dirichlet model.
#' Similarly, we consider any type of bivariate pair-copula families with a single parameter in \eqn{T_i}, \eqn{i\ge 2}.
#' (e.g. 0=Independence, 1=Gaussian, 3=Clayton, 4=Gumbel, 5=Frank, 6=Joe, 13=Survival Clayton, 14=Survival Gumbel, 16=Survival Joe)
#' We use the `AIC` for model selection and select the bivariate (tail) copula density families with the lowest `AIC` value for each edge in each tree sequentially.
#' Besides the `AIC`, we use two additional criteria to set bivariate parametric copulas to independence copulas: the effective sample size falls below 10 (\eqn{n_{D_e} < 10})
#'  or the empirical Kendall's tau is less than 0.05 (\eqn{\hat{\tau} < 0.05}) for each edge \eqn{e\in E_i}, \eqn{i=2,\ldots,d-1}.
#'  After that, the algorithm estimates the associated parameter(s) for the selected copula family.
#' For more details on the selection of bivariate parametric families along the vine tree sequence, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#' 
#' @param data An \eqn{n\times d} data matrix from the inverted Pareto distribution. In this case, `Rank=T` is not recommended as the rank transformation
#' does not properly work for the data already on uniform scale.
#' If either the original data or samples from multivariate Pareto distribution on Pareto scale are used, then `Rank=T` is required.
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#' @param famset_tc Numeric vector; the model type of bivariate tail copulas.
#' Possible tail copula models are:
#' * 1=Husler-Reiss model
#' * 2=Negative logistic model
#' * 3=Logistic model
#' * 4=Dirichlet model
#' @param famset_cop Numeric vector; the class of bivariate copula families with a single parameter
#' @param selectioncrit Character; the information criteria for model selections (default: `AIC`)
#' @param method Character; parameter estimation method (Default; `mle`)
#' @param weights Logical; weights for missing values (Default; `NA`)
#' @param effsampsize Integer; indicates the effective sample size (Default: 10) for the independence copula
#' @param tau_threshold Numeric; indicates the value of the Kendall's tau (Default: 0.05) for the independence copula
#' @param Rank Logical; whether rank transformation is performed or not (\code{Rank=T}; default).
#' @param qt Numeric; a lower enough threshold used in the rank transformation (e.g. switches from Pareto scale to uniform scale)
#'
#' @return A list of matrix components:
#' 1. `Params`: A \eqn{d\times d} strict upper triangular matrix including ML estimates
#' 2. `famsel`: A \eqn{d\times d} strict upper triangular matrix including the selected families based on the selection criteria
#' 3. `count`: A \eqn{d\times d} strict upper triangular matrix indicating 1 if correctly selected 0 otherwise
#' 4. `DepMeasure`: A \eqn{d\times d} strict upper triangular matrix containing estimates for dependence measures
#' 5. `ll_selected`: An upper triangular matrix that contains log-likelihood values for selected families
#' 6. `AIC_selected`: An upper triangular matrix that contains (averaged) AIC values for selected families
#' 7. `BIC_selected`: An upper triangular matrix that contains (averaged) BIC values for selected families
#' 8. `EffectSamp`: A \eqn{d\times d} strict upper triangular matrix containing effective sample sizes \eqn{n_{D_e}} for selected families
#' @export
#' 
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#'
#' @examples 
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#' 0, 2, 1, 3, 2,
#' 0, 0, 3, 1, 3,
#' 0, 0, 0, 4, 1,
#' 0, 0, 0, 0, 5),5,byrow = TRUE)
#' ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                 0, 0, 2, 2.5, 0.7,
#'                 0, 0, 0, 0.4, -0.3,
#'                 0, 0, 0, 0, 0.1,
#'                 0, 0, 0, 0, 0),5,byrow = TRUE)
#' FamMtx <- matrix(c(0, 1, 2, 3, 4,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 3, 1,
#'                    0, 0, 0, 0, 1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#' # X-vine specification
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' # Pareto random samples
#' Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
#' # Specify the class of families
#' familyset_tc=c(1:4)
#' familyset_cop=c(0,1,3,4,5,6,13,14,16)
#' # Sequential model selection of bivariate parametric families
#' FamSelOut=XVineFamSel(data = Dat_P, Rank = TRUE, qt = 0.05, XVS = XVS
#' , famset_tc = familyset_tc, famset_cop = familyset_cop)
XVineFamSel <- function(data, Rank=TRUE, qt=0.2, XVS, famset_tc, famset_cop, selectioncrit="AIC", method = "mle", effsampsize=10, tau_threshold=0.05, weights = NA)
{
  d <- ncol(data)
  ##  Noe that the 'XVineFamSel' function uses multivariate 'inverted' Pareto samples
  ##  If you directly use samples from the limiting distribution, they must be ones from Pareto distribution with Pareto margin.
  if(Rank){
    data <- ParetoTransRank(data = data, u_quan = qt, scaleType = "U") 
  }
  
  # If data are generated from XVineSim, variables are always put in increasing order 1:d
  # But, if the real data are used, then need to permute columns (relabel the variable indices)
  dior <- diag(XVS$xmat[,,1])
  order.new=varIndexloc(Diag = dior)
  data <- data[,order.new]
  
  Rmtx <- standardStrMtx(XVS$xmat[,,1])
  MaxMat <- createMaxMtx(Matrix = Rmtx)
  fam1 <- XVS$fmat[,,1]
  
  Params <- matrix(0, d, d)
  #Params2 <- matrix(0, d, d)
  DepMeasure <- matrix(0, d, d)
  Vdir <- list()
  Vindir <- list()
  optiout <- list()
  lls <- rep(NA, max(length(famset_cop),length(famset_tc)))
  AICs <- rep(NA, max(length(famset_cop),length(famset_tc)))
  BICs <- rep(NA, max(length(famset_cop),length(famset_tc)))
  famsel <- matrix(0, d, d) 
  count <- matrix(NA, d, d) 
  ll_selected <- matrix(NA, d, d)
  AIC_selected <- matrix(NA, d, d)
  BIC_selected <- matrix(NA, d, d)
  EffectSamp <- matrix(NA, d, d)
  #k=1,2,3,4
  for(k in 1:(d-1)){
    if(k==1){
      for(i in d:2){ #k=1 / i=5,4,3,2
        m <- MaxMat[1, i]
        Z=data[,c(i,m)] # 54
        Z1.2=Z[Z[,2]<=1,] #lambda5.4 given X4<1
        Z2.1=Z[Z[,1]<=1,c(2,1)] #lambda4.5 given X5<1
        Zstar=Z[apply(Z<=1,1,any),]
        #chi <- sum(Z1.2[,1] < 1)/nrow(Z1.2)
        chi1 <- sum(Z1.2[,1] < 1)/nrow(Z1.2) # X1<1 | X2<1
        chi2 <- sum(Z2.1[,1] < 1)/nrow(Z2.1) # X2<1 | X1<1
        chi <- mean(c(chi1,chi2))
        DepMeasure[1,i]=chi  # Empirical Chi
        Effect_k=nrow(Z1.2)+nrow(Z2.1)-nrow(Z[apply(Z<=1,1,all),])
        EffectSamp[1,i]=Effect_k
        
        for(j in 1:length(famset_tc)){
          ml1.2=mleBiTC(ft = LL.BiTC, family=famset_tc[j], data = Z1.2, range = ParRangeTC(famset_tc[j]))
          ml2.1=mleBiTC(ft = LL.BiTC, family=famset_tc[j], data = Z2.1, range = ParRangeTC(famset_tc[j]))
          MLest=(ml1.2+ml2.1)/2
          
          optiout[[j]]=list(family=j, par=MLest, par2=0)
          lls[j] <- -0.5*(LL.BiTC(par = MLest, x = Z1.2, family = famset_tc[j])+LL.BiTC(par = MLest, x = Z2.1, family = famset_tc[j]))
          #lls[j] <- -LL.BiTC(par = MLest, x = Z1.2, family = famset_tc[j])
          #lls[j] <- LL.BiExp(par = MLest, x = Zstar, family = famset_tc[j])
          npars <- 1
          AICs[j] <- -2 * lls[j] + 2 * npars
          BICs[j] <- -2 * lls[j] + log(Effect_k) * npars
        }
        sel <- switch(selectioncrit, logLik = which.max(lls), AIC = which.min(AICs), BIC = which.min(BICs))
        famsel[1,i] <- optiout[[sel]]$family
        Params[1,i]=optiout[[sel]]$par
        ll_selected[1,i] <- lls[sel]
        AIC_selected[1,i] <- AICs[sel]
        BIC_selected[1,i] <- BICs[sel]
        if(sel==fam1[1,i]){
          count[1,i]=1
        }else{
          count[1,i]=0
        }
        # h2(u1|u2) 
        direct <- CondTC(x1 = Z1.2[,1],x2 = Z1.2[,2],par = Params[1,i], family = famsel[1,i])
        indirect <- CondTC(x1 = Z2.1[,1],x2 = Z2.1[,2],par = Params[1,i], family = famsel[1,i])
        Vdir[[i]]=direct
        Vindir[[i]]=indirect
      }
    }
    if(k==2){
      for(i in d:3){ #k=2 / i=5,4,3
        m <- MaxMat[2, i]
        zr1 <- Vdir[[i]] 
        zr2 <- if (m == Rmtx[2, i]) {
          Vdir[[m]]
        }else {
          Vindir[[m]]
        }
        Z=cbind(zr1,zr2)
        ##  Effective sample size for each edge
        Eff_k=nrow(Z)
        EffectSamp[2,i]=Eff_k
        ##  Kendall's tau
        KendallTau=cor(x = zr1,y = zr2,method = 'kendall')
        DepMeasure[2,i]=KendallTau
        ##  Model selection
        for(j in 1:length(famset_cop)){
          #MLest=mle_BiCop(ft = loglik.BiCop, family=famset_cop[j], data = Z, range = ParRangeCop(famset_cop[j]))$ML
          #lls[j] <- -loglik.BiCop(par = MLest, u = Z, family = famset_cop[j])
          MLest=BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$par
          optiout[[j]]=list(family=famset_cop[j], par=MLest, par2=0)
          lls[j] <- BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$logLik
          AICs[j] <- BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$AIC
          BICs[j] <- BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$BIC
        }
        
        if(Eff_k < effsampsize | abs(KendallTau) < tau_threshold){
          sel <- 1
        }else{
          sel <- switch(selectioncrit, logLik = which.max(lls), 
                        AIC = which.min(AICs), BIC = which.min(BICs))
        }
        famsel[2,i] <- optiout[[sel]]$family
        Params[2,i] <- optiout[[sel]]$par
        ll_selected[2,i] <- lls[sel]
        AIC_selected[2,i] <- AICs[sel]
        BIC_selected[2,i] <- BICs[sel]
        count[2,i]=ifelse(famsel[2,i]==fam1[2,i],1,0)
        ##  Store conditional bivariate copulas
        # h2(u1|u2) 
        #direct <- CondCop(u1 = zr1,u2 = zr2,par = Params[2,i],family = famsel[2,i])
        #indirect <- CondCop(u1 = zr2,u2 = zr1,par = Params[2,i],family = famsel[2,i])
        direct <- BiCopHfunc2(zr1, zr2, family = famsel[2,i], par = Params[2,i], par2 = 0, check.pars = FALSE)
        indirect <- BiCopHfunc2(zr2, zr1, family = famsel[2,i], par = Params[2,i], par2 = 0, check.pars = FALSE)
        # h1(u2|u1)
        Vdir[[i]]=direct
        Vindir[[i]]=indirect
      }  
    }
    if(k > 2){
      for(i in d:(k+1)){ #k=3 / i=5,4
        Out=PseudoCop(j = k-1,cind = i,data = data,StrMtx = Rmtx,fam = famsel,MaxMtr = MaxMat,Par = Params)
        m <- MaxMat[k, i]
        zr1 <- Out$Vdir[[i]] 
        zr2 <- if (m == Rmtx[k, i]) {
          Out$Vdir[[m]]
        }else {
          Out$Vindir[[m]] 
        }
        Z=cbind(zr1,zr2)
        ##  Effective sample size
        Eff_k=nrow(Z)
        EffectSamp[k,i]=Eff_k 
        ##  Dependence measure
        KendallTau=cor(x = zr1,y = zr2,method = 'kendall')
        DepMeasure[k,i]=KendallTau
        ##  Model selection
        for(j in 1:length(famset_cop)){
          #MLest=mle_BiCop(ft = loglik.BiCop, family=famset_cop[j], data = Z, range = ParRangeCop(famset_cop[j]))$ML
          #lls[j] <- -loglik.BiCop(par = MLest, u = Z, family = famset_cop[j])
          MLest=BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$par
          optiout[[j]]=list(family=famset_cop[j], par=MLest, par2=0)
          lls[j] <- BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$logLik
          AICs[j] <- BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$AIC
          BICs[j] <- BiCopEst(u1 = zr1,u2 = zr2,family = famset_cop[j],method = "mle",se = FALSE)$BIC
        }
        
        if(Eff_k < effsampsize | abs(KendallTau) < tau_threshold){
          sel <- 1
        }else{
          sel <- switch(selectioncrit, logLik = which.max(lls), 
                        AIC = which.min(AICs), BIC = which.min(BICs))
        }
        famsel[k,i] <- optiout[[sel]]$family
        Params[k,i] <- optiout[[sel]]$par
        ll_selected[k,i] <- lls[sel]
        AIC_selected[k,i] <- AICs[sel]
        BIC_selected[k,i] <- BICs[sel]
        count[k,i]=ifelse(famsel[k,i]==fam1[k,i],1,0)
      }
    }
  }
  return(list("Params"=Params,"famsel"=famsel,"count"=count,"DepMeasure"=DepMeasure,"ll_selected"=ll_selected,"AIC_selected"=AIC_selected,"BIC_selected"=BIC_selected,"EffectSamp"=EffectSamp))
}

