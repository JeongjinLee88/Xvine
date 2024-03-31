#' Sequential parameter estimation for X-Vine models
#'
#' @description
#' `XVineSeqEst` estimates parameters associated with bivariate (tail) copula densities tree by tree sequentially from \eqn{T_1} to \eqn{T_{d-1}}.
#' Information about the X-Vine specification is contained in the object from \code{XVineSpec()}.
#' For more details on sequential parameter estimation, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#' 
#' @param data An \eqn{n\times d} data matrix from the inverted Pareto distribution. In this case, `Rank=T` is not recommended as the rank transformation
#' does not properly work for the data already on uniform scale.
#' If either the original data or samples from multivariate Pareto distribution on Pareto scale are used, then `Rank=T` is required.
#' 
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#' @param method Character; indicates the parameter estimation method (\code{method="mle"}; default).
#' @param Rank Logical; whether rank transformation is performed or not (\code{Rank=T}; default). It switches from Pareto scale to uniform scale.
#' @param qt Numeric; a lower threshold for the rank transformation. 
#'
#' @return 
#' A list of four components: the matrix of sequentially estimated parameters are stored in `Params`,
#' the estimated dependence measures are stored in the matrix `DepMeasure` where the first row has empirical Chi values and subsequent trees have empirical Kendall's Tau,
#' , the effective sample size for each edge in the vine tree is stored in the matrix `EffectSamp`,
#' and the log-likelihood value for each edge in each tree is stored in the matrix `logLik`.
#' @export
#'
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
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
#' Dat_P=ParetoSim(n = 5000, XVS = XVS) # Pareto scale
#' # Sequential parameter estimation
#' SeqEstOut=XVineSeqEst(data = Dat_P, Rank = TRUE, qt = 0.05, XVS=XVS, method = 'mle')
XVineSeqEst <- function(data, Rank=TRUE, qt=0.2, XVS, method = "mle")
{
  d <- ncol(data)
  ##  Noe that the 'XVineSeqEst' function uses multivariate 'inverted' Pareto samples
  ##  If you directly use samples from the limiting distribution, they must be ones from Pareto distribution with Pareto margin.
  if(Rank){
    data <- ParetoTransRank(data = data, u_quan = qt, scaleType = "U") 
  }
  
  # If data are generated from XVineSim, variables are always put in increasing order 1:d
  # But, if the real data are used, then need to permute columns (relabel the variable indices)
  dior <- diag(XVS$xmat[,,1])
  order.new=varIndexloc(Diag = dior)
  data <- data[,order.new]
  
  M <- standardStrMtx(XVS$xmat[,,1]) # reconstruct the structure matrix with diagonal elements in increasing order
  MaxMat <- createMaxMtx(Matrix = M) # define a max-matrix
  fam1 <- XVS$fmat[,,1] # store the family matrix corresponding to the structure matrix 'M'
  
  Params <- matrix(0, d, d)
  #Params2 <- matrix(0, d, d) # to be used when copula families with two parameters considered
  DepMeasure <- matrix(0, d, d)
  EffectSamp <- matrix(0, d, d)
  logLiks <- matrix(0, d, d)
  Vdir <- list() # store the pseudo data corresponding to the first argument
  Vindir <- list() # store the pseudo data corresponding to the second argument
  
  for(k in 1:(d-1)){
    if(k==1){
      for(i in d:2){ #k=1 / i=5,4,3,2
        m <- MaxMat[1, i]
        Z=data[,c(i,m)] # ex: (i,m)=(5,4)
        Z1.2=Z[Z[,2]<1,] #r(x5,x4) given x4<1. To coincide with the specified percentile, the equality is included.
        Z2.1=Z[Z[,1]<1,c(2,1)] #r(x4,x5) given X5<1
        ml1.2=mleBiTC(ft = LL.BiTC, family=fam1[k,i], data = Z1.2, range = ParRangeTC(fam1[k,i]))
        ml2.1=mleBiTC(ft = LL.BiTC, family=fam1[k,i], data = Z2.1, range = ParRangeTC(fam1[k,i]))
        MLest=(ml1.2+ml2.1)/2 # averaged MLE
        
        #chi <- sum(Z1.2[,1] < 1)/nrow(Z1.2)  # under the rank trans, chi1=chi2
        chi1 <- sum(Z1.2[,1] < 1)/nrow(Z1.2)
        chi2 <- sum(Z2.1[,1] < 1)/nrow(Z2.1)
        chi <- mean(c(chi1,chi2))
        
        # Lambda(x1|x2)
        direct <- CondTC(x1 = Z1.2[,1],x2 = Z1.2[,2],par = MLest,family = fam1[k,i])
        # Lambda(x2|x1)
        indirect <- CondTC(x1 = Z2.1[,1],x2 = Z2.1[,2],par = MLest,family = fam1[k,i])
        Vdir[[i]]=direct
        Vindir[[i]]=indirect
        Params[k,i]=MLest
        DepMeasure[k,i]=chi
        #EffectSamp[k,i]=nrow(Z[apply(Z<1,1,any),])
        EffectSamp[k,i]=nrow(Z1.2)+nrow(Z2.1)-nrow(Z[apply(Z<=1,1,all),])
        logLiks[k,i]=LL.BiTC(par = MLest,x = Z1.2,family = fam1[k,i])
      }
    }
    if(k==2){
      for(i in d:3){ #k=2 / i=5,4,3
        m <- MaxMat[2, i]
        zr1 <- Vdir[[i]] 
        zr2 <- if (m == M[2, i]) {
          Vdir[[m]]
        }else {
          Vindir[[m]]
        }
        Z=cbind(zr1,zr2)
        # ML estimate
        MLest=BiCopEst(u1 = zr1,u2 = zr2,family = fam1[k,i],method = "mle", se = FALSE)$par
        # C_1|2(u1|u2)
        direct <- BiCopHfunc2(zr1, zr2, family = fam1[k,i], par = MLest, par2 = 0, check.pars = FALSE)
        # C_2|1(u2|u1)
        indirect <- BiCopHfunc2(zr2, zr1, family = fam1[k,i], par = MLest, par2 = 0, check.pars = FALSE)
        Vdir[[i]]=direct
        Vindir[[i]]=indirect
        Params[k,i]=MLest
        DepMeasure[k,i]=cor(Z[,1],Z[,2], method = "kendall")
        #DepMeasure[k,i]=VineCopula:::TauMatrix(Z)[1,2] # same as above
        EffectSamp[k,i]=nrow(Z)
        logLiks[k,i]=BiCopEst(u1 = zr1,u2 = zr2,family = fam1[k,i],method = "mle", se = FALSE)$logLik
      }  
    }
    if(k > 2){
      for(i in d:(k+1)){ #s=3 / i=5,4
        Out=PseudoCop(j = k-1, cind = i, data = data, StrMtx = M, fam = fam1, MaxMtr = MaxMat, Par = Params)
        ##  'PseudoCop' defines the pseudo data with the conditioned support
        ##  j: # of conditioning variables ,cind: fix the column
        m <- MaxMat[k, i]
        ##  For T_k, k>2, we use the pseudo data from 'Out'
        zr1 <- Out$Vdir[[i]] 
        zr2 <- if (m == M[k, i]) {
          Out$Vdir[[m]]
        }else {
          Out$Vindir[[m]] 
        }
        Z=cbind(zr1,zr2)
        MLest=BiCopEst(u1 = zr1,u2 = zr2,family = fam1[k,i],method = "mle",se = FALSE)$par
        Params[k,i]=MLest
        DepMeasure[k,i]=cor(Z[,1],Z[,2], method = "kendall")
        EffectSamp[k,i]=nrow(Z)
        logLiks[k,i]=BiCopEst(u1 = zr1,u2 = zr2,family = fam1[k,i],method = "mle", se = FALSE)$logLik
      }
    }
    }
  return(list("Params"=Params,"DepMeasure"=DepMeasure,"EffectSamp"=EffectSamp,"logLiks"=logLiks))
}
 