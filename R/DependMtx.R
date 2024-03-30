#' Dependence measure matrix for X-vine models
#' 
#' @description
#' `DependMtx` provides a summary matrix of bivariate dependence measures for X-vine models, using parametric relationships between
#' parameters associated with (tail) copula models and (tail) dependence measures.
#'  The first row includes the (lower) tail dependence coefficient \eqn{\chi_{L,e}=R(1,1)} for \eqn{e\in T_1}.
#'  Possible tail copula families include:
#'  Husler-Reiss, Negative logistic, logistic, and Dirichlet models. 
#'  For explicit formulas, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#'   and the remaining rows include the measure of Kendall's tau, \eqn{\tau_{e}}, for \eqn{e\in T_i}, \eqn{i=2,\ldots,d-1}.
#'  Any copula family with a single parameter can be considered. For more details, refer to Czado (2019, Section 3.5).
#' 
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#'
#' @return A \eqn{d\times d} upper triangular matrix of dependence measures 
#' using the parameter relationships.
#' 
#' @export
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#'  Czado, C. (2019). Analyzing dependent data with vine copulas. Lecture Notes in Statistics, Springer, 222.
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
#' # X-vine speicification
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' # Dependence matrix for X-vine models
#' DependMtx(XVS)
DependMtx<- function(XVS)
{
  d=dim(XVS$xmat[,,1])[1]
  par <- XVS$pmat[,,1]
  fam1 <- XVS$fmat[,,1]
  DepMtx=matrix(0,d,d)
  
  for(i in 1:(d-1)){
    if(i==1){
      for(j in (i+1):d){
        if(fam1[1,j]==1){ # HR
          DepMtx[i,j]=2-2*pnorm(sqrt(par[1,j])/2)
        }
        if(fam1[1,j]==2){ # NL
          DepMtx[i,j]=2^(-1/par[1,j])
        }
        if(fam1[1,j]==3){ # L
          DepMtx[i,j]=2-2^(1/par[1,j])
        }
        if(fam1[1,j]==4){ # Diri
          DepMtx[i,j]=integrate(Vectorize(function(x)pbeta(1/(1+x),par[1,j]+1,par[1,j])),lower = 0,upper = 1)$value
        }
      }
    }
    if(i > 1){
      for(j in (i+1):d){
        if(fam1[i,j]==1){ # Gauss
          DepMtx[i,j]=2/pi*asin(par[i,j])
        }
        if(fam1[i,j]==3){ # Clayton
          DepMtx[i,j]=par[i,j]/(par[i,j]+2)
        }
        if(fam1[i,j]==4){ # Gumbel
          DepMtx[i,j]=1-1/par[i,j]
        }else{
          DepMtx[i,j]=VineCopula::BiCopPar2Tau(family = fam1[i,j],par = par[i,j])
        }
      }
    }
  }
  return(DepMtx)
}
