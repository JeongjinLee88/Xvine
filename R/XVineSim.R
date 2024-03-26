#' Simulation from X-Vine models
#'
#' @description
#' `XVineSim()` generates multivariate samples from X-Vine models 
#' given one of the variables is less than 1.
#'  The simulation algorithm selects the X-vine specification 
#'  associated with the conditioned variable and uses
#'  recursive approaches to return an \eqn{N \times d} data matrix.
#' For more details on simulation algorithms, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023)
#' @param N Integer; the number of samples to simulate.
#' @param XVS A list of three matrix components for each conditioning variable:
#'  * reproduced structure matrices
#'  * family matrices
#'  * parameter matrices.
#'  To specify the argument `XVS`, see [XVineSpec()].
#' @param k Integer; the \eqn{k}th, \eqn{k\in {1,\ldots,d}}, conditioning index to condition on.
#'
#' @return An \eqn{N\times d} data matrix used to calculate extremal functions in [ParetoSim()].
#' @export
#'
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#' @examples
#' ##  A 5-dim X-vine model
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#'                   0, 2, 1, 3, 2,
#'                   0, 0, 3, 1, 3,
#'                   0, 0, 0, 4, 1,
#'                   0, 0, 0, 0, 5),5,byrow = TRUE)
#' ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                    0, 0, 2, 2.5, 0.7,
#'                    0, 0, 0, 0.4, -0.3,
#'                    0, 0, 0, 0, 0.1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#' FamMtx <- matrix(c(0, 1, 2, 3, 4,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 3, 1,
#'                    0, 0, 0, 0, 1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#' ##  X-Vine specification
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' ##  Simulate from X-vine model when X_1 is less than 1.
#' XVineSim(N = 100, XVS = XVS, k = 1)
XVineSim <- function(N, XVS, k){
  
  d <- dim(XVS$xmat)[1]
  Rmat <- XVS$xmat[,,k]
  family <- XVS$fmat[,,k]
  par <- XVS$pmat[,,k]
  
  vdirect <- array(dim = c(d, d, N)) # stacked matrices
  vindirect <- array(dim = c(d, d, N))
  id <- 1:d + d * (0:(d - 1)) + rep((0:(N - 1)) * d^2, each = d)
  vdirect[id] = stats::runif(N * d)
  Diag=diag(Rmat)
  vdirect=vdirect[Diag,Diag,]
  vindirect[1, 1, ] <- vdirect[1, 1, ]
  
  CondMat=SecondArg(d = d,Rmat = Rmat) # above diagonal elements
  CondMat2=SecondArg2(d = d,Rmat = Rmat)
  if(d > 3){
    CondMat3=SecondArg3(d = d,Rmat = Rmat)
    CondMat4=SecondArg4(d = d,Rmat = Rmat)
  }
  
  for (i in 2:d) { #i=4;k=1
    for (k in (i-1):1) {
      diff=i-k
      ##  Choose the second arg
      if(Rmat[k,i]==0){
        u1 <- vdirect[which(Diag==i),which(Diag==i),]
      }else{
        if(k==1 & i==2){
          u1 <- vdirect[1,1,]
        }else if(k==1 & i > 2){
          u1 <- vdirect[1,which(Diag==Rmat[1,i]),]
        }else if(k > 1 & diff==1){ 
          u1 <- (if(CondMat[k,i])
            vdirect
            else vindirect)[k, i-1, ] 
        }else if(k > 1 & diff > 1){
          if(CondMat[k,i]){
            u1 <- vdirect[k, i-1, ]
          }else if(CondMat[k,i]==F & CondMat2[k,i]==T){ 
            u1 <- vindirect[k, i-1, ]
          }else if(CondMat[k,i]==F & CondMat2[k,i]==F){
            Loc=which(CondMat3[k,,i]==T)
            u1 <- (if(CondMat4[k,Loc,i])
              vdirect
              else vindirect)[k, Loc, ]
          }
        }
      }
      
      ##  Evaluate quantile functions and conditional dist
      if(k==1){
          vdirect[k,i,] <- InvTC(vdirect[k+1,i,], u1, par[k,i],family[k,i]) #lambda^{-1}_{2|1}(w2|u1=x1;gamma)
        if (i < d) {
            vindirect[k+1,i,] <- CondTC(u1, vdirect[k,i,], par[k,i],family[k,i]) #lambda_{1|2}(u1=x1|x2;gamma)
        } 
      }else if(k > 1){
        #vdirect[k,i,] <- InvCop(vdirect[k+1,i,], u1, par[k,i],family[k,i]) #C^{-1}_{2|1}(u2|u1;par)
        if(Rmat[k,i]==0){
          vdirect[k,i,] <- u1
        }else{
          vdirect[k,i,] <- BiCopHinv1(u1, vdirect[k+1,i,], family = family[k,i], par = par[k,i], par2 = 0, check.pars = FALSE)
        }
        if (i < d) {
          #vindirect[k+1,i,] <- CondCop(u1, vdirect[k,i,], par[k,i],family[k,i]) #C_{1|2}(u1|u2;par)
          if(Rmat[k,i]==0){
            vindirect[k+1,i,] <- u1
          }else{
            vindirect[k+1,i,] <- BiCopHfunc2(u1, vdirect[k,i,], family = family[k,i], par = par[k,i], par2 = 0, check.pars = FALSE)
          }
        }
      }
    }
  }
  vdirect_unordered=t(vdirect[1, , ])
  order.new=varIndexloc(Diag = Diag)
  return(vdirect_unordered[,order.new])
}
