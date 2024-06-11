#' Simulation for the selection of bivariate parametric families
#'
#' @description
#' `XVineFamSelSimul` returns simulation results of selecting bivariate parametric families
#' by calling the function [XVineFamSel()].
#' All information about the simulation result is saved in the list of `FamSelOut`. 
#'  
#' @param N Numeric; the sample size.
#' @param Rank Logical; whether rank transformation is performed or not (\code{Rank=T}; default).
#' @param qt Numeric; a lower enough threshold used in the rank transformation (e.g. switches from Pareto scale to uniform scale)
#' @param ite Numeric; the number of iterations for simulation.
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#' @param familyset_tc Numeric vector; the model type of bivariate tail copulas.
#' Possible tail copula models are:
#' * 1=Husler-Reiss model
#' * 2=Negative logistic model
#' * 3=Logistic model
#' * 4=Dirichlet model
#' @param familyset_cop Numeric vector; the class of bivariate copula families with a single parameter.
#' @param selectioncrit Character; the information criteria for model selections (default: `AIC`).
#' @param effsampsize Integer; indicates the effective sample size (Default: 10) for the independence copula.
#' @param tau_threshold Numeric; indicates the value of the Kendall's tau (Default: 0.05) for the independence copula.
#' @param edgeIndex Numeric vector; the index set that indicates the location of the matrix for the corresponding edge.
#'
#' @return 
#' 1. `FamSelOut`: A list of \eqn{d\times d} strict upper triangular matrices from the function [XVineFamSel()] for iterations.
#' 2. `EffSizePercent`: A \eqn{d\times d} strict upper triangular matrix including the percentages of the effective sample size.
#' 3. `GlobalProp`: Numeric; the proportion of correctly selected families over all trees.
#' 4. `TreeProp`: Numeric vector; the proportion of correctly selected families for each level tree.
#' 5. `EdgeProp`: A \eqn{d\times d} strict upper triangular matrix including the proportion of correctly selected families for each edge in each tree.
#' 6. `SelectIndCop`: Numeric; the number of times the algorithm chooses the independence copula for the specified edge by `edgeIndex`.
#' @export
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
#' # Approx 50 secs for N=4000
#' XVineFamSelSimul(N = 1000, qt = 0.05, ite = 100, XVS = XVS
#' ,edgeIndex = c(4,5), familyset_tc=c(1:4), familyset_cop = c(0,1,3,4,5,6,13,14,16))
XVineFamSelSimul <- function(N, Rank=TRUE, qt, ite, XVS, edgeIndex, familyset_tc, familyset_cop
                             , selectioncrit = 'AIC', effsampsize = 10, tau_threshold = 0.05)
{
  d <- dim(XVS$xmat)[1]
  EdgeTotal=sum(1:(d-1))
  FamSelOut <- list()
  
  ##  Repeat generating Pareto samples and performing family selections
  for(i in 1:ite){
    Dat_P=ParetoSim(n = N, XVS = XVS)
    FamSelOut[[i]] <- XVineFamSel(data = Dat_P, Rank = Rank, qt = qt, XVS = XVS, famset_tc = familyset_tc, famset_cop = familyset_cop
                                  , selectioncrit = selectioncrit,effsampsize = effsampsize, tau_threshold = tau_threshold)
  }
  
  ##  Extract the effective sample size matrices
  EffSizeComb=lapply(1:ite, function(i)FamSelOut[[i]]$EffectSamp)
  EffSizeComb=array(unlist(EffSizeComb),dim=c(d,d,ite))
  ##  Convert to the percentage of effective sample size with respect to the sample size N
  Avg.EffectSize=apply(X = EffSizeComb, c(1,2), function(x) mean(x,na.rm=T))
  EffSizePercent=(Avg.EffectSize/N)*100 # percentages of effective sample size
  #mean(EffSizeComb[3,4,] < 20) # percentage of effective sample size being less than a certain cutoff over iterations
  
  ####  Proportion of correctly selected families
  
  ##  1. Proportion of correctly selected families over all trees
  CountComb=lapply(1:ite,function(i)FamSelOut[[i]]$count)
  CountComb=array(unlist(CountComb),dim=c(dim(XVS$xmat)[1],dim(XVS$xmat)[1],ite))
  GlobalProp=sum(apply(X = CountComb,c(1,2),function(x)sum(x)),na.rm = T)/(EdgeTotal*ite)
  
  ##  2. The overall proportion of correctly selected families for each tree level
  Num=((d-1):1)*(ite)
  TreeProp=rowSums(apply(X = CountComb,c(1,2),function(x)sum(x)),na.rm=T)[1:(d-1)]/Num
  
  ##  3. Proportion of correctly selected families for each edge in each tree
  EdgeProp=apply(X = CountComb, c(1,2), function(x) mean(x,na.rm=T))  
  
  ##  4. How many times the algorithm selects the independence copula for a specific edge
  FamComb=lapply(1:ite,function(i)FamSelOut[[i]]$famsel)
  FamComb=array(unlist(FamComb),dim=c(d,d,ite))
  SelectIndCop=sum(FamComb[edgeIndex[1],edgeIndex[2],]==0)/ite 
  # if edge_index=(4,5), how many times the algorithm chooses the independence copula for the deepest edge.
  
  return(list("FamSelOut"=FamSelOut,"EffSizePercent"=EffSizePercent,"GlobalProp"=GlobalProp,"TreeProp"=TreeProp,"EdgeProp"=EdgeProp,"SelectIndCop"=SelectIndCop))
}

