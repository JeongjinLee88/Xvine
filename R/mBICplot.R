#' Sensitive analysis for the mBIC-optimal truncation level across low quantiles.
#'
#' @description
#' `mBICplot()` plots the mBIC-optimal truncation levels across low quantiles
#' , which can be used for sensitive analysis in data application.
#' 
#' @param data A \eqn{n\times d} data matrix that is passed into the function `XVineModelFit_lightversion()`
#' performs the rank transformation.
#' @param Rank Logical; whether rank transformation is performed or not (\code{Rank=T}; default).
#' @param quan A numeric vector; a sequence of lower quantiles in ascending order.
#' @param famset_tc Numeric vector; the class of bivariate exponent measures.
#' @param famset_pc Numeric vector; the class of bivariate pair-copula models.
#' @param selectioncrit Character string; indicates the selection criteria for the class of bivariate exponent measures or pair-copulas.
#' @param trunclevel Numeric; indicates the specified truncation level (\code{trunclevel=NULL}; default).
#' @param treecritT1 A character string, indicating the tree criterion for the first tree.
#' @param treecritT2 A character string, indicating the tree criterion for subsequent trees.
#' @param effsampsize Numeric; the specified effective sample size for the independence copula (\eqn{n_{D_e}}<10; default).
#' @param tau_threshold Numeric; the specified Kendall's tau value for the independence copula (\eqn{\hat{\tau}_e}<0.05; default).
#' @param si Numeric; a tuning parameter for mBIC (\eqn{\si_0=0.9}; default).
#'
#' @return A `mBICvsQuan` object
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
#' Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
#' quan=seq(0.05,0.15,by=0.01)
#' mBICplot(Dat_P, Rank=TRUE, quan
#'         , famset_tc=c(1:4), famset_pc=c(0,1,3,4,5,6,13,14,16)
#'         , selectioncrit="AIC",trunclevel=FALSE
#'         , treecritT1="chi",treecritT2="tau"
#'         , effsampsize=10, tau_threshold=0.05, si=0.9)
mBICplot <- function(data, Rank=TRUE, quan, famset_tc=c(1:4), famset_pc=c(0,1,3,4,5,6,13,14,16)
                     , selectioncrit="AIC",trunclevel=FALSE,treecritT1="chi",treecritT2="tau"
                     , effsampsize=10, tau_threshold=0.05, si=0.9)
{
  Quantile<-mBIC<-NULL
  k=length(quan)
  mBICout=rep(NA,k)
  for(i in 1:k){
    Out=XVineModelFit_lightversion(data = data, Rank = Rank, quan=quan[i] ,BIC_graph = TRUE
                      , tcfamset = famset_tc, pcfamset = famset_pc, MST1_HR=FALSE, si=si
                      , selectioncrit = selectioncrit, trunclevel = trunclevel, progress = TRUE
                      , treecritT1 = treecritT1, treecritT2 = treecritT2
                      , effsampsize = effsampsize, tau_threshold = tau_threshold
                      , weights = NA, cores=1)
    mBICout[i]=Out$TruncLevel_mBICmin
  }
  #mBICout=rev(x = mBICout)
  mBICmtx=cbind(quan,mBICout)
  colnames(mBICmtx)=c("Quantile","mBIC")
  mBICmtx=as.data.frame(mBICmtx)
  mBICvsQuan <- ggplot(data = mBICmtx,aes(x=Quantile,y=mBIC)) +
                labs(y="mBIC values") +
                geom_point(size=2.5) +
                theme(legend.title = element_blank(), axis.title=element_text(size=14)
                 ,axis.title.y = element_text(face="bold")
                 ,axis.text.x = element_text(color="black", face="bold", size=14)
                 ,axis.text.y = element_text(face="bold", size=14)
                 ,panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "solid")
                 ,panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey")
                 ,panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',colour = "grey")
                 ,plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),aspect.ratio = 1)
          
  mBICvsQuan
}
