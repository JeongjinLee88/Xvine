#' Box-plots of ML estimates for a specific edge across different quantiles
#' 
#' @description
#' Given the vine tree sequence and the families of bivariate parametric models,
#' `XVineBoxPlotEachPar()` creates box-plots of maximum likelihood (ML) estimates for a specific edge across three different quantiles.
#' The functions takes as input a list of parameter estimate matrices from the function `XVineBoxPlot()`.
#' The purpose of this function is reproduce Figure 8 in Kiriliouk, A., Lee, J., & Segers, J. (2023).
#'
#' 
#' @param BoxOut1 A list of parameter estimate matrices for repeated simulations with a specific \eqn{(k,n)} from the output of `XVineBoxPlot()`
#' @param BoxOut2 A list of parameter estimate matrices for repeated simulations with a specific \eqn{(k,n)} from the output of `XVineBoxPlot()`
#' @param BoxOut3 A list of parameter estimate matrices for repeated simulations with a specific \eqn{(k,n)} from the output of `XVineBoxPlot()`
#' @param RowInd Numeric; A location of a parameter matrix in row
#' @param ColInd Numeric; A location of a parameter matrix in column
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#'
#' @return A `BoxOnePar` object returns the boxplot of ML estimates for a specific edge across three different quantiles that the user specifies.
#'  
#' @export
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.

#' @examples
#' # Specify the number of iterations
#' ite=200
#' # Specify the structure matrix
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#' 0, 2, 1, 3, 2,
#' 0, 0, 3, 1, 3,
#' 0, 0, 0, 4, 1,
#' 0, 0, 0, 0, 5),5,byrow = TRUE)
#' # Specify the parameter matrix
#' ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                 0, 0, 2, 2.5, 0.7,
#'                 0, 0, 0, 0.4, -0.3,
#'                 0, 0, 0, 0, 0.1,
#'                 0, 0, 0, 0, 0),5,byrow = TRUE)
#' # Specify the family matrix                 
#' FamMtx <- matrix(c(0, 1, 2, 3, 4,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 3, 1,
#'                    0, 0, 0, 0, 1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#' # X-vine specification
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' # n = 1000, k/n = 0.02, 0.05, 0.1
#' # Not to run (less than 1 min for a 5-dim X-vine model)
#' # BoxOut1K02=XVineBoxplot(N = 1000, qt = 0.02, ite = ite, XVS = XVS, RankT = T)
#' # BoxOut1K05=XVineBoxplot(N = 1000, qt = 0.05, ite = ite, XVS = XVS, RankT = T)
#' # BoxOut1K1=XVineBoxplot(N = 1000, qt = 0.1, ite = ite, XVS = XVS, RankT = T)
#' # BoxPar1K=XVineBoxPlotEachPar(BoxOut1 = BoxOut1K02$MLE_WithoutFamSel,
#' # BoxOut2 = BoxOut1K05$MLE_WithoutFamSel,
#' # BoxOut3 = BoxOut1K1$MLE_WithoutFamSel,
#' # RowInd = 1, ColInd = 2, XVS = XVS)
XVineBoxPlotEachPar <- function(BoxOut1,BoxOut2,BoxOut3,RowInd=1,ColInd=2,XVS)
{
  A=c(BoxOut1[RowInd,ColInd,],BoxOut2[RowInd,ColInd,],BoxOut3[RowInd,ColInd,])
  B=rep(1:3,each=dim(BoxOut1)[3])
  MTX=cbind(A,B)
  MTX=as.data.frame(MTX)
  MTX$B=as.factor(MTX$B)
  X<-Y<-Xend<-Yend<-NULL
  
  y0s=XVS$pmat[RowInd,ColInd,1] # Specified MLEs
  dLines <- data.frame(X =1:3 - 0.4,
                       Y = y0s,
                       Xend = 1:3 + 0.4,
                       Yend = y0s,
                       color = "red")
  BoxOnePar <- ggplot(data = MTX,aes(x=B,y=A)) +
    geom_boxplot(lwd=0.8) +
    geom_segment(data = dLines, color = "red",aes(x = X, y = Y, xend = Xend, yend=Yend), inherit.aes = FALSE) +
    labs(title="",x="n=1000", y = "") +
    theme(axis.text.x = element_text(color="black", face="bold", size=13),
          axis.text.y = element_text(face="bold", size=14),
          axis.title=element_text(size=15,face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey"),
          plot.margin = margin(0, 1, 0, 0, "cm")) +
    scale_y_continuous(breaks = seq(0.5,3,1),limits = c(0.5,3,1)) +
    scale_x_discrete(labels=c("k=20","k=50","k=100"))
  return(list("BoxOnePar"=BoxOnePar))
}

