% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XVineModelFit.R
\name{XVineModelFit}
\alias{XVineModelFit}
\title{Fit X-vine models}
\usage{
XVineModelFit(
  data,
  Rank = TRUE,
  Rank_chiL = FALSE,
  Rank_Vario = FALSE,
  Rank_chi3 = TRUE,
  Chi3_graph = FALSE,
  Vario_graph = FALSE,
  MST1_HR = FALSE,
  quan = 0.2,
  N = 2000,
  XVS,
  tcfamset = c(1, 2, 3, 4),
  pcfamset = c(0, 1, 3, 4, 5, 6, 13, 14, 16),
  selectioncrit = "AIC",
  BIC_graph = TRUE,
  trunclevel = FALSE,
  progress = TRUE,
  treecritT1 = "chi",
  treecritT2 = "tau",
  si = 0.9,
  effsampsize = 10,
  tau_threshold = 0.05,
  se = FALSE,
  weights = NA,
  cores = 1
)
}
\arguments{
\item{data}{A \eqn{n\times d} data matrix of multivariate Pareto samples.}

\item{Rank}{Logical; whether rank transformation is performed or not (\code{Rank=T}; default).}

\item{Rank_chiL}{Logical; whether rank-based samples are used or not for creating the chi-plot of the lower pairwise tail dependence measure.
If \code{Rank_chiL=FALSE}, then \code{ChiLPlot} compares model-based chi's with fitted chi's via Monte Carlo simulation.}

\item{Rank_Vario}{Logical; whether rank-based samples are used or not for the plot of the empirical pairwise variogram versus the fitted pairwise variogram.}

\item{Rank_chi3}{Logical; whether rank-based samples are used or not for creating the chi-plot of the trivariate tail dependence measure.
If \code{Rank_chi3=FALSE}, then \code{ChiLPlot} compares model-based chi's with fitted chi's via Monte Carlo simulation.}

\item{Chi3_graph}{Logical; whether the chi-plot of the trivariate tail dependence measure is plotted or not.}

\item{Vario_graph}{Logical; whether the plot of the empirical pairwise variogram versus the fitted pairwise variogram is created or not.}

\item{MST1_HR}{Logical; the minimum spanning tree for the Husler-Reiss model is plotted or not.}

\item{quan}{Numeric; a lower threshold for the rank transformation. It switches from Pareto scale to uniform scale.}

\item{N}{Numeric; the sample size to draw samples from the limiting distribution of multivariate Pareto distribution.}

\item{XVS}{A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:\code{\link[=XVineSpec]{XVineSpec()}}.}

\item{tcfamset}{Numeric vector; the class of bivariate exponent measures.}

\item{pcfamset}{Numeric vector; the class of bivariate pair-copula models.}

\item{selectioncrit}{Character string; indicates the selection criteria for the class of bivariate exponent measures or pair-copulas.}

\item{BIC_graph}{Logical; whether the BIC graph over the tree level is plotted or not.}

\item{trunclevel}{Numeric; indicates the specified truncation level (\code{trunclevel=NULL}; default).}

\item{progress}{Logical; whether the progress of selecting vine tree structures is printed.}

\item{treecritT1}{A character string, indicating the tree criterion for the first tree.}

\item{treecritT2}{A character string, indicating the tree criterion for subsequent trees.}

\item{si}{Numeric; a tuning parameter for mBIC (\eqn{\si_0=0.9}; default).}

\item{effsampsize}{Numeric; the specified effective sample size for the independence copula (\eqn{n_{D_e}}<10; default).}

\item{tau_threshold}{Numeric; the specified Kendall's tau value for the independence copula (\eqn{\hat{\tau}_e}<0.05; default)}

\item{se}{Logical; whether standard errors for ML estimators are reported.}

\item{weights}{Logical; whether weights should be assigned to observations when missing values exist.}

\item{cores}{Numeric; indicates the number of cores for parallel computing (optional).}
}
\value{
A nested list object containing:
\itemize{
\item MST: Maximum spanning trees
\item VineTree: Fitted X-vine models
\item mBIC_g: the plot of mBIC values across tree levels
\item BIC_g: the plot of BIC values across tree levels
\item XVS_spec: the X-vine specification from fitted X-vine models
\item emp_chimat: the matrix of empirical pairwise chi's
\item XVine_chimat: the matrix of fitted pairwise chi's
\item ChiLPlot: the chi-plot of the lower pairwise tail dependence measure
\item Chi3Plot: the chi-plot of the trivariate tail dependence measure
\item VarioPlot: the plot of the pairwise empirical varigoram versus fitted variogram
\item TruncLevelStar: the optimal truncation level determined as the minimum between the tree level with the lowest mBIC value
and the tree level such that all pair-copulas are set to independence copula
\item TruncLevel_mBICmin: the truncation level with the lowest mBIC value
}
}
\description{
Fits X-Vine models to a d-dimensional inverted-Pareto data set. The function selects vine tree structures, using maximum spanning tree.
and select bivariate exponent measure classes and pair-copula classes and estimate parameters for the corresponding classes
Users can specify the type of edge weights for the first tree and subsequent trees separately and the class of bivariate exponent measures and pair-copulas
Note that if you do not use samples from the limiting distribution, then the \code{Rank} transformation is required to create chi-plots or variogram-plots.
That is, if you fit the X-vine model to the real dataset, then the rank transformation is required throughout.
}
\examples{
StrMtx <- matrix(c(1, 1, 2, 2, 4,
0, 2, 1, 3, 2,
0, 0, 3, 1, 3,
0, 0, 0, 4, 1,
0, 0, 0, 0, 5),5,byrow = TRUE)
ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
                0, 0, 2, 2.5, 0.7,
                0, 0, 0, 0.4, -0.3,
                0, 0, 0, 0, 0.1,
                0, 0, 0, 0, 0),5,byrow = TRUE)
FamMtx <- matrix(c(0, 1, 2, 3, 4,
                   0, 0, 3, 4, 1,
                   0, 0, 0, 3, 1,
                   0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0),5,byrow = TRUE)
# X-vine specification
XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
XVineFitOut=XVineModelFit(data = 1/Dat_P, N = 2000, XVS = XVS
                         , Rank = FALSE, Rank_chiL = FALSE, Rank_chi3 = FALSE
                         , Chi3_graph = TRUE, Rank_Vario = FALSE, Vario_graph = FALSE
                         , tcfamset = c(1:4), pcfamset = c(0,1,3,4,5,6,13,14,16)
                         , selectioncrit = "AIC", trunclevel = FALSE, progress = TRUE
                         , treecritT1 = "chi", treecritT2 = "tau"
                         , effsampsize = 10, tau_threshold = 0.05
                         , weights = NA, cores=1)
}
