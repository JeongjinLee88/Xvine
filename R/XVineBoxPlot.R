#' Box-plots of parameter estimates in terms of dependence measures
#'
#' @description
#' `XVineBoxplot` creates box-plots of parameter estimates for X-Vine models.
#'  We implement the maximum likelihood method for estimating parameters.
#'  The function provides ML estimates as well as dependence measure coefficients converted from the ML estimates.
#'  Given the regular vine structure, the function considers estimating parameters by either selecting the family of bivariate (tail) copula via the AIC
#'  or using the specified family.
#'  Four possible tail copula families are available:
#'  Husler-Reiss, Negative logistic, logistic, and Dirichlet models.
#'  For explicit formulas for the (lower) tail dependence coefficient \eqn{\chi_{L,e}=R(1,1)} for \eqn{e\in T_1} from the ML estimate, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#'  The function converts the ML estimate to the measure of Kendall's tau, \eqn{\tau_{e}}, for \eqn{e\in T_i}, \eqn{i=2,\ldots,d-1}.
#'  Any copula family with a single parameter can be considered. For more details, refer to Czado (2019, Section 3.5).
#' 
#' @param N Numeric; the sample size.
#' @param qt Numeric; a lower quantile as a threshold.
#' @param ite Numeric; the number of iterations for generating samples.
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#' @param RankT Logical; select either rank-based data (default) or data directly from the limiting distribution. 
#' @param familyset_tc Numeric vector; the model type of bivariate tail copulas.
#' Possible tail copula models are:
#' * 1=Husler-Reiss model
#' * 2=Negative logistic model
#' * 3=Logistic model
#' * 4=Dirichlet model
#' @param familyset_cop Numeric vector; the class of bivariate copula families with a single parameter
#' 
#' @return A list of parameter estimates matrices for repeated simulations and boxplots of the parameter estimates where the first d-1 plots are for the first level tree \eqn{T_1} and the next d-2 plots are for \eqn{T_2} and so on.
#' @export
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#'  Czado, C. (2019). Analyzing dependent data with vine copulas. Lecture Notes in Statistics, Springer, 222.
#'  
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
#' BoxOut2K05=XVineBoxplot(N = 2000, qt = 0.05, ite = ite, XVS = XVS, RankT = TRUE)
XVineBoxplot <- function(N, qt, ite, XVS, RankT=TRUE, familyset_tc=c(1:4), familyset_cop=c(0,1,3,4,5,6,13,14,16)){
  
  #  Matrices to be stored
  d <- dim(XVS$xmat)[1]
  MLE_WithoutFamSel <- array(NA,dim=c(d,d,ite)) # matrices of ML estimates without family selections
  MLE_WithFamSel <- array(NA,dim=c(d,d,ite)) # matrices of ML estimates with family selections
  MLE2Dep_WithoutFamSel <- array(NA,dim=c(d,d,ite)) # matrices of dep measures without family selections
  MLE2Dep_WithFamSel <- array(NA,dim=c(d,d,ite)) # matrices of dep measures with family selections
  FamMtxXVS=XVS$fmat[,,1]
  #DepMtx_Cop <- array(NA,dim=c(d,d,ite)) # matrices of parameter estimates
  #AvgEstMtx <- matrix(NA,d,d) # average matrix
  #AvgEstMtx=apply(EstMtx,1:2,mean) # sample average
  
  #  Repeat simulation 'ite' times
  if(RankT){
    for(i in 1:ite){
      tryCatch({
        Pout=ParetoSim(n = N, XVS=XVS)
        ##  ML estimates given bivariate parametric families
        OutWithoutFamSel=XVineSeqEst(data = Pout, Rank = RankT, qt = qt, XVS = XVS, method = 'mle')
        MLE_WithoutFamSel[,,i]=OutWithoutFamSel$Params
        MLE2Dep_WithoutFamSel[,,i]=ML2DepMtx(FamMtx = FamMtxXVS, ParMtx = OutWithoutFamSel$Params)
        ##  ML estimates after selecting bivariate parametric families
        OutWithFamSel=XVineFamSel(data = Pout, Rank = RankT, qt = qt, XVS = XVS, famset_tc = familyset_tc, famset_cop = familyset_cop,selectioncrit = 'AIC', effsampsize = 10, tau_threshold = 0.05)
        MLE_WithFamSel[,,i]=OutWithFamSel$Params
        MLE2Dep_WithFamSel[,,i]=ML2DepMtx(FamMtx = OutWithFamSel$famsel, ParMtx = OutWithFamSel$Params)
      }, error=function(e){})
    }
  }
  if(!isTRUE(RankT)){
    for(i in 1:ite){
      Pout=ParetoSim(n = N,XVS=XVS)
      ##  ML estimates given bivariate parametric families
      OutWithoutFamSel=XVineSeqEst(data = Pout, Rank = F, XVS = XVS, method = 'mle')
      MLE_WithoutFamSel[,,i]=OutWithoutFamSel$Params
      MLE2Dep_WithoutFamSel[,,i]=ML2DepMtx(FamMtx = FamMtxXVS, ParMtx = OutWithoutFamSel$Params)
      ##  ML estimates after selecting bivariate parametric families
      OutWithFamSel=XVineFamSel(data = Pout, Rank = F, XVS = XVS, famset_tc = familyset_tc, famset_cop = familyset_cop,selectioncrit = 'AIC', effsampsize = 10, tau_threshold = 0.05)
      MLE_WithFamSel[,,i]=OutWithFamSel$Params
      MLE2Dep_WithFamSel[,,i]=ML2DepMtx(FamMtx = OutWithFamSel$famsel, ParMtx = OutWithFamSel$Params)
    }
  }
  
  ##  Create boxplots
  ##  1. MLE via sequential parameter estimation (SPE)
  all.pairs <- combn(1:d, 2)
  MLest4eachpair=apply(all.pairs, 2, function(ind) rbind(MLE_WithoutFamSel[ind[1],ind[2],]))  
  data_longformat=melt(MLest4eachpair)
  data_longformat$Var2 <- as.factor(data_longformat$Var2)
  for(i in 1:(d-1)){
    a=cumsum((d-1):1)
    b=cumsum((d-1):1)+1
    if(i==1){
      data_longformat$Var1[data_longformat$Var2%in%(1:a[i])]=1  
    }else{
      data_longformat$Var1[data_longformat$Var2%in%(b[i-1]:a[i])]=i  
    }
  }
  data_longformat$Var1 <- as.character(data_longformat$Var1)
  TrueMtx=t(XVS$pmat[,,1]) # Specified MLEs
  y0s <- c(TrueMtx[lower.tri(TrueMtx)])
  dLines <- data.frame(X =1:choose(d,2) - 0.4,
                       Y = y0s,
                       Xend = 1:choose(d,2) + 0.4,
                       Yend = y0s,
                       Group = c("typeA"),
                       color = c("red"))
  bp_mle_SPE <- ggplot(data = data_longformat,aes(x=Var2,y=value,fill=Var1)) +
    geom_boxplot() +
    labs(title="",x="", y = "") +
    scale_fill_manual(labels=c(expression(T[1]),expression(T[2]),expression(T[3]),expression(T[4])),values=c("lightskyblue","lightgreen","lavender","papayawhip")) + 
    # legend labels
    theme(legend.title = element_blank(),
          legend.text = element_text(face="bold", size=15),
          legend.key.size = unit(1, 'cm'),
          legend.position="none",
          aspect.ratio = 1,
          axis.text.x = element_text(color="black", face="bold", size=9, angle = 20, hjust = 0.5),
          axis.text.y = element_text(face="bold", size=9),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey"),
          plot.margin = margin(0, 1, 0, 0, "cm")) +
    scale_y_continuous(breaks = seq(-1,4,1),limits = c(-1,4,1)) +
    geom_segment(data = dLines, color = "red",aes(x = X, y = Y, xend = Xend, yend=Yend), inherit.aes = FALSE) +
    scale_x_discrete(labels=c(expression(atop(chi[12],"\n(HR)")),expression(atop(chi[23],"\n(NL)")),expression(atop(chi[24],"\n(L)")),expression(atop(chi[45],"\n(Diri)")),expression(atop(tau[13~";"~2],"\n(Clay)")),expression(atop(tau[34~";"~2],"\n(Gum)")),expression(atop(tau[25~";"~4],"\n(Ga)")),expression(atop(tau[14~";"~23],"\n(Clay)")),expression(atop(tau[35~";"~24],"\n(Ga)")),expression(atop(tau[15~";"~234],"\n(Ga)"))))
  
  ## 2. MLE after selecting bivariate parametric families
  all.pairs <- combn(1:d, 2)
  MLest4eachpair=apply(all.pairs, 2, function(ind) rbind(MLE_WithFamSel[ind[1],ind[2],]))
  data_longformat=melt(MLest4eachpair)
  data_longformat$Var2 <- as.factor(data_longformat$Var2)
  for(i in 1:(d-1)){
    a=cumsum((d-1):1)
    b=cumsum((d-1):1)+1
    if(i==1){
      data_longformat$Var1[data_longformat$Var2%in%(1:a[i])]=1  
    }else{
      data_longformat$Var1[data_longformat$Var2%in%(b[i-1]:a[i])]=i  
    }
  }
  data_longformat$Var1 <- as.character(data_longformat$Var1)
  TrueMtx=t(XVS$pmat[,,1]) # Specified MLEs
  y0s <- c(TrueMtx[lower.tri(TrueMtx)])
  dLines <- data.frame(X =1:choose(d,2) - 0.4,
                       Y = y0s,
                       Xend = 1:choose(d,2) + 0.4,
                       Yend = y0s,
                       Group = c("typeA"),
                       color = c("red"))
  bp_mle_Fam <- ggplot(data = data_longformat,aes(x=Var2,y=value,fill=Var1)) +
    geom_boxplot() +
    labs(title="",x="", y = "") +
    scale_fill_manual(labels=c(expression(T[1]),expression(T[2]),expression(T[3]),expression(T[4])),values=c("lightskyblue","lightgreen","lavender","papayawhip")) + 
    # legend labels
    theme(legend.title = element_blank(),
          legend.text = element_text(face="bold", size=15),
          legend.key.size = unit(1, 'cm'),
          legend.position="none",
          aspect.ratio = 1,
          axis.text.x = element_text(color="black", face="bold", size=9, angle = 20, hjust = 0.5 ),
          axis.text.y = element_text(face="bold", size=9),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey"),
          plot.margin = margin(0, 1, 0, 0, "cm")) +
    scale_y_continuous(breaks = seq(-4,8.5,1.5),limits = c(-4,8.5,1.5)) +
    geom_segment(data = dLines, color = "red",aes(x = X, y = Y, xend = Xend, yend=Yend), inherit.aes = FALSE) +
    scale_x_discrete(labels=c(expression(theta[12]),expression(theta[23]),expression(theta[24]),expression(theta[45]),
                              expression(theta[13~";"~2]),expression(theta[34~";"~2]),expression(theta[25~";"~4]),expression(theta[14~";"~23]),
                              expression(theta[35~";"~24]),expression(theta[15~";"~234])))
  
  
  ##  3. Dependence measures through SPE
  all.pairs <- combn(1:d, 2)
  MLest4eachpair=apply(all.pairs, 2, function(ind) rbind(MLE2Dep_WithoutFamSel[ind[1],ind[2],]))
  data_longformat=melt(MLest4eachpair)
  data_longformat$Var2 <- as.factor(data_longformat$Var2)
  for(i in 1:(d-1)){
    a=cumsum((d-1):1)
    b=cumsum((d-1):1)+1
    if(i==1){
      data_longformat$Var1[data_longformat$Var2%in%(1:a[i])]=1  
    }else{
      data_longformat$Var1[data_longformat$Var2%in%(b[i-1]:a[i])]=i  
    }
  }
  data_longformat$Var1 <- as.character(data_longformat$Var1)
  TrueMtx_dep=t(DependMtx(XVS)) # Dependence measures
  y0s <- c(TrueMtx_dep[lower.tri(TrueMtx_dep)])
  dLines <- data.frame(X =1:choose(d,2) - 0.4,
                       Y = y0s,
                       Xend = 1:choose(d,2) + 0.4,
                       Yend = y0s,
                       Group = c("typeA"),
                       color = c("red"))
  bp_dep_SPE <- ggplot(data = data_longformat,aes(x=Var2,y=value,fill=Var1)) +
    geom_boxplot() +
    labs(title="",x="", y = "") +
    scale_fill_manual(labels=c(expression(T[1]),expression(T[2]),expression(T[3]),expression(T[4])),values=c("lightskyblue","lightgreen","lavender","papayawhip")) + 
    # legend labels
    theme(legend.title = element_blank(),
          legend.text = element_text(face="bold",size=15),
          legend.key.size = unit(1, 'cm'),
          legend.position="none",
          aspect.ratio = 1,
          axis.text.x = element_text(color="black", face="bold", size=9, angle = 20, hjust = 0.5 ),
          axis.text.y = element_text(face="bold", size=9),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey"),
          plot.margin = margin(0, 1, 0, 0, "cm")) +
    scale_x_discrete(labels=c(expression(atop(chi[12],"\n(HR)")),expression(atop(chi[23],"\n(NL)")),expression(atop(chi[24],"\n(L)")),expression(atop(chi[45],"\n(Diri)")),expression(atop(tau[13~";"~2],"\n(Clay)")),expression(atop(tau[34~";"~2],"\n(Gum)")),expression(atop(tau[25~";"~4],"\n(Ga)")),expression(atop(tau[14~";"~23],"\n(Clay)")),expression(atop(tau[35~";"~24],"\n(Ga)")),expression(atop(tau[15~";"~234],"\n(Ga)")))) +
    scale_y_continuous(breaks = seq(-0.4,0.8,0.2),limits = c(-0.4,0.8)) +
    geom_segment(data = dLines, color = "red",aes(x = X, y = Y, xend = Xend, yend=Yend), inherit.aes = FALSE)
  
  ##  4. Dependence measures after selecting bivariate parameteric families
  all.pairs <- combn(1:d, 2)
  MLest4eachpair=apply(all.pairs, 2, function(ind) rbind(MLE2Dep_WithFamSel[ind[1],ind[2],]))
  data_longformat=melt(MLest4eachpair)
  data_longformat$Var2 <- as.factor(data_longformat$Var2)
  for(i in 1:(d-1)){
    a=cumsum((d-1):1)
    b=cumsum((d-1):1)+1
    if(i==1){
      data_longformat$Var1[data_longformat$Var2%in%(1:a[i])]=1  
    }else{
      data_longformat$Var1[data_longformat$Var2%in%(b[i-1]:a[i])]=i  
    }
  }
  data_longformat$Var1 <- as.character(data_longformat$Var1)
  TrueMtx_dep=t(DependMtx(XVS)) # Dependence measures
  y0s <- c(TrueMtx_dep[lower.tri(TrueMtx_dep)])
  dLines <- data.frame(X =1:choose(d,2) - 0.4,
                       Y = y0s,
                       Xend = 1:choose(d,2) + 0.4,
                       Yend = y0s,
                       Group = c("typeA"),
                       color = c("red"))
  bp_dep_Fam <- ggplot(data = data_longformat,aes(x=Var2,y=value,fill=Var1)) +
    geom_boxplot() +
    labs(title="",x="", y = "") +
    scale_fill_manual(labels=c(expression(T[1]),expression(T[2]),expression(T[3]),expression(T[4])),values=c("lightskyblue","lightgreen","lavender","papayawhip")) + 
    # legend labels
    theme(legend.title = element_blank(),
          legend.text = element_text(face="bold",size=15),
          legend.key.size = unit(1, 'cm'),
          legend.position="none",
          aspect.ratio = 1,
          axis.text.x = element_text(color="black", face="bold", size=9, angle = 20, hjust = 0.5 ),
          axis.text.y = element_text(face="bold", size=9),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "solid"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey"),
          plot.margin = margin(0, 1, 0, 0, "cm")) +
    scale_x_discrete(labels=c(expression(chi[12]),expression(chi[23]),expression(chi[24]),expression(chi[45]),
                              expression(tau[13~";"~2]),expression(tau[34~";"~2]),expression(tau[25~";"~4]),expression(tau[14~";"~23]),
                              expression(tau[35~";"~24]),expression(tau[15~";"~234]))) +
    scale_y_continuous(breaks = seq(-0.4,0.8,0.2),limits = c(-0.4,0.8)) +
    geom_segment(data = dLines, color = "red",aes(x = X, y = Y, xend = Xend, yend=Yend), inherit.aes = FALSE)
  
  
  return(list("MLE_WithoutFamSel"=MLE_WithoutFamSel,"MLE_WithFamSel"=MLE_WithFamSel,"MLE2Dep_WithoutFamSel"=MLE2Dep_WithoutFamSel,"MLE2Dep_WithFamSel"=MLE2Dep_WithFamSel,'bp_mle_SPE'=bp_mle_SPE, 'bp_mle_Fam'=bp_mle_Fam, 'bp_dep_SPE'=bp_dep_SPE, 'bp_dep_Fam'=bp_dep_Fam))
}
