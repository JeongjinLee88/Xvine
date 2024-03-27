# X-vine models for multivariate extremes


The `Xvine` package is developed to analyze extremal dependence in multivariate extremes using the graphical structure of regular vines.
This modelling approach provides flexibility and allows for exploring sparsity. 

Computations proceed via recursive approaches in terms of bivariate building blocks. Some main components in the package contains: simulation from X-vine models, sequential parameter estimation, bivariate parametric family selection, and vine structure selection.

##  Installation

The [GitHub](https://github.com/JeongjinLee88/Xvine.git) version can be downloaded via:

```{r}
#install.packages("devtools")
devtools::install_github("JeongjinLee88/Xvine")
```

##  Workflow

1. X-vine specification

```{r}
XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
```

2. Simulation from Pareto distribution

```{r}
Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
```

3. Sequential parameter estimation via recursive approaches

```{r}
XVineSeqEst(data = Dat_P, Rank = T, qt = 0.05, XVS=XVS, method = 'mle')
```

4. Parametric family selection

```{r}

```

5. Vine structure selection

```{r}

```

