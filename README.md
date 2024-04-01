
# X-vine models for multivariate extremes

The `Xvine` package is designed to analyze extremal dependence in
multivariate extremes using the graphical structure of regular vines, as
described in Kiriliouk et al. (2023). This modelling approach provides
flexibility and enables exploration of sparsity.

Computationally, the package employs recursive approaches based on
bivariate building blocks. The key components contain simulation from
X-vine models, sequential parameter estimation, selection of bivariate
parametric families, and vine structure selection.

## Installation

The [GitHub](https://github.com/JeongjinLee88/Xvine.git) version can be
downloaded via:

``` r
#install.packages("devtools")
devtools::install_github("JeongjinLee88/Xvine")
```

## Workflow

1.  X-vine specification

``` r
XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
```

2.  Simulation from Pareto distribution

``` r
Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
```

3.  Sequential parameter estimation via recursive approaches

``` r
XVineSeqEst(data = Dat_P, Rank = T, qt = 0.05, XVS=XVS, method = 'mle')
```

4.  Choosing bivariate parametric families along an X-vine sequence

``` r
# Specify the class of families
familyset_tc=c(1:4)
familyset_cop=c(0,1,3,4,5,6,13,14,16)
# Sequential model selection of bivariate parametric families
FamSelOut=XVineFamSel(data = Dat_P, Rank = TRUE, qt = 0.05, XVS = XVS
                      , famset_tc=familyset_tc, famset_cop=familyset_cop)
```

5.  Selection of X-vine tree sequence

## References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-kiriliouk2023x" class="csl-entry">

Kiriliouk, A., Lee, J., and Segers, J. (2023). X-vine models for
multivariate extremes. *arXiv Preprint arXiv:2312.15205*.

</div>

</div>
