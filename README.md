
<!-- README.md is generated from README.Rmd. Please edit that file -->

# accmeta

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Penalised estimation of meta-analysis models for diagnostic accuracy, in
terms of sensitivity, specificity and prevalence. Three estimators are
available: a trivariate linear mixed model (TLMM) using a within-study
normal approximation, the trivariate generalized linear mixed model
(TGLMM), and an iterative bootstrap that corrects the bias of the TLMM
by simulating from the TGLMM.

## Installation

You can install the package via

``` r
devtools::install_github("giuseppealfonzetti/accmeta")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo giuseppealfonzetti/accmeta@HEAD
#> Rcpp    (1.1.1  -> 1.1.2 ) [CRAN]
#> TMB     (1.9.20 -> 1.9.22) [CRAN]
#> ucminf  (1.2.2  -> 1.2.3 ) [CRAN]
#> statmod (1.5.1  -> 1.5.2 ) [CRAN]
#> Installing 4 packages: Rcpp, TMB, ucminf, statmod
#> Installing packages into '/private/var/folders/4z/s_6y54qn43xdykvrmqg13kh40000gn/T/RtmpRdWZMH/temp_libpathd133129d66d3'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/4z/s_6y54qn43xdykvrmqg13kh40000gn/T//RtmpMsby4r/downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/4z/s_6y54qn43xdykvrmqg13kh40000gn/T/RtmpMsby4r/remotesd589358a5ff/giuseppealfonzetti-accmeta-4563fbf50ec7ec43d6619453d86ac66a9dda3343/DESCRIPTION’ ... OK
#> * preparing ‘accmeta’:
#> * checking DESCRIPTION meta-information ... OK
#> * cleaning src
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘accmeta_0.0.0.9000.tar.gz’
#> Installing package into '/private/var/folders/4z/s_6y54qn43xdykvrmqg13kh40000gn/T/RtmpRdWZMH/temp_libpathd133129d66d3'
#> (as 'lib' is unspecified)
```

## Example

Setup true parameter value

``` r
library(accmeta)

th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
theta2list(th)
#> $MU
#>   eta    xi gamma 
#>  2.94 -2.20 -0.40 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   1.2099754 0.4399955 0.3299966
#> xi    0.4399955 0.5200184 0.2400031
#> gamma 0.3299966 0.2400031 0.3800236
min(eigen(theta2list(th)$SIGMA)$values)
#> [1] 0.2000188
```

Simulate data from true model

``` r
set.seed(1)
x <- sim_data(15, th, sample(40:200, 15, TRUE))
head(x)
#>      n11 n10 n01 n00
#> [1,]  19  11   3  74
#> [2,]  32   9   3 124
#> [3,]  27   2   1  52
#> [4,]  22   1   2  28
#> [5,]  51   4   1  34
#> [6,]  37   6   8  73
```

Counts go in through `set_meta_data()`, which checks them and computes
the summaries the models need.

``` r
d <- set_meta_data(x)
#> empty cell in study 10: 'est' and 'wvar' are infinite. Rebuild with CC > 0 if needed.
d
#> <accmeta_data> 15 studies, continuity correction 0 
#>               SE        SP       PREV   n
#> min    0.7500000 0.6481481 0.08053691  46
#> median 0.9642857 0.8947368 0.45283019 112
#> max    1.0000000 0.9655172 0.68141593 168
```

A continuity correction is needed by the approximate model whenever a
cell is empty. You can use the CC argument:

``` r
d_cc <- set_meta_data(x, CC = 0.5)
d_cc
#> <accmeta_data> 15 studies, continuity correction 0.5 
#>               SE        SP       PREV   n
#> min    0.7307692 0.6454545 0.08609272  46
#> median 0.9482759 0.8846154 0.45454545 112
#> max    0.9807692 0.9545455 0.67826087 168
```

Now we can fit the models. Note that the TGLMM can be directly fitted on
the original data.

``` r
tlmm <- fit_tlmm(d_cc)
tglmm <- fit_tglmm(d)
ib <- fit_ib(d_cc)
ib
#> <accmeta_ib> 5 iterations, stopped on tol 
#>    SE    SP  PREV 
#> 0.949 0.891 0.383
```

`coef()` extract estimated parameters:

``` r
coef(tlmm)
#> $MU
#>        eta         xi      gamma 
#>  2.5940772 -1.8894812 -0.4257916 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.6457029 0.4089330 0.5307133
#> xi    0.4089330 0.4199930 0.3447851
#> gamma 0.5307133 0.3447851 0.5422777
coef(tglmm)
#> $MU
#>        eta         xi      gamma 
#>  2.9348201 -2.0309756 -0.4404185 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.8868478 0.4709505 0.6434649
#> xi    0.4709505 0.4849326 0.3519928
#> gamma 0.6434649 0.3519928 0.5859120
coef(ib)
#> $MU
#>       eta        xi     gamma 
#>  2.919326 -2.105983 -0.476320 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.8795409 0.6208340 0.6861275
#> xi    0.6208340 0.6249061 0.4015962
#> gamma 0.6861275 0.4015962 0.7210829
```

All models are estimated, by default, with the recommended Wishart
prior. You can change it via

``` r
tglmm_mle <- fit_tglmm(d, PRIOR = set_prior(DEGREES = 4, SCALE = Inf))
coef(tglmm_mle)
#> $MU
#>        eta         xi      gamma 
#>  2.8757287 -2.0176921 -0.4397477 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.6882836 0.4539968 0.5911395
#> xi    0.4539968 0.4211910 0.3304126
#> gamma 0.5911395 0.3304126 0.5367958
min(eigen(coef(tglmm_mle)$SIGMA)$values)
#> [1] 1.289971e-10
min(eigen(coef(tglmm)$SIGMA)$values)
#> [1] 0.07366225
```
