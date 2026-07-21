
<!-- README.md is generated from README.Rmd. Please edit that file -->

# accmeta

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Meta-analysis of diagnostic accuracy, in terms of sensitivity,
specificity and prevalence. Three estimators are available: a trivariate
linear mixed model (TLMM) using a within-study normal approximation, the
trivariate generalized linear mixed model (TGLMM), and an iterative
bootstrap that corrects the bias of the TLMM by simulating from the
TGLMM.

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

Raw estimates:

``` r

# raw estimates
res <- rbind(
  truth = th, 
  TLMM = tlmm$THETA, 
  TGLMM = tglmm$THETA, 
  IB = ib$THETA)
res
#>           [,1]      [,2]       [,3]        [,4]      [,5]       [,6]      [,7]
#> truth 2.940000 -2.200000 -0.4000000  0.09530000 0.4000000 -0.5108000 0.3000000
#> TLMM  2.594077 -1.889481 -0.4257916 -0.21870787 0.5089039 -0.9131450 0.6604556
#> TGLMM 2.934820 -2.030976 -0.4404185 -0.06004097 0.5000930 -0.7244263 0.6832826
#> IB    2.919326 -2.105983 -0.4763200 -0.06417762 0.6619840 -0.8391709 0.7316053
#>              [,8]       [,9]
#> truth  0.20000000 -0.6931000
#> TLMM   0.02162351 -1.1240082
#> TGLMM  0.02122976 -1.0660574
#> IB    -0.19143915 -0.9512751
```
