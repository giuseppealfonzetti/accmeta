
<!-- README.md is generated from README.Rmd. Please edit that file -->

# accmeta

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)[![R-CMD-check](https://github.com/giuseppealfonzetti/accmeta/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/giuseppealfonzetti/accmeta/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Penalised estimation of meta-analysis models for diagnostic accuracy, in
terms of sensitivity, specificity and prevalence. Three estimators are
available: a trivariate linear mixed model (TLMM) using a within-study
normal approximation, the trivariate generalized linear mixed model
(TGLMM), and an iterative bootstrap (IB) that corrects the bias of the
TLMM by simulating from the TGLMM. IB fits come with a covariance matrix
and Wald or basic confidence intervals, built from a simulated Jacobian
of the TGLMM score.

## Installation

You can install the package via

``` r
pak::pkg_install("giuseppealfonzetti/accmeta")
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
#> [1,]  18   9   4  76
#> [2,]  39   6   2 121
#> [3,]  20   7   0  55
#> [4,]  30   1   3  19
#> [5,]  44   4   1  41
#> [6,]  33   7   7  77
```

Counts go in through `set_meta_data()`, which checks them and computes
the summaries the models need.

``` r
d <- set_meta_data(x)
#> empty cell in studies 3, 7, 10, 13: 'est' and 'wvar' are infinite. Rebuild with CC > 0 if needed.
d
#> <accmeta_data> 15 studies, continuity correction 0 
#>               SE        SP      PREV   n
#> min    0.7272727 0.7500000 0.0738255  46
#> median 0.9516129 0.9111111 0.4166667 112
#> max    1.0000000 0.9722222 0.6814159 168
```

A continuity correction is needed by the approximate model whenever a
cell is empty. You can use the CC argument:

``` r
d_cc <- set_meta_data(x, CC = 0.5)
d_cc
#> <accmeta_data> 15 studies, continuity correction 0.5 
#>               SE        SP      PREV   n
#> min    0.7083333 0.7456140 0.0794702  46
#> median 0.9444444 0.9021739 0.4193548 112
#> max    0.9893617 0.9594595 0.6782609 168
```

Now we can fit the models. Note that the TGLMM can be directly fitted on
the original data.

``` r
tlmm <- fit_tlmm(d_cc)
tglmm <- fit_tglmm(d)
ib <- fit_ib(d_cc)
ib
#> <accmeta_ib> 5 iterations, stopped on converged 
#>    SE    SP  PREV 
#> 0.939 0.906 0.372
```

The IB stops when the Mahalanobis distance of its last step falls below
`PRECISION / H` Monte Carlo variances (`"converged"`), when it stops
improving (`"stall"`), when the next update would leave the box
(`"boundary"`) or after `MAX_ITER` iterations (`"maxit"`). These
settings, and the box bounds, are set through `set_ib_control()`.

`coef()` extract estimated parameters:

``` r
coef(tlmm)
#> $MU
#>        eta         xi      gamma 
#>  2.4408070 -2.0612074 -0.4659787 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.5211082 0.2431792 0.4878588
#> xi    0.2431792 0.2943992 0.2372783
#> gamma 0.4878588 0.2372783 0.5972153
coef(tglmm)
#> $MU
#>        eta         xi      gamma 
#>  2.8737928 -2.2437993 -0.4812452 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.8546781 0.3150288 0.5893510
#> xi    0.3150288 0.3717160 0.1852082
#> gamma 0.5893510 0.1852082 0.6539028
coef(ib)
#> $MU
#>       eta        xi     gamma 
#>  2.741613 -2.269342 -0.522596 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.7715933 0.3834236 0.6289695
#> xi    0.3834236 0.4744769 0.2652382
#> gamma 0.6289695 0.2652382 0.6964915
```

All models are estimated, by default, with the recommended Wishart
prior. You can change it via

``` r
tglmm_mle <- fit_tglmm(d, PRIOR = set_prior(DEGREES = 4, SCALE = Inf))
coef(tglmm_mle)
#> $MU
#>        eta         xi      gamma 
#>  2.8031210 -2.2304939 -0.4804977 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.6191482 0.2841822 0.5459757
#> xi    0.2841822 0.3081583 0.1815984
#> gamma 0.5459757 0.1815984 0.5999210
min(eigen(coef(tglmm_mle)$SIGMA)$values)
#> [1] 0.03578618
min(eigen(coef(tglmm)$SIGMA)$values)
#> [1] 0.1286891
```

## Inference for the IB estimate

With `VCOV = TRUE`, `fit_ib()` also computes the Jacobian at the
estimate. `NCORES` spreads both the simulated fits of each iteration and
the Jacobian draws for parallel computations

``` r
ib <- fit_ib(d_cc, VCOV = TRUE, CONTROL = set_ib_control(NCORES = 4))
```

`update(ib, VCOV = TRUE)` adds the Jacobian to an existing fit, using
the fit’s `CONTROL`. `vcov()` and `confint()` work in the working
coordinates of the IB, i.e the means, the log standard deviations and
the Fisher-z partial correlations.

``` r
sqrt(diag(vcov(ib)))
#>       mu_eta        mu_xi     mu_gamma   log_sd_eta    log_sd_xi log_sd_gamma 
#>    0.3353733    0.2292392    0.2359522    0.3761369    0.3245307    0.2233054 
#>      z_rho12      z_rho13    z_rho23_1 
#>    0.5361700    0.6752042    0.6776142
confint(ib)
#>                   2.5 %      97.5 %
#> mu_eta        2.1489738  3.46361320
#> mu_xi        -2.7007043 -1.80210302
#> mu_gamma     -0.9528341 -0.02791858
#> log_sd_eta   -0.8532515  0.62117797
#> log_sd_xi    -1.0284955  0.24364161
#> log_sd_gamma -0.6118428  0.26349844
#> z_rho12      -0.4598593  1.64188842
#> z_rho13      -0.2099535  2.43679828
#> z_rho23_1    -1.3399962  1.31620282
confint(ib, METHOD = "basic")
#>                    2.5 %      97.5 %
#> mu_eta        2.11498171  3.42125809
#> mu_xi        -2.67502680 -1.79722682
#> mu_gamma     -0.94537195 -0.03217173
#> log_sd_eta   -0.76055283  0.68747094
#> log_sd_xi    -0.95459164  0.36555991
#> log_sd_gamma -0.55943702  0.30966978
#> z_rho12      -0.37215876  1.67155949
#> z_rho13      -0.08307779  2.41676793
#> z_rho23_1    -1.41323684  1.08607792
```

Intervals for the means map to the natural scale; specificity is one
minus the inverse logit of `mu_xi`, so its bounds swap.

``` r
ci <- confint(ib, parm = 1:3)
rbind(
  sensitivity = plogis(ci["mu_eta", ]),
  specificity = unname(1 - plogis(ci["mu_xi", 2:1])),
  prevalence = plogis(ci["mu_gamma", ])
)
#>                 2.5 %    97.5 %
#> sensitivity 0.8955728 0.9696345
#> specificity 0.8584047 0.9370682
#> prevalence  0.2783152 0.4930208
```
