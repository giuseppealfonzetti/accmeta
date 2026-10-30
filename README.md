
<!-- README.md is generated from README.Rmd. Please edit that file -->

# accmeta

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)[![R-CMD-check](https://github.com/giuseppealfonzetti/accmeta/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/giuseppealfonzetti/accmeta/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Penalised estimation of meta-analysis models for diagnostic accuracy, in
terms of sensitivity and specificity (bivariate model) or sensitivity,
specificity and prevalence (trivariate model), combining cohort and
case-control studies. Three estimators are available: a linear mixed
model (LMM) using a within-study normal approximation, the generalized
linear mixed model (GLMM), and an iterative bootstrap (IB) that corrects
the bias of the LMM by simulating from the GLMM. All are fitted through
`fit_meta_model()`, with the estimator’s options set by `set_control()`.
IB fits come with a covariance matrix and Wald or basic confidence
intervals, built from a simulated Jacobian of the TGLMM score.

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
tlmm <- fit_meta_model(d_cc, MODEL = "trivariate", ESTIMATOR = "lmm")
tglmm <- fit_meta_model(d, MODEL = "trivariate", ESTIMATOR = "glmm")
ib <- fit_meta_model(d_cc, MODEL = "trivariate", ESTIMATOR = "ib")
ib
#> <accmeta_ib> 6 iterations, stopped on converged 
#>    SE    SP  PREV 
#> 0.940 0.906 0.372
```

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
#> eta   0.8546781 0.3150288 0.5893511
#> xi    0.3150288 0.3717160 0.1852082
#> gamma 0.5893511 0.1852082 0.6539028
coef(ib)
#> $MU
#>        eta         xi      gamma 
#>  2.7466954 -2.2689605 -0.5218288 
#> 
#> $SIGMA
#>             eta        xi     gamma
#> eta   0.7834359 0.3877441 0.6318065
#> xi    0.3877441 0.4743261 0.2673085
#> gamma 0.6318065 0.2673085 0.6958735
```

All models are estimated, by default, with the recommended Wishart
prior, whose `DEGREES = NULL` means $p + 2$ for $p$ random effects. You
can change it via

``` r
tglmm_mle <- fit_meta_model(
  d,
  MODEL = "trivariate",
  ESTIMATOR = "glmm",
  PRIOR = set_prior(DEGREES = 4, SCALE = Inf)
)
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
#> [1] 0.03578619
min(eigen(coef(tglmm)$SIGMA)$values)
#> [1] 0.1286891
```

## Inference for the IB estimate

With `VCOV = TRUE` in `set_control()`, the IB also computes the Jacobian
at the estimate. `NCORES` spreads both the simulated fits of each
iteration and the Jacobian draws for parallel computations

``` r
ib <- fit_meta_model(
  d_cc,
  MODEL = "trivariate",
  ESTIMATOR = "ib",
  CONTROL = set_control("ib", VCOV = TRUE, NCORES = 4)
)
```

`update(ib, VCOV = TRUE)` adds the Jacobian to an existing fit, using
the fit’s `CONTROL`. `vcov()` and `confint()` work in the working
coordinates of the IB, i.e the means, the log standard deviations and
the Fisher-z partial correlations.

``` r
sqrt(diag(vcov(ib)))
#>       mu_eta        mu_xi     mu_gamma   log_sd_eta    log_sd_xi log_sd_gamma 
#>    0.3379689    0.2290867    0.2359493    0.3726743    0.3242787    0.2229079 
#>      z_rho12      z_rho13    z_rho23_1 
#>    0.5312978    0.6633973    0.6541761
confint(ib)
#>                   2.5 %      97.5 %
#> mu_eta        2.1510149  3.47582878
#> mu_xi        -2.7017213 -1.80371804
#> mu_gamma     -0.9524149 -0.02751072
#> log_sd_eta   -0.8388008  0.62205561
#> log_sd_xi    -1.0271236  0.24402566
#> log_sd_gamma -0.6109018  0.26288113
#> z_rho12      -0.4623547  1.62029435
#> z_rho13      -0.2093973  2.39107249
#> z_rho23_1    -1.2723168  1.29200628
confint(ib, METHOD = "basic")
#>                   2.5 %      97.5 %
#> mu_eta        2.1466417  3.42988049
#> mu_xi        -2.6803782 -1.80032506
#> mu_gamma     -0.9509505 -0.02475862
#> log_sd_eta   -0.7523939  0.66694223
#> log_sd_xi    -0.9457089  0.36499269
#> log_sd_gamma -0.5606429  0.31394127
#> z_rho12      -0.3835294  1.64579306
#> z_rho13      -0.1039431  2.35580691
#> z_rho23_1    -1.3324174  1.07129913
```

Intervals mapped onto the natural scale via

``` r
confint(ib, SCALE = "natural")
#>                    2.5 %      97.5 %
#> mu_eta         2.1510149  3.47582878
#> mu_xi         -2.7017213 -1.80371804
#> mu_gamma      -0.9524149 -0.02751072
#> sd_eta         0.4322286  1.86275320
#> sd_xi          0.3580353  1.27637709
#> sd_gamma       0.5428611  1.30067209
#> rho_eta_xi    -0.4320014  0.92466691
#> rho_eta_gamma -0.2063895  0.98338320
#> rho_xi_gamma  -0.3171718  0.84160083
```

Intervals mapped onto the reporting scale via

``` r
confint(ib, SCALE = "reporting")
#>                    2.5 %    97.5 %
#> sensitivity    0.8957636 0.9699921
#> specificity    0.8586009 0.9371281
#> prevalence     0.2783994 0.4931228
#> sd_eta         0.4322286 1.8627532
#> sd_xi          0.3580353 1.2763771
#> sd_gamma       0.5428611 1.3006721
#> rho_eta_xi    -0.4320014 0.9246669
#> rho_eta_gamma -0.2063895 0.9833832
#> rho_xi_gamma  -0.3171718 0.8416008
```
