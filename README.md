
<!-- README.md is generated from README.Rmd. Please edit that file -->
ConservativeEstimates
=====================

ConservativeEstimates is a R package to efficiently compute conservative estimates of excursion sets of functions under Gaussian random field priors. It implements an efficient estimator for orthant probabilites of high-dimensional Gaussian vectors. See the paper [Azzimonti, D. and Ginsbourger D. (2016)](https://hal.archives-ouvertes.fr/hal-01289126) for more details.

### Features

The package main functions are:

-   `conservativeEstimate` : the main function for conservative estimates computation. Requires the mean and covariance of the posterior field at a discretization design.

-   `ProbaMax`: the main function for high dimensional othant probabilities. Computes *P(max X \> t)*, where *X* is a Gaussian vector and *t* is the selected threshold. The function computes the probability with the decomposition explained [here](https://hal.archives-ouvertes.fr/hal-01289126). It implements both the `GMC` and `GANMC` algorithms.

### Installation

To install the package run the following code from a R console:

``` r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("dazzimonti/ConservativeEstimates")
```

### References

Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at [hal-01289126](https://hal.archives-ouvertes.fr/hal-01289126)
