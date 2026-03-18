# boot.pval

This R package provides functions for computing bootstrap p-values based
on `boot` objects, and convenience functions for bootstrap confidence
intervals and p-values for various regression models.

## Installation

To install the package from CRAN:

    install.packages("boot.pval")

To install the development version from Github:

    library(devtools)
    install_github("mthulin/boot.pval")

## Features

`boot.pval` can be used to:

- [Obtain p-values and confidence intervals for coefficients of a wide
  range of regression
  models.](https://mthulin.github.io/boot.pval/articles/boot_summary.html)
- [Perform one and two-sample bootstrap tests of location, with
  accompanying confidence
  intervals.](https://mthulin.github.io/boot.pval/articles/boot_t_test.html)
- [Create custom bootstrap
  tests.](https://mthulin.github.io/boot.pval/articles/custom_tests.html)

## Background

p-values can be computed by inverting the corresponding confidence
intervals, as described in Section 14.2 of [Thulin
(2024)](https://www.modernstatisticswithr.com/mathschap.html#confintequal)
and Section 3.12 in [Hall
(1992)](https://link.springer.com/book/10.1007/978-1-4612-4384-7). This
package contains functions for computing bootstrap p-values in this way.
The approach relies on the fact that:

- The p-value of the two-sided test for the parameter theta is the
  smallest alpha such that theta is not contained in the corresponding
  1-alpha confidence interval,
- For a test of the parameter theta with significance level alpha, the
  set of values of theta that aren’t rejected by the two-sided test
  (when used as the null hypothesis) is a 1-alpha confidence interval
  for theta.
