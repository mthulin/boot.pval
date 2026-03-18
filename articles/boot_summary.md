# Bootstrap p-values and confidence intervals for regression models

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

## Summaries for regression models

Summary tables with confidence intervals and p-values for the
coefficients of regression models can be obtained using the
`boot_summary` (most models) and `censboot_summary` (models with
censored response variables) functions. Currently, the following models
are supported:

- Linear models fitted using `lm`,
- Generalised linear models fitted using `glm` or
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html),
- Nonlinear models fitted using `nls`,
- Robust linear models fitted using
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html),
- Ordered logistic or probit regression models fitted (without weights)
  using [`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html),
- Linear mixed models fitted using
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) or
  `lmerTest::lmer`,
- Generalised linear mixed models fitted using
  [`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html).
- Cox PH regression models fitted using
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  (using `censboot_summary`).
- Accelerated failure time models fitted using
  [`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  or [`rms::psm`](https://rdrr.io/pkg/rms/man/psm.html) (using
  `censboot_summary`).
- Any regression model such that: `residuals(object, type="pearson")`
  returns Pearson residuals; `fitted(object)` returns fitted values;
  `hatvalues(object)` returns the leverages, or perhaps the value 1
  which will effectively ignore setting the hatvalues. In addition, the
  `data` argument should contain no missing values among the columns
  actually used in fitting the model.

A number of examples are available in Chapters 8 and 9 of [Modern
Statistics with R](https://www.modernstatisticswithr.com/).

Here are some simple examples with a linear regression model for the
`mtcars` data:

``` r
library(boot.pval)

# Bootstrap summary of a linear model for mtcars:
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model)
#>                Estimate Lower.bound Upper.bound p.value
#> (Intercept) 26.96300111  21.4837009 33.13272436  <0.001
#> hp          -0.05453412  -0.0847155 -0.02624529  <0.001
#> vs           2.57622314  -1.3354573  6.23773863   0.225

# Use 9999 bootstrap replicates and adjust p-values for
# multiplicity using Holm's method:
boot_summary(model, R = 9999, adjust.method = "holm")
#>                Estimate Lower.bound Upper.bound p.value Adjusted p-value
#> (Intercept) 26.96300111 21.29253243  32.6863499  <1e-04           0.0003
#> hp          -0.05453412 -0.08217782  -0.0256178   4e-04           0.0008
#> vs           2.57622314 -1.42802199   6.4742569  0.2012           0.2012

# Use case resampling instead of residual resampling:
boot_summary(model, method = "case")
#>                Estimate Lower.bound Upper.bound p.value
#> (Intercept) 26.96300111  21.5588473 35.07007975  <0.001
#> hp          -0.05453412  -0.1024433 -0.02828226  <0.001
#> vs           2.57622314  -1.4885872  6.50911625   0.256

# Export results to a gt table:
boot_summary(model, R = 9999) |>
  summary_to_gt()
```

|             | Estimate | 95 % CI          | p-value |
|-------------|----------|------------------|---------|
| (Intercept) | 26.963   | (21.397, 32.735) | \<1e-04 |
| hp          | −0.055   | (−0.083, −0.026) | \<1e-04 |
| vs          | 2.576    | (−1.408, 6.325)  | 0.1987  |

See Davison and Hinkley (1997) for details about residual resampling
(the default) and case resampling.

``` r
# Export results to a Word document:
library(flextable)
boot_summary(model, R = 9999) |>
  summary_to_flextable() |> 
  save_as_docx(path = "my_table.docx")
```

And a toy example for a generalised linear mixed model (using a small
number of bootstrap repetitions):

``` r
library(lme4)
model <- glmer(TICKS ~ YEAR + (1|LOCATION),
           data = grouseticks, family = poisson)
boot_summary(model, R = 99)
```

## Speeding up computations

For complex models, speed can be greatly improved by using
parallelisation. For `lmer` and `glmer` models, this is set using the
`parallel` (available options are `"multicore"` and `"snow"`). The
number of CPUs to use is set using `ncpus`.

``` r
model <- glmer(TICKS ~ YEAR + (1|LOCATION),
           data = grouseticks, family = poisson)
boot_summary(model, R = 999, parallel = "multicore", ncpus = 10)
```

For other models, use `ncores`:

``` r
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model, R = 9999, ncores = 10)
```

## Survival models

Survival regression models should be fitted using the argument
`model = TRUE`. A summary table can then be obtained using
`censboot_summary`. By default, the table contains exponentiated
coefficients (i.e. hazard ratios, in the case of a Cox PH model).

``` r
library(survival)
# Weibull AFT model:
model <- survreg(formula = Surv(time, status) ~ age + sex, data = lung,
                dist = "weibull", model = TRUE)
# Table with exponentiated coefficients:
censboot_summary(model)
#> Using exponentiated coefficients.
#>                Estimate Lower.bound Upper.bound p.value
#> (Intercept) 531.0483429  214.499098 1343.961015  <0.001
#> age           0.9878178    0.973827    1.001999   0.089
#> sex           1.4653368    1.168961    1.875647   0.003

# Cox PH model:
model <- coxph(formula = Surv(time, status) ~ age + sex, data = lung,
               model = TRUE)
# Table with hazard ratios:
censboot_summary(model)
#> Using exponentiated coefficients.
#>     Estimate Lower.bound Upper.bound p.value
#> age 1.017191   0.9996458    1.035802   0.054
#> sex 0.598566   0.4303380    0.815483  <0.001
# Table with original coefficients:
censboot_summary(model, coef = "raw")
#> Using raw coefficients.
#>        Estimate  Lower.bound Upper.bound p.value
#> age  0.01704533 -0.002462382  0.03643551   0.086
#> sex -0.51321852 -0.850467317 -0.20129079   0.005
```

To speed up computations using parallelisation, use the `parallel` and
`ncpus` arguments:

``` r
censboot_summary(model, parallel = "multicore", ncpus = 10)
```

## References

- Davison, A.C. and Hinkley, D.V. (1997) *Bootstrap Methods and Their
  Application*. Cambridge University Press.
- Hall P (1992). *The Bootstrap and Edgeworth Expansion*. Springer, New
  York.
- Thulin, M. (2024) *Modern Statistics with R*. Second edition. Chapman
  & Hall/CRC Press.
