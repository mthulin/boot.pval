# Summarising Survival Regression Models Using the Bootstrap

Summaries for Cox proportional hazards and accelerated failure time
models, using the bootstrap for p-values and confidence intervals.

## Usage

``` r
censboot_summary(
  model,
  type = "perc",
  sim = "ordinary",
  strata = NULL,
  coef = "exp",
  conf.level = 0.95,
  R = 999,
  pval_precision = NULL,
  adjust.method = "none",
  ...
)
```

## Arguments

- model:

  An object fitted using "survival::coxph", "survival::survreg", or
  "rms::psm".

- type:

  A vector of character strings representing the type of interval to
  base the test on. The value should be one of "norm", "basic", and
  "perc" (the default).

- sim:

  The method used for bootstrapping. See
  [`?boot::censboot`](https://rdrr.io/pkg/boot/man/censboot.html) for
  details. Currently only "ordinary" (case resampling) is supported.

- strata:

  The strata used in the calls to `survfit.` It can be a vector or a
  matrix with 2 columns. If it is a vector then it is assumed to be the
  strata for the survival distribution, and the censoring distribution
  is assumed to be the same for all observations. If it is a matrix then
  the first column is the strata for the survival distribution and the
  second is the strata for the censoring distribution. When
  `sim = "ordinary"`, only one set of strata is used to stratify the
  observations. This is taken to be the first column of `strata` when it
  is a matrix.

- coef:

  A string specifying whether to use exponentiated coefficients in the
  summary table. Either "exp" (for exponentiated coefficients, i.e.
  hazard ratios in the case of a Cox PH model) or "raw" (for
  coefficients on their original scale). The default is "exp".

- conf.level:

  The confidence level for the confidence intervals. The default is
  0.95.

- R:

  The number of bootstrap replicates. The default is 999.

- pval_precision:

  The desired precision for the p-value. The default is 1/R.

- adjust.method:

  Adjustment of p-values for multiple comparisons using `p.adjust`. The
  default is "none", in which case the p-values aren't adjusted. The
  other options are "holm", "hochberg", "hommel", "bonferroni", "BH",
  "BY", and "fdr"; see
  [`?p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for details on
  these methods.

- ...:

  Additional arguments passed to `censboot`, such as `parallel` for
  parallel computations. See
  [`?boot::censboot`](https://rdrr.io/pkg/boot/man/censboot.html) for
  details.

## Value

A data frame containing coefficient estimates, bootstrap confidence
intervals, and bootstrap p-values.

## Details

p-values can be computed by inverting the corresponding confidence
intervals, as described in Section 14.2 of Thulin (2024) and Section
3.12 in Hall (1992). This function computes p-values in this way from
"coxph" or "survreg" objects. The approach relies on the fact that:

- the p-value of the two-sided test for the parameter theta is the
  smallest alpha such that theta is not contained in the corresponding
  1-alpha confidence interval,

- for a test of the parameter theta with significance level alpha, the
  set of values of theta that aren't rejected by the two-sided test
  (when used as the null hypothesis) is a 1-alpha confidence interval
  for theta.

## References

Hall P (1992). *The Bootstrap and Edgeworth Expansion*. Springer, New
York. ISBN 9781461243847. Thulin M (2024). *Modern Statistics with R*.
Chapman & Hall/CRC Press, Boca Raton. ISBN 9781032512440,
<https://www.modernstatisticswithr.com/>.

## Examples

``` r
library(survival)
#> 
#> Attaching package: ‘survival’
#> The following object is masked from ‘package:boot’:
#> 
#>     aml
# Weibull AFT model:
# Note that model = TRUE is required for use with censboot_summary:
model <- survreg(formula = Surv(time, status) ~ age + sex, data = lung,
                 dist = "weibull", model = TRUE)
censboot_summary(model, R = 99)
#> Using exponentiated coefficients.
#>                Estimate Lower.bound Upper.bound p.value
#> (Intercept) 531.0483429 185.4615956 1309.466858   <0.01
#> age           0.9878178   0.9746921    1.002980    0.17
#> sex           1.4653368   1.1832707    1.913004   <0.01
# (Values for R greater than 99 are recommended for most applications.)

# Cox PH model:
model <- coxph(formula = Surv(time, status) ~ age + sex, data = lung,
               model = TRUE)
# Table with hazard ratios:
censboot_summary(model, R = 99)
#> Using exponentiated coefficients.
#>     Estimate Lower.bound Upper.bound p.value
#> age 1.017191   0.9986972   1.0382557    0.08
#> sex 0.598566   0.4022327   0.8631437   <0.01
censboot_summary(model, coef = "raw", R = 99)
#> Using raw coefficients.
#>        Estimate  Lower.bound Upper.bound p.value
#> age  0.01704533  0.001207305  0.03921424   <0.01
#> sex -0.51321852 -0.946193565 -0.24131217   <0.01
```
