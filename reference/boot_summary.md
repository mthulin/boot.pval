# Summarising Regression Models Using the Bootstrap

Summaries for regression models, including "lm", "glm", "glm.nb", "nls",
"rlm", "polr", and "merMod" ("lmer", "glmer") objects, using the
bootstrap for p-values and confidence intervals.

## Usage

``` r
boot_summary(
  model,
  type = "perc",
  method = NULL,
  conf.level = 0.95,
  R = 999,
  coef = "raw",
  pval_precision = NULL,
  adjust.method = "none",
  ...
)
```

## Arguments

- model:

  An object fitted using e.g. "lm", "glm", "glm.nb", "nls", "rlm",
  "polr", lmer", or "glmer".

- type:

  A vector of character strings representing the type of interval to
  base the test on. The value should be one of "norm", "basic", "bca",
  and "perc" (the default). "bca" is not supported for "lmer" and
  "glmer" models.

- method:

  The method used for bootstrapping. For "lm" and "nls" objects use
  either "residual" (for resampling of scaled and centred residuals, the
  default) or "case" (for case resampling). For "glm" objects, use
  "case" (the default). For "merMod" objects (mixed models) use either
  "parametric" (the default) or "semiparametric".

- conf.level:

  The confidence level for the confidence intervals. The default is
  0.95.

- R:

  The number of bootstrap replicates. The default is 999.

- coef:

  A string specifying whether to use exponentiated coefficients in the
  summary table. Either "exp" (for exponentiated coefficients, i.e. odds
  ratios in the case of a logistic regression model) or "raw" (for
  coefficients on their original scale). The default is "raw", which is
  recommended for linear models.

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

  Additional arguments passed to `Boot` or `bootMer`, such as `parallel`
  for parallel computations. See
  [`?car::Boot`](https://rdrr.io/pkg/car/man/Boot.html) and
  [`?lme4::bootMer`](https://rdrr.io/pkg/lme4/man/bootMer.html) for
  details.

## Value

A data frame containing coefficient estimates, bootstrap confidence
intervals, and bootstrap p-values.

## Details

p-values can be computed by inverting the corresponding confidence
intervals, as described in Section 14.2 of Thulin (2024) and Section
3.12 in Hall (1992). This function computes p-values for coefficients of
regression models in this way. The approach relies on the fact that:

- the p-value of the two-sided test for the parameter theta is the
  smallest alpha such that theta is not contained in the corresponding
  1-alpha confidence interval,

- for a test of the parameter theta with significance level alpha, the
  set of values of theta that aren't rejected by the two-sided test
  (when used as the null hypothesis) is a 1-alpha confidence interval
  for theta.

The function can be used with "lm", "glm", "glm.nb", "nls", "rlm", and
"merMod" ("lmer", "glmer") objects. In addition, it should work for any
regression model such that: `residuals(object, type="pearson")` returns
Pearson residuals; `fitted(object)` returns fitted values;
`hatvalues(object)` returns the leverages, or perhaps the value 1 which
will effectively ignore setting the hatvalues. In addition, the `data`
argument should contain no missing values among the columns actually
used in fitting the model.

## References

Hall P (1992). *The Bootstrap and Edgeworth Expansion*. Springer, New
York. ISBN 9781461243847. Thulin M (2024). *Modern Statistics with R*.
Chapman & Hall/CRC Press, Boca Raton. ISBN 9781032512440,
<https://www.modernstatisticswithr.com/>.

## See also

[`boot_t_test()`](https://mthulin.github.io/boot.pval/reference/boot_t_test.md)
for bootstrap t-tests,
[`boot_median_test()`](https://mthulin.github.io/boot.pval/reference/boot_median_test.md)
for bootstrap tests for medians.

## Examples

``` r
# Bootstrap summary of a linear model for mtcars:
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model, R = 99)
#>                Estimate Lower.bound Upper.bound p.value
#> (Intercept) 26.96300111 22.46858043 33.04643462   <0.01
#> hp          -0.05453412 -0.08515078 -0.03300863   <0.01
#> vs           2.57622314 -1.62537273  6.29452416    0.37
# (Values for R greater than 99 are recommended for most applications.)

# Adjust p-values for multiplicity using Holm's method:
boot_summary(model, R = 99, adjust.method = "holm")
#>                Estimate Lower.bound Upper.bound p.value Adjusted p-value
#> (Intercept) 26.96300111 20.64731661 33.01217023   <0.01             0.03
#> hp          -0.05453412 -0.08505702 -0.02204617   <0.01             0.03
#> vs           2.57622314 -1.10406035  7.69178446    0.23             0.23
```
