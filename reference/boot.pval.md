# Compute Bootstrap p-values

Compute bootstrap p-values through confidence interval inversion, as
described in Hall (1992) and Thulin (2024).

## Usage

``` r
boot.pval(
  boot_res,
  type = "perc",
  theta_null = 0,
  pval_precision = NULL,
  alternative = "two.sided",
  ...
)
```

## Arguments

- boot_res:

  An object of class "boot" containing the output of a bootstrap
  calculation.

- type:

  A vector of character strings representing the type of interval to
  base the test on. The value should be one of "norm", "basic", "stud",
  "perc" (the default), and "bca".

- theta_null:

  The value of the parameter under the null hypothesis.

- pval_precision:

  The desired precision for the p-value. The default is 1/R, where R is
  the number of bootstrap samples in `boot_res`.

- alternative:

  A character string specifying the alternative hypothesis. Must be one
  of "two.sided" (default), "greater", or "less".

- ...:

  Additional arguments passed to `boot.ci`.

## Value

A bootstrap p-value.

## Details

p-values can be computed by inverting the corresponding confidence
intervals, as described in Section 14.2 of Thulin (2024) and Section
3.12 in Hall (1992). This function computes p-values in this way from
"boot" objects. The approach relies on the fact that:

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

## See also

[`boot_t_test()`](https://mthulin.github.io/boot.pval/reference/boot_t_test.md)
for bootstrap t-tests,
[`boot_median_test()`](https://mthulin.github.io/boot.pval/reference/boot_median_test.md)
for bootstrap tests for medians,
[`boot_summary()`](https://mthulin.github.io/boot.pval/reference/boot_summary.md)
for bootstrap tests for coefficients of regression models.

## Examples

``` r
# Hypothesis test for the city data
# H0: ratio = 1.4
library(boot)
ratio <- function(d, w) sum(d$x * w)/sum(d$u * w)
city.boot <- boot(city, ratio, R = 99, stype = "w", sim = "ordinary")
boot.pval(city.boot, theta_null = 1.4)
#> [1] 0.3838384

# Studentized test for the two sample difference of means problem
# using the final two series of the gravity data.
diff.means <- function(d, f)
{
  n <- nrow(d)
  gp1 <- 1:table(as.numeric(d$series))[1]
  m1 <- sum(d[gp1,1] * f[gp1])/sum(f[gp1])
  m2 <- sum(d[-gp1,1] * f[-gp1])/sum(f[-gp1])
  ss1 <- sum(d[gp1,1]^2 * f[gp1]) - (m1 *  m1 * sum(f[gp1]))
  ss2 <- sum(d[-gp1,1]^2 * f[-gp1]) - (m2 *  m2 * sum(f[-gp1]))
  c(m1 - m2, (ss1 + ss2)/(sum(f) - 2))
}
grav1 <- gravity[as.numeric(gravity[,2]) >= 7, ]
grav1.boot <- boot(grav1, diff.means, R = 99, stype = "f",
                   strata = grav1[ ,2])
boot.pval(grav1.boot, type = "stud", theta_null = 0)
#> [1] 0.01010101
```
