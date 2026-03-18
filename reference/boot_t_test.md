# Bootstrap t-Test

Performs one- and two-sample bootstrap t-tests and computes the
corresponding bootstrap confidence interval.

## Usage

``` r
boot_t_test(x, ...)

# Default S3 method
boot_t_test(
  x,
  y = NULL,
  alternative = c("two.sided", "less", "greater"),
  mu = 0,
  paired = FALSE,
  var.equal = FALSE,
  conf.level = 0.95,
  R = 9999,
  type = "stud",
  ...
)

# S3 method for class 'formula'
boot_t_test(formula, data, subset, na.action, ...)

# S3 method for class 'data.frame'
boot_t_test(x, formula, ...)

# S3 method for class 'matrix'
boot_t_test(x, formula, ...)
```

## Arguments

- x:

  a (non-empty) numeric vector of data values.

- ...:

  Additional arguments passed to `boot`, such as `parallel` for parallel
  computations. See
  [`?boot::boot`](https://rdrr.io/pkg/boot/man/boot.html) for details.

- y:

  an optional (non-empty) numeric vector of data values.

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `"two.sided"` (default), `"greater"` or `"less"`. You can specify
  just the initial letter.

- mu:

  a number indicating the true value of the mean (or difference in means
  if you are performing a two sample test).

- paired:

  a logical indicating whether you want a paired t-test.

- var.equal:

  a logical variable indicating whether to treat the two variances as
  being equal. If `TRUE` then the pooled variance is used to estimate
  the variance otherwise the Welch (or Satterthwaite) approximation to
  the degrees of freedom is used.

- conf.level:

  confidence level of the interval.

- R:

  The number of bootstrap replicates. The default is 9999.

- type:

  A vector of character strings representing the type of interval to
  base the test on. The value should be one of "norm", "basic", "bca",
  perc", and "stud" (the default).

- formula:

  a formula of the form `lhs ~ rhs` where `lhs` is a numeric variable
  giving the data values and `rhs` either `1` for a one-sample or paired
  test or a factor with two levels giving the corresponding groups. If
  `lhs` is of class `"`[`Pair`](https://rdrr.io/r/stats/Pair.html)`"`
  and `rhs` is `1`, a paired test is done, see Examples.

- data:

  an optional matrix or data frame (or similar: see
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html)) containing
  the variables in the formula `formula`. By default the variables are
  taken from `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  [`NA`](https://rdrr.io/r/base/NA.html)s.

## Value

A list with class `"htest"`) containing the following components:

- statistic:

  the value of the t-statistic.

- R:

  the number of bootstrap replicates used.

- p.value:

  the bootstrap p-value for the test.

- conf.int:

  a bootstrap confidence interval for the mean appropriate to the
  specified alternative hypothesis.

- estimate:

  the estimated mean or difference in means depending on whether it was
  a one-sample test or a two-sample test.

- null.value:

  the specified hypothesized value of the mean or mean difference
  depending on whether it was a one-sample test or a two-sample test.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string indicating what type of t-test was performed.

- data.name:

  a character string giving the name(s) of the data.

## Details

p-values can be computed by inverting the corresponding confidence
intervals, as described in Section 14.2 of Thulin (2024) and Section
3.12 in Hall (1992). This function computes p-values for the t-test in
this way. The approach relies on the fact that:

- the p-value of the two-sided test for the parameter theta is the
  smallest alpha such that theta is not contained in the corresponding
  1-alpha confidence interval,

- for a test of the parameter theta with significance level alpha, the
  set of values of theta that aren't rejected by the two-sided test
  (when used as the null hypothesis) is a 1-alpha confidence interval
  for theta. Consequently, the p-value will be consistent with the
  confidence interval, in the sense that the null hypothesis is rejected
  if and only if the null parameter values is not contained in the
  confidence interval.

## References

Hall P (1992). *The Bootstrap and Edgeworth Expansion*. Springer, New
York. ISBN 9781461243847. Thulin M (2024). *Modern Statistics with R*.
Chapman & Hall/CRC Press, Boca Raton. ISBN 9781032512440,
<https://www.modernstatisticswithr.com/>.

## See also

[`boot_median_test()`](https://mthulin.github.io/boot.pval/reference/boot_median_test.md)
for bootstrap tests for medians,
[`boot_summary()`](https://mthulin.github.io/boot.pval/reference/boot_summary.md)
for bootstrap tests for coefficients of regression models.

## Examples

``` r
# Generate example data:
# x is the variable of interest
# y is the grouping variable
example_data <- data.frame(x = rnorm(40), y = rep(c(1,2), 20))

# Two-sample (Welch) test:
boot_t_test(x ~ y, data = example_data, R = 999)
#> 
#>  Welch Two Sample Bootstrap t-test (studentized)
#> 
#> data:  x by y
#> t = 0.99842, R = 999, p-value = 0.3123
#> alternative hypothesis: true difference in means between group 1 and group 2 is not equal to 0
#> 95 percent confidence interval:
#>  -0.3611907  1.1025029
#> sample estimates:
#> mean in group 1 mean in group 2 
#>      0.08955452     -0.26890631 
#> 

# Two-sample (Welch) test using the pipe:
example_data |> boot_t_test(x ~ y, R = 999)
#> 
#>  Welch Two Sample Bootstrap t-test (studentized)
#> 
#> data:  x by y
#> t = 0.99842, R = 999, p-value = 0.3824
#> alternative hypothesis: true difference in means between group 1 and group 2 is not equal to 0
#> 95 percent confidence interval:
#>  -0.3840846  1.0607978
#> sample estimates:
#> mean in group 1 mean in group 2 
#>      0.08955452     -0.26890631 
#> 

# With a directed alternative hypothesis:
example_data |> boot_t_test(x ~ y, R = 999, alternative = "greater")
#> 
#>  Welch Two Sample Bootstrap t-test (studentized)
#> 
#> data:  x by y
#> t = 0.99842, R = 999, p-value = 0.1552
#> alternative hypothesis: true difference in means between group 1 and group 2 is greater than 0
#> 95 percent confidence interval:
#>  -0.2215106        Inf
#> sample estimates:
#> mean in group 1 mean in group 2 
#>      0.08955452     -0.26890631 
#> 

# One-sample test:
boot_t_test(example_data$x, R = 999)
#> 
#>  One Sample Bootstrap t-test (studentized)
#> 
#> data:  example_data$x
#> t = -0.49957, R = 999, p-value = 0.6587
#> alternative hypothesis: true mean is not equal to 0
#> 95 percent confidence interval:
#>  -0.4128478  0.2713103
#> sample estimates:
#>  mean of x 
#> -0.0896759 
#> 

# One-sample test using the pipe:
example_data |> boot_t_test(x ~ 1, R = 999)
#> 
#>  One Sample Bootstrap t-test (studentized)
#> 
#> data:  x
#> t = -0.49957, R = 999, p-value = 0.6106
#> alternative hypothesis: true mean is not equal to 0
#> 95 percent confidence interval:
#>  -0.4808734  0.2779611
#> sample estimates:
#>  mean of x 
#> -0.0896759 
#> 

# With a directed alternative hypothesis:
example_data |> boot_t_test(x ~ 1, R = 999, mu = 0.5, alternative = "less")
#> 
#>  One Sample Bootstrap t-test (studentized)
#> 
#> data:  x
#> t = -3.285, R = 999, p-value = 0.002002
#> alternative hypothesis: true mean is less than 0.5
#> 95 percent confidence interval:
#>      -Inf 0.183044
#> sample estimates:
#>  mean of x 
#> -0.0896759 
#> 

# Paired test:
boot_t_test(example_data$x[example_data$y==1],
            example_data$x[example_data$y==2],
            paired = TRUE, R = 999)
#> 
#>  Paired Bootstrap t-test (studentized)
#> 
#> data:  example_data$x[example_data$y == 1] and example_data$x[example_data$y == 2]
#> t = 0.98995, R = 999, p-value = 0.3233
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -0.3972185  1.1647955
#> sample estimates:
#> mean of the differences 
#>               0.3584608 
#> 

# Paired test using the pipe (after reshaping to wide format):
example_data$id <- rep(1:20, rep(2, 20))
example_data2 <- reshape(example_data, direction = "wide",
                         idvar = "id", timevar = "y")
example_data2 |> boot_t_test(Pair(x.1, x.2) ~ 1)
#> 
#>  Paired Bootstrap t-test (studentized)
#> 
#> data:  Pair(x.1, x.2)
#> t = 0.98995, R = 9999, p-value = 0.3102
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -0.3584844  1.1611572
#> sample estimates:
#> mean of the differences 
#>               0.3584608 
#> 
```
