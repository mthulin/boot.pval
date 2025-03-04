# boot.pval <img src="man/figures/logo.png" align="right" width="350" />
This R package provides functions for computing bootstrap p-values based on `boot` objects, and convenience functions for bootstrap confidence intervals and p-values for various regression models.

## Installation
To install the package from CRAN:

```
install.packages("boot.pval")
```

To install the development version from Github:

```
library(devtools)
install_github("mthulin/boot.pval")
```

## Background
p-values can be computed by inverting the corresponding confidence intervals, as described in Section 14.2 of [Thulin (2024)](https://www.modernstatisticswithr.com/mathschap.html#confintequal) and Section 3.12 in [Hall (1992)](https://link.springer.com/book/10.1007/978-1-4612-4384-7). This package contains functions for computing bootstrap p-values in this way. The approach relies on the fact that:

- The p-value of the two-sided test for the parameter theta is the smallest alpha such that theta is not contained in the corresponding 1-alpha confidence interval,
- For a test of the parameter theta with significance level alpha, the set of values of theta that aren't rejected by the two-sided test (when used as the null hypothesis) is a 1-alpha confidence interval for theta.


## Summaries for regression models
Summary tables with confidence intervals and p-values for the coefficients of regression models can be obtained using the `boot_summary` (most models) and `censboot_summary` (models with censored response variables) functions. Currently, the following models are supported:

- Linear models fitted using `lm`,
- Generalised linear models fitted using `glm` or `glm.nb`,
- Nonlinear models fitted using `nls`,
- Robust linear models fitted using `MASS::rlm`,
- Ordered logistic or probit regression models fitted (without weights) using `MASS:polr`,
- Linear mixed models fitted using `lme4::lmer` or `lmerTest::lmer`,
- Generalised linear mixed models fitted using `lme4::glmer`.
- Cox PH regression models fitted using `survival::coxph` (using `censboot_summary`).
- Accelerated failure time models fitted using `survival::survreg` or `rms::psm` (using `censboot_summary`).
- Any regression model such that: `residuals(object, type="pearson")` returns Pearson residuals; `fitted(object)` returns fitted values; `hatvalues(object)` returns the leverages, or perhaps the value 1 which will effectively ignore setting the hatvalues. In addition, the `data` argument should contain no missing values among the columns actually used in fitting the model.

A number of examples are available in Chapters 8 and 9 of [Modern Statistics with R](https://www.modernstatisticswithr.com/).

Here are some simple examples with a linear regression model for the `mtcars` data:

```r
# Bootstrap summary of a linear model for mtcars:
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model)

# Use 9999 bootstrap replicates and adjust p-values for
# multiplicity using Holm's method:
boot_summary(model, R = 9999, adjust.method = "holm")

# Use case resampling instead of residual resampling:
boot_summary(model, method = "case")

# Export results to a gt table:
boot_summary(model, R = 9999) |>
  summary_to_gt()
```

See Davison and Hinkley (1997) for details about residual resampling (the default) and case resampling.

```r
# Export results to a Word document:
library(flextable)
boot_summary(model, R = 9999) |>
  summary_to_flextable() |> 
  save_as_docx(path = "my_table.docx")
```

And a toy example for a generalised linear mixed model (using a small number of bootstrap repetitions):

```r
library(lme4)
model <- glmer(TICKS ~ YEAR + (1|LOCATION),
           data = grouseticks, family = poisson)
boot_summary(model, R = 99)
```

## Speeding up computations
For complex models, speed can be greatly improved by using parallelisation. For `lmer` and `glmer` models, this is set using the `parallel` (available options are `"multicore"` and `"snow"`). The number of CPUs to use is set using `ncpus`.

```r
model <- glmer(TICKS ~ YEAR + (1|LOCATION),
           data = grouseticks, family = poisson)
boot_summary(model, R = 999, parallel = "multicore", ncpus = 10)
```

For other models, use `ncores`:

```r
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model, R = 9999, ncores = 10)
```

## Survival models
Survival regression models should be fitted using the argument `model = TRUE`. A summary table can then be obtained using `censboot_summary`. By default, the table contains exponentiated coefficients (i.e. hazard ratios, in the case of a Cox PH model).

```r
library(survival)
# Weibull AFT model:
model <- survreg(formula = Surv(time, status) ~ age + sex, data = lung,
                dist = "weibull", model = TRUE)
# Table with exponentiated coefficients:
censboot_summary(model)

# Cox PH model:
model <- coxph(formula = Surv(time, status) ~ age + sex, data = lung,
               model = TRUE)
# Table with hazard ratios:
censboot_summary(model)
# Table with original coefficients:
censboot_summary(model, coef = "raw")
```

To speed up computations using parallelisation, use the `parallel` and `ncpus` arguments:

```r
censboot_summary(model, parallel = "multicore", ncpus = 10)
```

## Tests of locations
Traditional versions of Student's t-test (`t.test` in R) rely on the assumption of normality. For non-normal data, this can lead to misleading p-values and confidence intervals. In such cases, it is often recommended to use the Wilcoxon-Mann-Whitney test (`wilcox.test` in R) instead. Despite being described as a test of location, or a test for differences of medians, the Wilcoxon-Mann-Whitney test is actually a test of equivalence of distributions, unless strict assumptions are met. In addition, `wilcox.test` does not provide a confidence interval for the difference of the medians.

In many cases, a better option is to use a bootstrap t-test (for inference about means) or a bootstrap median test (for inference about medians). These can be used without the normality assumption, and will provide confidence intervals for the parameters of interest.

### Two-sample bootstrap t-tests
To illustrate the use of bootstrap t-tests, we'll use the classic `sleep` data, which "show the effect of two soporific drugs (increase in hours of sleep compared to control) on 10 patients" (see `?sleep` for details).

We wish to test whether the mean value of the `extra` (increase in hours of sleep) variable differs between the two groups described by the `group` variable. The syntax for this is identical to that for `t.test`:

```r
boot_t_test(extra ~ group, data = sleep)
```

If you prefer, you can also use the `|>` pipe as follows:

```r
sleep |> boot_t_test(extra ~ group)
```

By default, the confidence interval and p-value are based on the studentized bootstrap confidence interval. Other options available are normal, basic, percentile and BCa intervals; see Chapter 5 of Davison and Hinkley (1997) for details. You can choose the method used using the `type` argument.

```r
sleep |> boot_t_test(extra ~ group, type = "perc") # Percentile interval
sleep |> boot_t_test(extra ~ group, type = "bca") # BCa interval
```

You can control the number of bootstrap replicates used (argument `R`; the default is 9999) or the direction of the alternative hypothesis (argument `alternative`):

```r
sleep |> boot_t_test(extra ~ group, R = 999, alternative = "less")
```

In this case, the data is actually paired, so it would make sense to perform a paired bootstrap t-test instead. We reshape the data to a wide format, so that the first measurements ends up in the variable `extra.1`, and the second measurement ends up in the variable `extra.2`. We can then run the test as follows:

```r
# Reshape to wide format:
sleep2 <- reshape(sleep, direction = "wide",
                  idvar = "ID", timevar = "group")

# Traditional interface:
boot_t_test(sleep2$extra.1, sleep2$extra.2, paired = TRUE)

# Using pipes:
sleep2 |> boot_t_test(Pair(extra.1, extra.2) ~ 1)
```

### One-sample bootstrap t-tests
For one sample bootstrap t-tests, we only need to provide a single vector containing the measurements. We can also specify the null value of the mean (argument `mu`):

```r
# Traditional interface:
boot_t_test(sleep$extra, mu = 1)

# Using pipes:
sleep |> boot_t_test(extra ~ 1, mu = 1)
```

### Bootstrap median tests
Running a bootstrap median test with the `boot_median_test` function is completely analogously to running a bootstrap t-test. The only difference is under the hood - medians are used instead of means. Because the studentized and BCa versions of this test use an inner bootstrap to estimate the variance of the statistic, these takes longer to run than other tests presented here.

```r
boot_median_test(extra ~ group, data = sleep, type = "perc")

sleep |> boot_median_test(extra ~ group, R = 999, alternative = "less")

boot_median_test(sleep$extra, mu = 1)
```



## Other hypothesis tests
Bootstrap p-values for hypothesis tests based on `boot` objects can be obtained using the `boot.pval` function. The following examples are extensions of those given in the documentation for `boot::boot`:

```r
# Hypothesis test for the city data
# H0: ratio = 1.4
library(boot)
ratio <- function(d, w) sum(d$x * w)/sum(d$u * w)
city.boot <- boot(city, ratio, R = 999, stype = "w", sim = "ordinary")
boot.pval(city.boot, theta_null = 1.4)

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
grav1.boot <- boot(grav1, diff.means, R = 999, stype = "f",
                   strata = grav1[ ,2])
boot.pval(grav1.boot, type = "stud", theta_null = 0)
```

## References
* Davison, A.C. and Hinkley, D.V. (1997) _Bootstrap Methods and Their Application_. Cambridge University Press.
* Thulin, M. (2024) _Modern Statistics with R_. Second edition. Chapman & Hall/CRC Press.
