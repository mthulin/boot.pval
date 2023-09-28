# boot.pval
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
p-values can be computed by inverting the corresponding confidence intervals, as described in Section 12.2 of [Thulin (2021)](https://www.modernstatisticswithr.com/mathschap.html#confintequal) and Section 3.12 in [Hall (1992)](https://link.springer.com/book/10.1007/978-1-4612-4384-7). This package contains functions for computing bootstrap p-values in this way. The approach relies on the fact that:

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

Here is an example with a linear regression model for the `mtcars` data:

```
# Bootstrap summary of a linear model for mtcars:
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model)

# Use 9999 bootstrap replicates and adjust p-values for
# multiplicity using Holm's method:
boot_summary(model, R = 9999, adjust.method = "holm")
```

And a toy example for a generalised linear mixed model (using a small number of bootstrap repetitions):

```
library(lme4)
model <- glmer(TICKS ~ YEAR + (1|LOCATION),
           data = grouseticks, family = poisson)
boot_summary(model, R = 99)
```

Survival regression models should be fitted using the argument `model = TRUE`. A summary table can then be obtained using `censboot_summary`. By default, the table contains exponentiated coefficients (i.e. hazard ratios, in the case of a Cox PH model).

```
library(survival)
# Weibull AFT model:
model <- survreg(formula = Surv(time, status) ~ age + sex, data = lung,
                dist = "weibull", model = TRUE)
censboot_summary(model)

# Cox PH model:
model <- coxph(formula = Surv(time, status) ~ age + sex, data = lung,
               model = TRUE)
# Table with hazard ratios:
censboot_summary(model)
# Table with original coefficients:
censboot_summary(model, coef = "raw")
```

## Other hypothesis tests
Bootstrap p-values for hypothesis tests based on `boot` objects can be obtained using the `boot.pval` function. The following examples are extensions of those given in the documentation for `boot::boot`:

```
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
