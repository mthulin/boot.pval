library(boot)

city_boot <- boot(
  city,
  function(d, w) sum(d$x * w) / sum(d$u * w),
  R = 99,
  stype = "w",
  sim = "ordinary"
)

test_that("boot.pval returns a single numeric between 0 and 1", {
  result <- boot.pval(city_boot, theta_null = 1.4)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("boot.pval works with norm, basic, and perc CI types", {
  for (type in c("norm", "basic", "perc")) {
    result <- boot.pval(city_boot, type = type, theta_null = 1.4)
    expect_type(result, "double")
    expect_length(result, 1)
    expect_gte(result, 0)
    expect_lte(result, 1)
  }
})

test_that("boot.pval works with one-sided alternatives", {
  result_less <- boot.pval(city_boot, theta_null = 1.4, alternative = "less")
  result_greater <- boot.pval(city_boot, theta_null = 1.4, alternative = "greater")

  for (result in list(result_less, result_greater)) {
    expect_type(result, "double")
    expect_length(result, 1)
    expect_gte(result, 0)
    expect_lte(result, 1)
  }
})

test_that("boot.pval errors on invalid alternative", {
  expect_snapshot(
    boot.pval(city_boot, alternative = "invalid"),
    error = TRUE
  )
})

test_that("boot.pval stud type works with a studentized boot object", {
  diff_means <- function(d, f) {
    n <- nrow(d)
    gp1 <- 1:table(as.numeric(d$series))[1]
    m1 <- sum(d[gp1, 1] * f[gp1]) / sum(f[gp1])
    m2 <- sum(d[-gp1, 1] * f[-gp1]) / sum(f[-gp1])
    ss1 <- sum(d[gp1, 1]^2 * f[gp1]) - (m1 * m1 * sum(f[gp1]))
    ss2 <- sum(d[-gp1, 1]^2 * f[-gp1]) - (m2 * m2 * sum(f[-gp1]))
    c(m1 - m2, (ss1 + ss2) / (sum(f) - 2))
  }
  grav1 <- gravity[as.numeric(gravity[, 2]) >= 7, ]
  grav1_boot <- boot(grav1, diff_means, R = 99, stype = "f", strata = grav1[, 2])
  result <- boot.pval(grav1_boot, type = "stud", theta_null = 0)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_gte(result, 0)
  expect_lte(result, 1)
})
