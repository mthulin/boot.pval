x <- c(2.1, 3.4, 1.9, 4.2, 2.8, 3.7, 2.5, 3.1)
y <- c(1.5, 2.8, 1.2, 3.1, 2.2, 2.9, 1.8, 2.4)
df_two <- data.frame(val = c(x, y), grp = rep(c(1, 2), each = length(x)))
df_one <- data.frame(x = x)

expected_htest_names <- c(
  "statistic", "parameters", "p.value", "conf.int",
  "estimate", "null.value", "alternative", "method", "data.name"
)

# ---- boot_t_test ----

test_that("boot_t_test one-sample returns an htest with correct structure", {
  result <- boot_t_test(x, R = 99)
  expect_s3_class(result, "htest")
  expect_setequal(names(result), expected_htest_names)
  expect_length(result$estimate, 1)
  expect_length(result$conf.int, 2)
  expect_equal(attr(result$conf.int, "conf.level"), 0.95)
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("boot_t_test two-sample returns an htest with two estimates", {
  result <- boot_t_test(x, y, R = 99)
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 2)
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("boot_t_test paired returns an htest with one estimate", {
  result <- boot_t_test(x, y, paired = TRUE, R = 99)
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 1)
})

test_that("boot_t_test formula interface works for two-sample test", {
  result <- boot_t_test(val ~ grp, data = df_two, R = 99)
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 2)
})

test_that("boot_t_test formula interface works for one-sample test", {
  result <- boot_t_test(x ~ 1, data = df_one, R = 99)
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 1)
})

test_that("boot_t_test data.frame interface works", {
  result <- df_two |> boot_t_test(val ~ grp, R = 99)
  expect_s3_class(result, "htest")
})

test_that("boot_t_test one-sided alternatives return correct alternative string", {
  for (alt in c("less", "greater")) {
    result <- boot_t_test(x, mu = 2, R = 99, alternative = alt)
    expect_equal(result$alternative, alt)
    expect_gte(result$p.value, 0)
    expect_lte(result$p.value, 1)
  }
})

test_that("boot_t_test works with norm, basic, and perc CI types", {
  for (type in c("norm", "basic", "perc")) {
    result <- boot_t_test(x, R = 99, type = type)
    expect_s3_class(result, "htest")
    expect_gte(result$p.value, 0)
    expect_lte(result$p.value, 1)
  }
})

test_that("boot_t_test conf.level is reflected in conf.int attribute", {
  result <- boot_t_test(x, R = 99, conf.level = 0.90)
  expect_equal(attr(result$conf.int, "conf.level"), 0.90)
})

# ---- boot_median_test ----
# Use type = "perc" to avoid the slow nested bootstrap in the default "stud" type.

test_that("boot_median_test one-sample returns an htest with correct structure", {
  result <- boot_median_test(x, R = 99, type = "perc")
  expect_s3_class(result, "htest")
  expect_setequal(names(result), expected_htest_names)
  expect_length(result$estimate, 1)
  expect_length(result$conf.int, 2)
  expect_equal(attr(result$conf.int, "conf.level"), 0.95)
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("boot_median_test two-sample returns an htest with two estimates", {
  result <- boot_median_test(x, y, R = 99, type = "perc")
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 2)
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("boot_median_test paired returns an htest with one estimate", {
  result <- boot_median_test(x, y, paired = TRUE, R = 99, type = "perc")
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 1)
})

test_that("boot_median_test formula interface works for two-sample test", {
  result <- boot_median_test(val ~ grp, data = df_two, R = 99, type = "perc")
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 2)
})

test_that("boot_median_test formula interface works for one-sample test", {
  result <- boot_median_test(x ~ 1, data = df_one, R = 99, type = "perc")
  expect_s3_class(result, "htest")
  expect_length(result$estimate, 1)
})

test_that("boot_median_test data.frame interface works", {
  result <- df_two |> boot_median_test(val ~ grp, R = 99, type = "perc")
  expect_s3_class(result, "htest")
})

test_that("boot_median_test one-sided alternatives return correct alternative string", {
  for (alt in c("less", "greater")) {
    result <- boot_median_test(x, mu = 2, R = 99, type = "perc", alternative = alt)
    expect_equal(result$alternative, alt)
    expect_gte(result$p.value, 0)
    expect_lte(result$p.value, 1)
  }
})
