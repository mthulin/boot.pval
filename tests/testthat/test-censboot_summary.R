library(survival)

lung_clean <- lung[complete.cases(lung[, c("time", "status", "age", "sex")]), ]
cox_model <- coxph(
  Surv(time, status) ~ age + sex,
  data = lung_clean,
  model = TRUE
)

test_that("censboot_summary returns a data frame with correct structure for coxph", {
  result <- suppressMessages(censboot_summary(cox_model, R = 99))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Estimate", "Lower.bound", "Upper.bound", "p.value"))
  expect_equal(nrow(result), length(coef(cox_model)))
})

test_that("censboot_summary CI bounds are finite numerics", {
  result <- suppressMessages(censboot_summary(cox_model, R = 99))
  expect_type(result$Lower.bound, "double")
  expect_type(result$Upper.bound, "double")
  expect_true(all(is.finite(result$Lower.bound)))
  expect_true(all(is.finite(result$Upper.bound)))
})

test_that("censboot_summary lower CI bound is below upper CI bound", {
  result <- suppressMessages(censboot_summary(cox_model, R = 99))
  expect_true(all(result$Lower.bound < result$Upper.bound))
})

test_that("censboot_summary with coef = 'raw' returns correct structure", {
  result <- suppressMessages(censboot_summary(cox_model, R = 99, coef = "raw"))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Estimate", "Lower.bound", "Upper.bound", "p.value"))
})

test_that("censboot_summary with adjust.method adds Adjusted p-value column", {
  result <- suppressMessages(
    censboot_summary(cox_model, R = 99, adjust.method = "holm")
  )
  expect_true("Adjusted p-value" %in% names(result))
})

test_that("censboot_summary errors without model = TRUE", {
  cox_no_model <- coxph(Surv(time, status) ~ age + sex, data = lung_clean)
  expect_snapshot(censboot_summary(cox_no_model, R = 99), error = TRUE)
})

test_that("censboot_summary works with survreg model", {
  survreg_model <- survreg(
    Surv(time, status) ~ age + sex,
    data = lung_clean,
    dist = "weibull",
    model = TRUE
  )
  result <- suppressMessages(censboot_summary(survreg_model, R = 99))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(coef(survreg_model)))
})
