model_lm <- lm(mpg ~ hp + vs, data = mtcars)

# car::Boot with residual resampling relies on .carEnv being available in the
# calling frame, which doesn't work inside test_that(). Use method = "case"
# (case resampling) for lm tests to avoid this limitation.

test_that("boot_summary returns a data frame with correct structure for lm", {
  result <- boot_summary(model_lm, R = 99, method = "case")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Estimate", "Lower.bound", "Upper.bound", "p.value"))
  expect_equal(nrow(result), length(coef(model_lm)))
})

test_that("boot_summary Estimate column matches model coefficients", {
  result <- boot_summary(model_lm, R = 99, method = "case")
  expect_equal(result$Estimate, unname(coef(model_lm)))
})

test_that("boot_summary CI bounds are finite numerics", {
  result <- boot_summary(model_lm, R = 99, method = "case")
  expect_type(result$Lower.bound, "double")
  expect_type(result$Upper.bound, "double")
  expect_true(all(is.finite(result$Lower.bound)))
  expect_true(all(is.finite(result$Upper.bound)))
})

test_that("boot_summary lower CI bound is below upper CI bound", {
  result <- boot_summary(model_lm, R = 99, method = "case")
  expect_true(all(result$Lower.bound < result$Upper.bound))
})

test_that("boot_summary with adjust.method adds Adjusted p-value column", {
  result <- boot_summary(model_lm, R = 99, method = "case", adjust.method = "holm")
  expect_true("Adjusted p-value" %in% names(result))
  expect_equal(nrow(result), length(coef(model_lm)))
})

test_that("boot_summary with coef = 'exp' returns positive estimates and CI bounds", {
  result <- suppressMessages(boot_summary(model_lm, R = 99, method = "case", coef = "exp"))
  expect_true(all(result$Estimate > 0))
  expect_true(all(result$Lower.bound > 0))
  expect_true(all(result$Upper.bound > 0))
})

test_that("boot_summary works with glm using case method", {
  model_glm <- glm(vs ~ mpg + hp, data = mtcars, family = binomial)
  result <- suppressWarnings(boot_summary(model_glm, R = 99))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Estimate", "Lower.bound", "Upper.bound", "p.value"))
  expect_equal(nrow(result), length(coef(model_glm)))
})

test_that("boot_summary errors when residual resampling is used with glm", {
  model_glm <- glm(vs ~ mpg + hp, data = mtcars, family = binomial)
  expect_snapshot(
    boot_summary(model_glm, R = 99, method = "residual"),
    error = TRUE
  )
})
