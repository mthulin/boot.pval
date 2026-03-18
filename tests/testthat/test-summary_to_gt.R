model <- lm(mpg ~ hp + vs, data = mtcars)
bs <- boot_summary(model, R = 99)

test_that("summary_to_gt returns a gt_tbl object", {
  result <- summary_to_gt(bs)
  expect_s3_class(result, "gt_tbl")
})

test_that("summary_to_gt works with custom decimals and conf label", {
  result <- summary_to_gt(bs, decimals = 2, conf = "90 % CI")
  expect_s3_class(result, "gt_tbl")
})

test_that("summary_to_flextable returns a flextable object", {
  result <- summary_to_flextable(bs)
  expect_s3_class(result, "flextable")
})

test_that("summary_to_flextable works with custom decimals and conf label", {
  result <- summary_to_flextable(bs, decimals = 2, conf = "90 % CI")
  expect_s3_class(result, "flextable")
})
