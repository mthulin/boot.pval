#' Summarising Regression Models Using the Bootstrap
#'
#' Summaries for "lm", "glm", and "nls" objects, using the bootstrap for p-values and confidence intervals.
#'
#' @param model An object of class "lm", "glm", or "nls".
#' @param type A vector of character strings representing the type of interval to base the test on. The value should be one of "norm", "basic", "stud", "perc" (the default), and "bca".
#' @param method Either "residual" (for resampling of scaled and centred residuals, the default) or "case" (for case resampling).
#' @param conf.level The confidence level for the confidence intervals. The default is 0.95.
#' @param R The number of bootstrap replicates. The default is 999.
#' @param pval_precision The desired precision for the p-value. The default is 1/R.
#' @param adjust.method Adjustment of p-values for multiple comparisons using \code{p.adjust}. The default is "none", in which case the p-value aren't adjusted. The other options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr"; see \code{?p.adjust} for details on these methods.
#' @param ... Additional arguments passed to \code{boot}.
#'
#' @return A data frame containing coefficient estimates, bootstrap confidence intervals, and bootstrap p-values.
#' @importFrom Rdpack reprompt
#' @examples
#' # Bootstrap summary of a linear model for mtcars:
#' model <- lm(mpg ~ hp + wt, data = mtcars)
#' boot_summary(model)
#'
#' # Adjust p-values for multiplicity using Holm's method:
#' boot_summary(model, adjust.method = "holm")
#' @export
boot_summary <- function(model,
                      type = "perc",
                      method = "residual",
                      conf.level = 0.95,
                      R = 999,
                      pval_precision = NULL,
                      adjust.method = "none",
                      ...)
{
  # Bootstrap the regression model:
  library(car) # Required for Boot to work inside of function
  .carEnv <- car:::.carEnv # Required for Boot to work inside of function
  boot_res <- car::Boot(model,
                        method = method,
                        R = R,
                        ...)

  # Create data frame to store results in:
  p <- length(model$coefficients)
  results <- data.frame(Estimate = model$coefficients,
                       Lower = rep(NA, p),
                       Upper = rep(NA, p),
                       p.value = rep(NA, p))

  # Compute confidence intervals:
  for(i in 1:p)
  {
      ci <- boot::boot.ci(boot_res, conf = conf.level, type = type, ...)
      results[i, 2:3] <- switch(type,
                       norm = ci$normal[,2:3],
                       basic = ci$basic[,4:5],
                       stud = ci$student[,4:5],
                       perc = ci$percent[,4:5],
                       bca = ci$bca[,4:5])
  }

  # Compute p-values:
  for(i in 1:p)
  {
    results[i, 4] <- boot.pval(boot_res,
                              type = type,
                              theta_null = 0,
                              pval_precision = pval_precision,
                              ...)
  }
  results[, 4] <- stats::p.adjust(results[, 4], method = adjust.method)

  return(results)
}


