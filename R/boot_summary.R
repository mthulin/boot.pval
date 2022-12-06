#' Summarising Regression Models Using the Bootstrap
#'
#' Summaries for regression models, including "lm", "glm", "glm.nb", nls", "rlm", and "merMod" ("lmer", "glmer") objects, using the bootstrap for p-values and confidence intervals.
#'
#' @param model An object fitted using e.g. "lm", "glm", "glm.nb", "nls", "rlm", "lmer", or "glmer".
#' @param type A vector of character strings representing the type of interval to base the test on. The value should be one of "norm", "basic", "stud", "perc" (the default), and "bca". "stud" and "bca" are not available for "lmer" and "glmer" models.
#' @param method The method used for bootstrapping. For "lm" and "nls" objects use either "residual" (for resampling of scaled and centred residuals, the default) or "case" (for case resampling). For "glm" objects, use "case" (the default). For "merMod" objects (mixed models) use either "parametric" (the default) or "semiparametric".
#' @param conf.level The confidence level for the confidence intervals. The default is 0.95.
#' @param R The number of bootstrap replicates. The default is 999.
#' @param coef A string specifying whether to use exponentiated coefficients in the summary table. Either "exp" (for exponentiated coefficients, i.e. odds ratios in the case of a logistic regression model) or "raw" (for coefficients on their original scale). The default is "raw", which is recommended for linear models.
#' @param pval_precision The desired precision for the p-value. The default is 1/R.
#' @param adjust.method Adjustment of p-values for multiple comparisons using \code{p.adjust}. The default is "none", in which case the p-values aren't adjusted. The other options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr"; see \code{?p.adjust} for details on these methods.
#' @param ... Additional arguments passed to \code{Boot} or \code{bootMer}, such as \code{parallel} for parallel computations. See \code{?car::Boot} and \code{?lme4::bootMer} for details.
#'
#' @return A data frame containing coefficient estimates, bootstrap confidence intervals, and bootstrap p-values.
#' @details p-values can be computed by inverting the corresponding confidence intervals, as described in Section 12.2 of Thulin (2021) and Section 3.12 in Hall (1992). This function computes p-values for coefficients of regression models in this way. The approach relies on the fact that:
#' - the p-value of the two-sided test for the parameter theta is the smallest alpha such that theta is not contained in the corresponding 1-alpha confidence interval,
#' - for a test of the parameter theta with significance level alpha, the set of values of theta that aren't rejected by the two-sided test (when used as the null hypothesis) is a 1-alpha confidence interval for theta.
#'
#' The function can be used with "lm", "glm", "glm.nb", "nls", "rlm", and "merMod" ("lmer", "glmer") objects. In addition, it should work for any regression model such that: \code{residuals(object, type="pearson")} returns Pearson residuals; \code{fitted(object)} returns fitted values; \code{hatvalues(object)} returns the leverages, or perhaps the value 1 which will effectively ignore setting the hatvalues. In addition, the \code{data} argument should contain no missing values among the columns actually used in fitting the model.
#' @importFrom Rdpack reprompt
#' @references
#'  \insertRef{hall92}{boot.pval}
#'
#'  \insertRef{thulin21}{boot.pval}
#' @examples
#' # Bootstrap summary of a linear model for mtcars:
#' model <- lm(mpg ~ hp + vs, data = mtcars)
#' boot_summary(model, R = 99)
#' # (Values for R greater than 99 are recommended for most applications.)
#'
#' # Adjust p-values for multiplicity using Holm's method:
#' boot_summary(model, R = 99, adjust.method = "holm")
#' @export
boot_summary <- function(model,
                      type = "perc",
                      method = NULL,
                      conf.level = 0.95,
                      R = 999,
                      coef = "raw",
                      pval_precision = NULL,
                      adjust.method = "none",
                      ...)
{
  # Bootstrap the regression model:
  # (Different functions are used depending on the type of model.)
  if(class(model)[1] %in% c("lmerMod", "glmerMod", "lmerModLmerTest"))
  {
    # Use lme4::bootMer for mixed models:
    if(is.null(method)) { method <- "parametric" }
    boot_res <- lme4::bootMer(model,
                              FUN = lme4::fixef,
                              type = method,
                              nsim = R,
                              ...)
  } else {
    # Use car::Boot for other objects, including lm, glm, and nls objects:
    if(is.null(method)) { if(class(model)[1] == "glm") { method <- "case" } else { method <- "residual" } }

    # Throw an error if the user attempts to use residual resampling with a GLM:
    if(class(model)[1] == "glm" & method == "residual") { stop("Residual resampling is not recommended for GLM's (see http://www.modernstatisticswithr.com/regression.html#bootstrap-confidence-intervals-1). Please use method = \"case\" instead.") }

    # For Boot to work inside this function, we must export the car environment
    # to the global environment:
    # (see https://cran.r-project.org/web/packages/car/vignettes/embedding.pdf).
    assign(".carEnv", car::.carEnv, envir = parent.frame())
    boot_res <- car::Boot(model,
                          method = method,
                          R = R,
                          ...)
  }

  # Create data frame to store results in:
  if(class(model)[1] %in% c("lmerMod", "glmerMod", "lmerModLmerTest")) {
    estimates <- lme4::fixef(model)
    } else estimates <- model$coefficients
  p <- length(estimates)
  results <- data.frame(Estimate = estimates,
                       `Lower bound` = rep(NA, p),
                       `Upper bound` = rep(NA, p),
                       `p-value` = rep(NA, p))

  # Compute confidence intervals:
  for(i in 1:p)
  {
      ci <- boot::boot.ci(boot_res, conf = conf.level, type = type, index = i, ...)
      results[i, 2:3] <- switch(type,
                       norm = ci$normal[,2:3],
                       basic = ci$basic[,4:5],
                       stud = ci$student[,4:5],
                       perc = ci$percent[,4:5],
                       bca = ci$bca[,4:5])
  }

  # Exponentiate coefficients if requested:
  if(coef == "exp") {
    cat("Using exponentiated coefficients.\n")
    results[,1:3] <- exp(results[,1:3])
  }

  # Compute p-values:
  for(i in 1:p)
  {
    results[i, 4] <- boot.pval(boot_res,
                              type = type,
                              theta_null = 0,
                              pval_precision = pval_precision,
                              index = i,
                              ...)
  }
  if(adjust.method != "none") {
    results$`Adjusted p-value` <- round(stats::p.adjust(results[, 4], method = adjust.method), round(log10(R)))
  }
  # Round p-values:
  results[,4] <- round(results[,4], round(log10(R)))


  return(results)
}


