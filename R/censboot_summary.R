# Internal functions used for computing coefficients of bootstrap samples:
exp_reg_coef_cox <- function(data, formula) {
    m_boot <- survival::coxph(formula, data = data)
    return(exp(stats::coef(m_boot)))
}

reg_coef_cox <- function(data, formula) {
  m_boot <- survival::coxph(formula, data = data)
  return(stats::coef(m_boot))
}

exp_reg_coef_survreg <- function(data, formula, dist) {
  m_boot <- survival::survreg(formula, data = data, dist = dist)
  return(exp(stats::coef(m_boot)))
}

reg_coef_survreg <- function(data, formula, dist) {
  m_boot <- survival::survreg(formula, data = data, dist = dist)
  return(stats::coef(m_boot))
}

#' Summarising Survival Regression Models Using the Bootstrap
#'
#' Summaries for "coxph" and "survreg" objects, using the bootstrap for p-values and confidence intervals.
#'
#' @param model An object fitted using "lm", "glm", "nls", "lmer" or "glmer".
#' @param type A vector of character strings representing the type of interval to base the test on. The value should be one of "norm", "basic", "stud", "perc" (the default), and "bca". "stud" and "bca" are not available for "lmer" and "glmer" models.
#' @param sim The method used for bootstrapping. See \code{?boot::censboot} for details. Currently only "ordinary" (case resampling) is supported.
#' @param strata The strata used in the calls to \code{survfit.} It can be a vector or a matrix with 2 columns. If it is a vector then it is assumed to be the strata for the survival distribution, and the censoring distribution is assumed to be the same for all observations. If it is a matrix then the first column is the strata for the survival distribution and the second is the strata for the censoring distribution. When \code{sim = "ordinary"}, only one set of strata is used to stratify the observations. This is taken to be the first column of \code{strata} when it is a matrix.
#' @param coef A string specifying whether to use exponentiated coefficients in the summary table. Either "exp" (for exponentiated coefficients, i.e. hazard ratios in the case of a Cox PH model) or "raw" (for coefficients on their original scale). The default is "exp".
#' @param conf.level The confidence level for the confidence intervals. The default is 0.95.
#' @param R The number of bootstrap replicates. The default is 999.
#' @param pval_precision The desired precision for the p-value. The default is 1/R.
#' @param adjust.method Adjustment of p-values for multiple comparisons using \code{p.adjust}. The default is "none", in which case the p-values aren't adjusted. The other options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr"; see \code{?p.adjust} for details on these methods.
#' @param ... Additional arguments passed to \code{censboot}, such as \code{parallel} for parallel computations. See \code{?boot::censboot} for details.
#'
#' @return A data frame containing coefficient estimates, bootstrap confidence intervals, and bootstrap p-values.
#' @details p-values can be computed by inverting the corresponding confidence intervals, as described in Section 12.2 of Thulin (2021) and Section 3.12 in Hall (1992). This function computes p-values in this way from "boot" objects. The approach relies on the fact that:
#' - the p-value of the two-sided test for the parameter theta is the smallest alpha such that theta is not contained in the corresponding 1-alpha confidence interval,
#' - for a test of the parameter theta with significance level alpha, the set of values of theta that aren't rejected by the two-sided test (when used as the null hypothesis) is a 1-alpha confidence interval for theta.
#' @importFrom Rdpack reprompt
#' @references
#'  \insertRef{hall92}{boot.pval}
#'
#'  \insertRef{thulin21}{boot.pval}
#' @examples
#' library(survival)
#' # Weibull AFT model:
#' # Note that model = TRUE is required for use with censboot_summary:
#' model <- survreg(formula = Surv(time, status) ~ age + sex, data = lung,
#'                  dist = "weibull", model = TRUE)
#' censboot_summary(model, R = 99)
#' # (Values for R greater than 99 are recommended for most applications.)
#'
#' # Cox PH model:
#' model <- coxph(formula = Surv(time, status) ~ age + sex, data = lung,
#'                model = TRUE)
#' # Table with hazard ratios:
#' censboot_summary(model, R = 99)
#  # Table with original coefficients:
#' censboot_summary(model, coef = "raw", R = 99)
#' @export
censboot_summary <- function(model,
                 type = "perc",
                 sim = "ordinary",
                 strata = NULL,
                 coef = "exp",
                 conf.level = 0.95,
                 R = 999,
                 pval_precision = NULL,
                 adjust.method = "none",
                 ...)
{
  # Check arguments:
  if(!(class(model) %in% c("coxph", "survreg"))) { stop("The model must be fitted using either coxph or survreg (see ?censboot_summary).") }
  if(is.null(model$model)) { stop("The model must be fitted using model=TRUE (see ?censboot_summary).") }
  if(is.null(strata)) { strata <- matrix(1, nrow(model$y), 2) }
  cox <- ifelse(class(model) == "coxph", TRUE, FALSE)

  if(cox) {
    # Cox PH regression:
    boot_res <- boot::censboot(data = cbind(data.frame(as.matrix(model$model[,1])),
                          model$model[,2:ncol(model$model)]),
             statistic = switch(coef,
                                exp = exp_reg_coef_cox,
                                raw = reg_coef_cox),
             F.surv = survival::survfit(model),
             strata = strata,
             sim = sim,
             formula = model$formula,
             R = R,
             ...)
  } else {
  # AFT models:
  boot_res <- boot::censboot(data = cbind(data.frame(as.matrix(model$model[,1])),
                        model$model[,2:ncol(model$model)]),
           statistic = switch(coef,
                              exp = exp_reg_coef_survreg,
                              raw = reg_coef_survreg),
           F.surv = survival::survfit(model),
           strata = strata,
           sim = sim,
           formula = stats::formula(model$terms),
           dist = model$dist,
           R = R,
           ...)
  }

  # Create data frame to store results in:
  estimates <- switch(coef,
                      exp = exp(model$coefficients),
                      raw = model$coefficients)
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


