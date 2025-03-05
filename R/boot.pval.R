#' Compute Bootstrap p-values
#'
#' Compute bootstrap p-values through confidence interval inversion, as described in Hall (1992) and Thulin (2024).
#'
#' @param boot_res An object of class "boot" containing the output of a bootstrap calculation.
#' @param type A vector of character strings representing the type of interval to base the test on. The value should be one of "norm", "basic", "stud", "perc" (the default), and "bca".
#' @param theta_null The value of the parameter under the null hypothesis.
#' @param pval_precision The desired precision for the p-value. The default is 1/R,  where R is the number of bootstrap samples in \code{boot_res}.
#' @param alternative A character string specifying the alternative hypothesis. Must be one of "two.sided" (default), "greater", or "less".
#' @param ... Additional arguments passed to \code{boot.ci}.
#'
#' @return A bootstrap p-value.
#' @details p-values can be computed by inverting the corresponding confidence intervals, as described in Section 14.2 of Thulin (2024) and Section 3.12 in Hall (1992). This function computes p-values in this way from "boot" objects. The approach relies on the fact that:
#' - the p-value of the two-sided test for the parameter theta is the smallest alpha such that theta is not contained in the corresponding 1-alpha confidence interval,
#' - for a test of the parameter theta with significance level alpha, the set of values of theta that aren't rejected by the two-sided test (when used as the null hypothesis) is a 1-alpha confidence interval for theta.
#' @importFrom Rdpack reprompt
#' @references
#'  \insertRef{hall92}{boot.pval}
#'  \insertRef{thulin21}{boot.pval}
#' @seealso [boot_t_test()] for bootstrap t-tests, [boot_median_test()] for bootstrap tests for medians, [boot_summary()] for bootstrap tests for coefficients of regression models.
#' @examples
#' # Hypothesis test for the city data
#' # H0: ratio = 1.4
#' library(boot)
#' ratio <- function(d, w) sum(d$x * w)/sum(d$u * w)
#' city.boot <- boot(city, ratio, R = 99, stype = "w", sim = "ordinary")
#' boot.pval(city.boot, theta_null = 1.4)
#'
#' # Studentized test for the two sample difference of means problem
#' # using the final two series of the gravity data.
#' diff.means <- function(d, f)
#' {
#'   n <- nrow(d)
#'   gp1 <- 1:table(as.numeric(d$series))[1]
#'   m1 <- sum(d[gp1,1] * f[gp1])/sum(f[gp1])
#'   m2 <- sum(d[-gp1,1] * f[-gp1])/sum(f[-gp1])
#'   ss1 <- sum(d[gp1,1]^2 * f[gp1]) - (m1 *  m1 * sum(f[gp1]))
#'   ss2 <- sum(d[-gp1,1]^2 * f[-gp1]) - (m2 *  m2 * sum(f[-gp1]))
#'   c(m1 - m2, (ss1 + ss2)/(sum(f) - 2))
#' }
#' grav1 <- gravity[as.numeric(gravity[,2]) >= 7, ]
#' grav1.boot <- boot(grav1, diff.means, R = 99, stype = "f",
#'                    strata = grav1[ ,2])
#' boot.pval(grav1.boot, type = "stud", theta_null = 0)
#' @export
boot.pval <- function(boot_res,
                      type = "perc",
                      theta_null = 0,
                      pval_precision = NULL,
                      alternative = "two.sided",
                      ...)
{
    if(is.null(pval_precision)) { pval_precision = 1/boot_res$R }

    # Create a sequence of alphas:
    alpha_seq <- seq(pval_precision, 1-pval_precision, pval_precision)

    # Compute the 1-alpha confidence intervals, and extract
    # their bounds:
    ci <- suppressWarnings(boot::boot.ci(boot_res,
            conf = 1-alpha_seq,
            type = type,
            ...))

    bounds <- switch(type,
           norm = ci$normal[,2:3],
           basic = ci$basic[,4:5],
           stud = ci$student[,4:5],
           perc = ci$percent[,4:5],
           bca = ci$bca[,4:5])

    # Find the smallest alpha such that theta_null is not contained in the 1-alpha
    # confidence interval:
    if(!(alternative %in% c("two.sided", "less", "greater"))) {
      stop('alternative must be either "two.sided", "less" or "greater"')
    }
    if(alternative == "two.sided") {
      alpha <- alpha_seq[which.min(theta_null >= bounds[,1] & theta_null <= bounds[,2])]
    }
    if(alternative == "greater") {
      if(sum(theta_null < bounds[,1]) == 0)
      { alpha <- 1 } else
      {
        alpha <- alpha_seq[which.min(theta_null >= bounds[,1])]/2
      }

    }
    if(alternative == "less") {
      if(sum(theta_null > bounds[,2]) == 0)
      { alpha <- 1 } else
      {
        alpha <- alpha_seq[which.min(theta_null <= bounds[,2])]/2
      }
    }

  # Return the p-value:
  return(alpha)
}
