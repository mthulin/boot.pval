t_stat <- function(data, i)
{
    data <- data[i,]
    sample1 <- data[data[, 2] == 1,]
    sample2 <- data[data[, 2] == 2,]
    mean0 <- mean(sample1[, 1]) - mean(sample2[, 1])
    variance <- stats::var(sample1[, 1])/nrow(sample1) + stats::var(sample2[, 1])/nrow(sample2)
    return(c(mean0, variance))
}

t_stat_1samp <- function(data, i)
{
    data <- data[i]
    mean0 <- mean(data)
    variance <- stats::var(data)/length(data)
    return(c(mean0, variance))
}

#' Bootstrap t-Test
#'
#' Performs one- and two-sample bootstrap t-tests and computes the corresponding bootstrap confidence interval.
#'
#' @param R The number of bootstrap replicates. The default is 9999.
#' @param type A vector of character strings representing the type of interval to base the test on. The value should be one of "norm", "basic", "bca", perc", and "stud" (the default).
#' @param ... Additional arguments passed to \code{boot}, such as \code{parallel} for parallel computations. See \code{?boot::boot} for details.
#' @inheritParams stats::t.test.default
#'
#' @details p-values can be computed by inverting the corresponding confidence intervals, as described in Section 14.2 of Thulin (2024) and Section 3.12 in Hall (1992). This function computes p-values for the t-test in this way. The approach relies on the fact that:
#' - the p-value of the two-sided test for the parameter theta is the smallest alpha such that theta is not contained in the corresponding 1-alpha confidence interval,
#' - for a test of the parameter theta with significance level alpha, the set of values of theta that aren't rejected by the two-sided test (when used as the null hypothesis) is a 1-alpha confidence interval for theta.
#' Consequently, the p-value will be consistent with the confidence interval, in the sense that the null hypothesis is rejected if and only if the null parameter values is not contained in the confidence interval.
#' @importFrom Rdpack reprompt
#' @return A list with class \code{"boot.htest"})
#' containing the following components:
#'   \item{statistic}{the value of the t-statistic.}
#'   \item{p.value}{the bootstrap p-value for the test.}
#'   \item{conf.int}{a bootstrap confidence interval for the mean appropriate to the
#'     specified alternative hypothesis.}
#'   \item{estimate}{the estimated mean or difference in means depending on
#'     whether it was a one-sample test or a two-sample test.}
#'   \item{null.value}{the specified hypothesized value of the mean or mean
#'     difference depending on whether it was a one-sample test or a
#'     two-sample test.}
#'   \item{stderr}{the standard error of the mean (difference), used as
#'     denominator in the t-statistic formula.}
#'   \item{alternative}{a character string describing the alternative
#'     hypothesis.}
#'   \item{method}{a character string indicating what type of t-test was
#'     performed.}
#'   \item{data.name}{a character string giving the name(s) of the data.}
#' @importFrom Rdpack reprompt
#' @references
#'  \insertRef{hall92}{boot.pval}
#'  \insertRef{thulin21}{boot.pval}
#' @examples
#' # Generate example data:
#' example_data <- data.frame(x = rnorm(40), y = rep(c(1,2), 20))
#'
#' # Two-sample (Welch) test:
#' boot_t_test(x ~ y, data = example_data, R = 999)
#'
#' # One-sample test:
#' boot_t_test(example_data$x, R = 999)
#'
#' # Paired test:
#' boot_t_test(example_data$x[example_data$y==1],
#'             example_data$x[example_data$y==2],
#'             paired = TRUE, R = 999)
#' @export
boot_t_test <- function(x, ...) UseMethod("boot_t_test")

#' @rdname boot_t_test
#' @export
boot_t_test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
           mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
           R = 9999, type = "stud",
           ...)
{
    alternative <- match.arg(alternative)

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
      stop("'mu' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")
    if( !is.null(y) ) {
      dname <- paste(deparse(substitute(x)),"and",
                     deparse(substitute(y)))
      if(paired)
        xok <- yok <- stats::complete.cases(x,y)
      else {
        yok <- !is.na(y)
        xok <- !is.na(x)
      }
      y <- y[yok]
    }
    else {
      dname <- deparse(substitute(x))
      if (paired) stop("'y' is missing for paired test")
      xok <- !is.na(x)
      yok <- NULL
    }
    x <- x[xok]
    if (paired) {
      x <- x-y
      y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- stats::var(x)
    if(is.null(y)) {
      if(nx < 2) stop("not enough 'x' observations")
      df <- nx-1
      stderr <- sqrt(vx/nx)
      if(stderr < 10 *.Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      tstat <- (mx-mu)/stderr
      method <- if(paired) "Paired Bootstrap t-test" else "One Sample Bootstrap t-test"
      estimate <-
        stats::setNames(mx, if(paired)"mean of the differences" else "mean of x")
      boot_res <- boot::boot(x, t_stat_1samp, R = R, stype = "i", ...)
    } else {
      ny <- length(y)
      if(nx < 1 || (!var.equal && nx < 2))
        stop("not enough 'x' observations")
      if(ny < 1 || (!var.equal && ny < 2))
        stop("not enough 'y' observations")
      if(var.equal && nx+ny < 3) stop("not enough observations")
      my <- mean(y)
      vy <- stats::var(y)
      method <- paste(if(!var.equal)"Welch", "Two Sample Bootstrap t-test")
      estimate <- c(mx,my)
      names(estimate) <- c("mean of x","mean of y")
      if(var.equal) {
        stop("Do not use the equal-variances t-test. Use the Welch t-test instead.")
      } else {
        stderrx <- sqrt(vx/nx)
        stderry <- sqrt(vy/ny)
        stderr <- sqrt(stderrx^2 + stderry^2)
        df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
      }
      if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
        stop("data are essentially constant")

      data <- data.frame(x = c(x, y), group = rep(c(1, 2), c(length(x), length(y))))
      boot_res <- boot::boot(data, t_stat, R = R, stype = "i", strata = data[, 2], ...)
      tstat <- (mx - my - mu)/stderr
    }
    if (alternative == "less") {
      pval <- NA
      cint <- NA
    }
    else if (alternative == "greater") {
      pval <- NA
      cint <- NA
    }
    else {
      pval <- boot.pval(boot_res, type = type, theta_null = mu)
      ci <- boot::boot.ci(boot_res, type = type, conf = conf.level)
      cint <- switch(type,
                     norm = ci$normal[,2:3],
                     basic = ci$basic[,4:5],
                     stud = ci$student[,4:5],
                     perc = ci$percent[,4:5],
                     bca = ci$bca[,4:5])
    }

    method <- paste(method, switch(type,
                                   norm = "(normal)",
                                   basic = "(basic)",
                                   stud = "(studentized)",
                                   perc = "(percentile)",
                                   bca = "BCa)"))

    names(tstat) <- "t"
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, #parameter = df, # uncomment to make gtsummary::tbl_regression() work with output
                 p.value = pval,
                 conf.int = cint, estimate = estimate, null.value = mu,
                 alternative = alternative,
                 method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}

#' @rdname boot_t_test
#' @export
boot_t_test.formula <- function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
      stop("grouping factor must have exactly 2 levels")
    DATA <- stats::setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("boot_t_test", c(DATA, list(...)))
    y$data.name <- DNAME
    if(length(y$estimate) == 2L)
      names(y$estimate) <- paste("mean in group", levels(g))
    y
}
