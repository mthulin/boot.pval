#' Convert Bootstrap Summary Tables to gt Objects
#'
#' Converts tables created using \code{boot_summary} and \code{censboot_summary} to nicely formatted \code{gt} tables.
#'
#' @param summary_table A table created using \code{boot_summary} or \code{censboot_summary}.
#' @param decimals The number of decimals to print for estimates and confidence intervals. The default is 3.
#' @param p_threshold p-values below this value will be printed as "<p_threshold", e.g. "<0.001". The default is 0.001.
#' @param conf The text at the top of the confidence interval column in the gt table. The default is "95 % CI".
#'
#' @return A gt table.
#' @examples
#' # Bootstrap summary of a linear model for mtcars:
#' model <- lm(mpg ~ hp + vs, data = mtcars)
#' boot_summary(model, R = 99) |> summary_to_gt()
#' @export
summary_to_gt <- function(summary_table,
                          decimals = 3,
                          p_threshold = 0.001,
                          conf = "95 % CI")
{
  summary_table |>
    gt::gt(rownames_to_stub = TRUE) |>
    gt::fmt_number(columns = c("Estimate", "Lower.bound", "Upper.bound"),
               decimals = decimals) |> # Display estimates and CI with 3 decimals
    gt::sub_small_vals(columns = "p.value",
                   threshold = p_threshold) |> # Show small p-values as "<0.001"
    gt::cols_merge(columns = c("Lower.bound", "Upper.bound"),
               pattern = "({1}, {2})") |> # Merge the CI bounds into a single column
    gt::cols_label(Lower.bound = conf,
               p.value = "p-value") -> bs_gt # Change the column names

  return(bs_gt)
}


#' Convert Bootstrap Summary Tables to flextable Objects
#'
#' Converts tables created using \code{boot_summary} and \code{censboot_summary} to nicely formatted \code{flextable} tables.
#'
#' @param summary_table A table created using \code{boot_summary} or \code{censboot_summary}.
#' @param decimals The number of decimals to print for estimates and confidence intervals. The default is 3.
#' @param conf The text at the top of the confidence interval column in the gt table. The default is "95 % CI".
#'
#' @return A flextable object.
#' @examples
#' # Bootstrap summary of a linear model for mtcars:
#' model <- lm(mpg ~ hp + vs, data = mtcars)
#' boot_summary(model, R = 99) |> summary_to_flextable()
#'
#' # Export to Word:
#' \dontrun{
#' boot_summary(model, R = 99) |>
#'    summary_to_flextable() |>
#'    flextable::save_as_docx(path = "my_table.docx")
#' }
#' @export
summary_to_flextable <- function(summary_table,
                          decimals = 3,
                          conf = "95 % CI")
{
  summary_table$Lower.bound <- round(summary_table$Lower.bound, 3)
  summary_table$Upper.bound <- round(summary_table$Upper.bound, 3)
  summary_table$Lower.bound <- paste0("(", summary_table$Lower.bound, ", ", summary_table$Upper.bound, ")")
  summary_table <- summary_table[,-3]
  summary_table$` ` <- row.names(summary_table)
  summary_table <- summary_table[,c(4, 1, 2, 3)]

  summary_table |>
    flextable::flextable() |>
    flextable::colformat_double(j=1:3, digits = decimals) |>
    flextable::set_header_labels(Lower.bound = conf,
                                 p.value = "p-value") |>
    flextable::vline(j = 1) |>
    flextable::set_table_properties(layout = "autofit") -> bs_fl

  return(bs_fl)
}

