# Convert Bootstrap Summary Tables to gt Objects

Converts tables created using `boot_summary` and `censboot_summary` to
nicely formatted `gt` tables.

## Usage

``` r
summary_to_gt(summary_table, decimals = 3, conf = "95 % CI")
```

## Arguments

- summary_table:

  A table created using `boot_summary` or `censboot_summary`.

- decimals:

  The number of decimals to print for estimates and confidence
  intervals. The default is 3.

- conf:

  The text at the top of the confidence interval column in the gt table.
  The default is "95 % CI".

## Value

A gt table.

## Examples

``` r
# Bootstrap summary of a linear model for mtcars:
model <- lm(mpg ~ hp + vs, data = mtcars)
boot_summary(model, R = 99) |> summary_to_gt()


  
```
