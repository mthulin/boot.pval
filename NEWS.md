# boot.pval version 0.7
Added support for one-sided tests. Added the boot_t_test and boot_median_test functions for carrying out bootstrap tests of location. Updated documentation for boot_summary and censboot_summary to not mention studentized intervals (which aren't supported by upstream packages). Added a new vignette.

# boot.pval version 0.6
Added automatic handling of missing values, so that these don't have to be removed manually from the data prior to using boot_summary(). Added support for BCa intervals again, and improved performance for these intervals. Improved presentation of p-values in regression models. Added a vignette.

# boot.pval version 0.5
Added support for AFT models fitted using rms::psm. Added functions for creating publication-ready summary tables using the gt and flextable packages. Removed support for BCa intervals for some regression models, as these are not longer supported by upstream packages.

# boot.pval version 0.4.1
Fixed a bug that caused boot_summary to throw an error when used with GLM's in R version >= 4.2.

# boot.pval version 0.4
Two bugs in censboot_summary() have been fixed: a) a bug which caused incorrect p-values for exponentiated coefficients, and b) a bug that cause the function to fail if the variables in the Surv object weren't named time and status. An options for creating a table for exponentiated coefficients using boot_summary() has also been added (useful e.g. for logistic regression models).

# boot.pval version 0.3
boot_summary() now also works for mixed linear models fitted using the lmerTest package. In previous releases, such models weren't correctly identified by boot_summary().

# boot.pval version 0.2
Changes have been made in the documentation, and additional error messages have been added.

# boot.pval version 0.1
This is the first release of the package.
