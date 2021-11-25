# boot.pval version 0.4
Two bugs in censboot_summary() have been fixed: a) a bug which caused incorrect p-values for exponentiated coefficients, and b) a bug that cause the function to fail if the variables in the Surv object weren't named time and status. An options for creating a table for exponentiated coefficients using boot_summary() has also been added (useful e.g. for logistic regression models).

# boot.pval version 0.3
boot_summary() now also works for mixed linear models fitted using the lmerTest package. In previous releases, such models weren't correctly identified by boot_summary().

# boot.pval version 0.2
Changes have been made in the documentation, and additional error messages have been added.

# boot.pval version 0.1
This is the first release of the package.
