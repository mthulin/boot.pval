## Test environments
* local Ubuntu 18.04 install, R 4.0.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * Found the following assignments to the global environment:
  File ‘boot.pval/R/boot_summary.R’:
  assign(".carEnv", car::.carEnv, envir = .GlobalEnv)

boot_summary() makes a call to car::Boot. The latter relies on the car environment being exported, which must be done manually (see https://cran.r-project.org/web/packages/car/vignettes/embedding.pdf). There is a check to see if there already is an object called .carEnv in the global environment, to avoid problems with overwriting an existing object.

## Downstream dependencies
None.
