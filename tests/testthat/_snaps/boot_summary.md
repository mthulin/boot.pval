# boot_summary errors when residual resampling is used with glm

    Code
      boot_summary(model_glm, R = 99, method = "residual")
    Condition
      Error in `boot_summary()`:
      ! Residual resampling is not recommended for GLM's (see http://www.modernstatisticswithr.com/regression.html#bootstrap-confidence-intervals-1). Please use method = "case" instead.

