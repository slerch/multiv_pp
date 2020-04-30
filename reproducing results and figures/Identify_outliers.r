# -------------------------------------------------------------------------- #

# Helping Function to identify and output outliers (modified boxplot.stats())
# The argument coef can be set to rescale the range of the "fences",
# Standard approach is do define a value an outlier if it is outside
# 1.5 times Interquartile range. For this standard setting, set coef=1.5

hf <- function (x, coef = 3.0, do.conf = TRUE, do.out = TRUE)
{
  if (coef < 0)
    stop("'coef' must not be negative")
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- stats::fivenum(x, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0)
    do.out <- FALSE
  else {
    out <- if (!is.na(iqr)) {
      x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef *
                                            iqr)
    }
    else !is.finite(x)
    if (any(out[nna], na.rm = TRUE))
      stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
  }
  conf <- if (do.conf)
    stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
  out <- if (do.out)  (out & nna) else numeric()
  out
}