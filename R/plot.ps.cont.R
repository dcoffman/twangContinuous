#' Plot the `ps.cont` object.
#'
#' @param x ps.cont object
#' @param plots An indicator of which type of plot is desired. The options are
#'   * `"optimize"` A plot of the balance criteria as a function of the GBM
#'     iteration.
#'   * `"es"` Plots of the standardized effect size of the pre-treatment
#'     variables before and after weighting
#' @param subset Used to restrict which of the `stop.method`s will be used
#'   in the figure. For example `subset = c(1,3)` would indicate that the
#'   first and third `stop.method`s (in alphabetical order of those specified
#'   in the original call to the ps.cont function) should be included in the
#'   figure.
#' @param ... Additional arguments.
#'
#' @method plot ps.cont
#' @export

plot.ps.cont <- function(x, plots="optimize", subset = NULL, ...) {
  pt1 <- diag.plot(x, plots, subset, ...)
  return(pt1)
}
