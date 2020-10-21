#' Displays a useful description of a `ps.cont` object.
#'
#' @param object A `ps.cont` object
#' @param ... Additional arguments.
#'
#' @method summary ps.cont
#' @export

summary.ps.cont <- function(object,...){
  summary.names <- c("n", "ess", "max.wcor", "mean.wcor", "rms.wcor", "n.trees")

  summary.tab.list <- lapply(object$desc, function(x) {
    do.call("c", x[summary.names])
  })

  summary.tab <- do.call("rbind", summary.tab.list)

  colnames(summary.tab)[colnames(summary.tab) == "n.trees"] <- "iter"

  #class(summary.tab) <- "summary.ps"
  return(summary.tab)
}
