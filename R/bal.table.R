#' Compute the balance table for a ps.cont object.
#'
#' @param x A `ps.cont` object
#' @param digits Number of digits to round to
#' @param ... Additional arguments.
#' @method bal.table ps.cont
#' @export

bal.table <-
  function(x, digits = 3, ...) {
    unw <- round(x$desc$unw$bal.tab$results, digits)
    wcor <- round(x$desc$AAC$bal.tab$results, digits)
    bal.tab <- cbind(unw, wcor)
    colnames(bal.tab) <- c("unw", "wcor")
    return(bal.tab)
  }
