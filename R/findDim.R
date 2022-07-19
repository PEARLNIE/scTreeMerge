
## @.@.@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@.@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@

#' Return index with a specific value
#'
#' This function is to return the indexes with a user-defined value.
#'
#' @param x an object of class \code{matrix}.
#' @param object a user-defined value.
#'
#' @return an object of class \code{matrix} of which the 1st column represents the row index and the second represents the column index.
#' @export
#'
#' @examples
#' data(GSE45719_268_count)
#' res <- findDim(x = GSE45719_268_count, object = 516)
#' res


findDim <- function(x, object) {
  z <- NULL
  for (i in 1:nrow(x)) {
    e <- x[i,]
    tmp <- which(e == object)
    z <- c(z, tmp)
  }
  col <- as.integer(z)

  v <- NULL
  for(u in 1:ncol(x)) {
    e <- x[, u]
    tmp <- which(e == object)
    v <- c(v, tmp)
  }
  row <- as.integer(v)

  o <- matrix(c(row, col),ncol = 2, dimnames = list(NULL, c("row", "col")))
  o
}
