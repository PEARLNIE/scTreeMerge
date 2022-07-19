


## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@



#' @title Basic Partitions Producted by Hierarchical Clustering
#' @description This function creates basic partitions with different combinations of distance measures and hierarchical methods.
#'
#' @param d a character \code{dist} or \code{list} which stores a number of objects of class \code{dist}.
#' @param method a character vector represents the chosen methods to compute. It must be one or combination of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"}, \code{"centroid"}. If users want to use 8 methods, \code{"all"} could be selected in this function.
#'
#' @return an object of class \code{list}:
#' \itemize{
#' \item partition : an object of class \code{multiPhylo} which includes numbers of basic partitions.
#' \item name : the combinations of specific distances and algorithms.
#' }
#' @export
#' @importFrom ape rtree as.phylo
#' @importFrom stats hclust
#' @examples
#' data(GSE45719_268_count)
#' processed_data <- getPPdata(GSE45719_268_count)
#'
#' d <- getDlist(x = t(processed_data), mtd = "euclidean")
#' class(d)
#' b <- getBasicPartitions(d, method = "complete")
#' class(b)
#' b1 <- getBasicPartitions(d, method = c("complete", "average"))
#' b2 <- getBasicPartitions(d, method = "all")
#'
#' d1 <- getDlist(x = t(processed_data), mtd = c("euclidean", "manhattan"))
#' b3 <- getBasicPartitions(d1, method = "all")

getBasicPartitions <- function(d, method = "complete") {
  # Error checking.
  if (!inherits(d, c("dist", "list"))) {
    stop("d must be object of class 'dist' or 'list'. ")
  } else {
    if (inherits(d, "dist")) {
      if (any(is.nan(d)))
        stop("d should not be used, because it has NaNs.")
      dname <- attr(d, "method")
      d <- list(dname = d)
    } else {
      if (inherits(d, "list")) {
        er <- NULL
        for (e in 1:length(d)) {
          if (any(is.nan(d[[e]])))
            er <- c(er, e)
          if (!is.null(er) & e == length(d))
            stop("d should not be used, because trees with the above indexes ", print(er), " have NaNs.")
        }
        dname <- names(d)
      }
    }
  }
  if (all(length(method) > 1 & "all" %in% method))
    stop("'all' covers all methods in this function. Once 'all' is chosen, there's no need to choose any of the 8 methods.")


  m <- c("ward.D",
         "ward.D2",
         "single",
         "complete",
         "average",
         "mcquitty",
         "median",
         "centroid") # 8

  res <- rtree(10) # A placeholder tree.
  ns <- c("a")

  u <- ifelse(inherits(d, "dist"), 1, length(d))
  unames <- dname
  if ("all" %in% method) {
    v <- m
  } else {
     v <- method
  }

  for (i in 1:u) {
    message("****************")
    message(paste("------", i, "-------", sep = " "))
    message("****************")

    for (g in v) {
      print(g)
      hc <- hclust(d[[i]], g)
      hc_tmp <- as.phylo(hc)
      res <- c(res, hc_tmp)
      ns <- c(ns, paste(unames[[i]], g, sep = " + "))
    }
  }
  names(res) <- ns
  # Deleting the 1st
  res <- res[2:length(res)]
  ns <- names(res)

  return(list("partition" = res, "name" = ns))
}








