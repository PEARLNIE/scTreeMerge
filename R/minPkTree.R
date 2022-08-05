
## @.@.@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@.@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@





#' Selecting one tree with a specific k from a set of phylo trees
#'
#' This function is to find the tree with the optimal k identified by the averaged P value. The averaged P value depends on how different the k clusters are. The smaller the averaged P value is, the more separated the clusters are.

#' @param x an object of class \code{'matrix'}. Each column corresponds to a sample and each row to a variable.
#' @param tree an object of class \code{phylo} or \code{multiPhylo}.
#' @param min_k integer value. minimum cluster number to evaluate.
#' @param max_k integer value. maximum cluster number to evaluate.
#' @param iteration integer value. number of iterations. Default is 10.
#'
#' @return an object of class \code{Stree} which is a list.
#' @export
#' @importFrom methods setOldClass setClass new
#' @importFrom ape as.hclust.phylo
#' @importFrom stats cutree
#' @examples
#' data(GSE45719_268_count)
#' processed_data <- getPPdata(GSE45719_268_count)
#' d <- getDlist(x = t(processed_data), mtd = c("maximum", "euclidean", "manhattan",
#'                                            "minkowski", "chebyshev", "sorensen",
#'                                            "gower", "soergel", "kulczynski_d",
#'                                            "canberra", "lorentzian", "intersection",
#'                                            "non-intersection", "wavehedges", "czekanowski"))
#' b <- getBasicPartitions(d, method = "all")
#' m <- findBitrees(b$partition, "all")
#' s <- minPkTree(x = processed_data, tree = m, min_k = 2, max_k = 3, iteration = 10)
#' class(s)
#' str(s)


minPkTree <- function(x, tree, min_k = NULL, max_k = NULL, iteration = 10) {


  # Error checking.
  # if (!inherits(x, c("data.frame", "matrix")))
  #   stop("x must be object of class 'data.frame' or 'matrix'.")
  if (!inherits(x, "matrix"))
    stop("x must be object of class or 'matrix'.")
  if (!class(tree) %in% c("phylo", "multiPhylo"))
    stop("tree must be object of class 'phylo' or 'multiPhylo'.")


  if (is.null(min_k))
    min_k <- 2
  if (is.null(max_k))
    max_k <- floor(sqrt(ncol(x)))

  r <- length(tree)
  c <- max_k - min_k + 1
  dname <- list(paste("iter", rep(1:iteration), sep = ""), 2:max_k)
  res_tmp <- matrix(NA, nrow = iteration, ncol = c, dimnames = dname)
  res <- matrix(NA, nrow = r, ncol = c, dimnames = list(paste("tree", rep(1:r), sep = ""), 2:max_k))

  for(u in 1:r) {
    tr <- as.hclust.phylo(tree[[u]])
    p <- NULL
    k <- NULL
    for(i in min_k:max_k) {
      cl <- cutree(tr, i)
      tmp <- calPvalue(x = x, clusters = cl, nsim = 1000)
      p <- c(p, tmp)
    }
    res[u, ] <- p
  }
  # t <- res
  res <- round(res, digits = 4)
  res_min <- min(res)
  ind <- findDim(x = res, object = res_min)
  col <- ind[1, 1]
  row <- ind[1, 2]

  return(new("Stree", raw = res, tree = tree[[row]],
             k = as.integer(colnames(res)[row]), pvalue = res_min))
}




