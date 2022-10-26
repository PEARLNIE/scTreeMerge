
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
#' @param mcores the number of CPU cores.
#'
#' @return an object of class \code{Stree} which is a list.
#' @export
#' @importFrom methods setOldClass setClass new
#' @importFrom ape as.hclust.phylo
#' @importFrom phangorn acctran
#' @importFrom stats cutree
#' @importFrom phylogram as.dendrogram.phylo
#' @importFrom dendextend cutree_1k.dendrogram
#' @importFrom parallel detectCores makeCluster parLapply clusterExport stopCluster



minPkTree <- function(x, tree, min_k = NULL, max_k = NULL, mcores = NULL) {


  # Error checking.
  # if (!inherits(x, c("data.frame", "matrix")))
  #   stop("x must be object of class 'data.frame' or 'matrix'.")
  if (!inherits(x, "matrix"))

    x <- as.matrix(x)

  if (!class(tree) %in% c("phylo", "multiPhylo"))

    stop("tree must be object of class 'phylo' or 'multiPhylo'.")


  if (is.null(min_k))

    min_k <- 2

  if (is.null(max_k))

    max_k <- floor(sqrt(ncol(x)))

  if (is.null(mcores))

    mcores <- parallel::detectCores() - 1

  r <- length(tree)

  c <- max_k - min_k + 1


  mcl <- parallel::makeCluster(getOption("cl.cores", mcores))

  parallel::clusterExport(mcl, varlist = c("tree", "min_k", "max_k", "x", "calPvalue"), envir = environment())

  res <- parallel::parSapply(mcl, 1:r, function(u) {

    # print(names(tree)[u])

    p <- NULL

    tr1 <- tree[[u]]

    tr2 <- phylogram::as.dendrogram.phylo(tr1)

    for (i in min_k:max_k) {

      cl <- dendextend::cutree_1k.dendrogram(tr2, i, warn = FALSE)

      tmp <- ifelse(any(is.na(cl)), NA, calPvalue(x, cl, nsim = 1000))

      p <- c(p, tmp)
    }

    p

  })

  parallel::stopCluster(mcl)

  colnames(res) <- paste("tree", rep(1:r), sep = "")

  rownames(res) <- 2:max_k

  res <- t(res)

  res <- round(res, digits = 4)

  res_min <- min(res, na.rm = TRUE)

  ind <- findDim(x = res, object = res_min)

  row <- ind[1, 1]

  col <- ind[1, 2]

  return(new("Stree",
             raw = res,
             tree = tree[[row]],
             k = as.integer(colnames(res)[col]),
             pvalue = res_min))
}

# if (!is.null(data)) {
#   dat <- data
#   class(dat) <- "phyDat"
#   attr(dat, "weight") <- rep(1, nrow(data))
# } else {
#   stop("'data' should not be NULL, because 'tree[[", u, "]]' ", "has no branch length.")
# }
#
# # add banch length.
# tr1 <- phangorn::acctran(tr1, dat)


# data(GSE45719_268_count)
# processed_data <- getPPdata(GSE45719_268_count)
# d <- getDlist(x = t(processed_data), mtd = c("maximum", "euclidean", "manhattan",
#                                            "minkowski", "chebyshev", "sorensen",
#                                            "gower", "soergel", "kulczynski_d",
#                                            "canberra", "lorentzian", "intersection",
#                                            "non-intersection", "wavehedges", "czekanowski"))
# b <- getBasicPartitions(d, method = "all")
# m <- findBitrees(b$partition, "all")
# s <- minPkTree(x = processed_data, tree = m)
# class(s)
# str(s)







