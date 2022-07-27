

## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@


#' Tress from bipartition data
#'
#' This function is to find the optimal trees from a bipartition data.
#'
#' @param treelist an object of class \code{multiPhylo}.
#' @param format a \code{chanracter} vector represents the method of generating the ultimate tree. \code{"dtree"} represents the method based on the distance matrix. \code{"ptree"} represents the method based on maximum parisomy. \code{"ltree"} represents the method based on likelihood. It must be one or combinations of \code{"dtree"}, \code{"ptree"}, \code{"ltree"}. If users want to generate three formats of trees, \code{"all"} could be selected in this function.
#'
#' @return an object of class \code{multiPhylo} which is a list. The elements of this list are phylo trees.
#' @export
#' @importFrom parallelDist parDist
#' @importFrom phangorn as.phyDat optim.parsimony pratchet parsimony pml optim.pml upgma
#' @importFrom ape nj
#' @importFrom stats anova AIC
#' @examples
#' data(GSE45719_268_count)
#' processed_data <- getPPdata(GSE45719_268_count)
#' d <- getDlist(x = t(processed_data), mtd = c("maximum", "euclidean"))
#' b <- getBasicPartitions(d, method = "all")
#'
#' m <- findBitrees(b$partition, "dtree")
#' m1 <- findBitrees(b$partition, c("dtree", "ptree"))
#' m2 <- findBitrees(b$partition, "all")

findBitrees <- function(treelist, format = "all") {

  # Errors checking
  if (!inherits(treelist, "multiPhylo"))
    stop("treelist must be object of class 'multiPhylo'.")
  if (all(length(format) > 1 & "all" %in% format))
    stop("'all' covers all formats in this function. Once 'all' is chosen, there's no need to choose any of the 3 formats.")

  tree_format <- vector(mode = "character", length = 4)
  tree_format <- c("dtree", "ptree", "ltree", "all")
  if (!all(is.element(format, tree_format)))
    stop("Format '", format, "' is not implemented in this function.")

  # Extracting infos into a binary matrix.
  dat <- tree2Bidata(tree = treelist) # data.frame
  x <- dat
  class(x) <- "phyDat"
  attr(x, "weight") <- rep(1, nrow(dat))
  res <- rtree(length(treelist[[1]]$tip.label))
  res_name <- NULL

  d <- parallelDist::parDist(t(dat), method = "euclidean")
  # d <- getDlist(x = t(dat), mtd = "euclidean")

  if (any(format == "dtree")|"all" %in% format) {
    # ---- Distance-based ----
    dtree <- upgma(d)
    # str(dtree)
    attr(dtree, "distance") <- "euclidean"
    attr(dtree, "mtd") <- "upgma"
    res <- c(res, dtree)
    res_name <- c(res_name, "dtree")
  }

  tre_ini <- nj(d)
  if (any(format == "ptree")|"all" %in% format){
    # ---- Maximum parsimony ----
    tre_pars <- optim.parsimony(tree = tre_ini,
                                data = x,
                                method = "fitch",
                                cost = NULL,
                                trace = 1,
                                rearrangements = "SPR") # rearrangements: SPR or NNI
    # Or
    tre_p <- pratchet(data = x,
                      start = NULL,
                      method = "fitch",
                      maxit = 1000,
                      minit = 10,
                      k = 10,
                      trace = 1,
                      all = FALSE,
                      rearrangements = "SPR",
                      perturbation = "ratchet") # pratchet implements the parsimony ratchet, and is the preferred way to search for the best tree.


    # Returning the parsimony score of a tree using either the sankoff or the fitch algorithm.
    score <- NULL
    trees <- c(tre_ini, tre_pars, tre_p)
    for(i in trees) {
      tmp <- parsimony(tree = i,
                       data = x,
                       method = "fitch",
                       cost = NULL,
                       site = "pscore")
      score <- c(score, tmp)
    }
    mi <- which.min(score)
    ptree <- trees[[mi]]
    attr(ptree, "distance") <- "hamming"
    attr(ptree, "mtd") <- NULL
    res <- c(res, ptree)
    res_name <- c(res_name, "ptree")
  }

  if (any(format == "ltree")|"all" %in% format) {
    # ---- Maximum likelihood ----
    fit_ini <- pml(tree = tre_ini,
                   data = x,
                   bf = NULL,
                   Q = NULL,
                   inv = 0,
                   k = 1,
                   shape = 1,
                   rate = 1,
                   model = NULL,
                   site.rate = "gamma")
    # table(as.character(x))

    fit <- optim.pml(object = fit_ini,
                     optNni = TRUE, # optimized tree topology,
                     optBf = TRUE, # optimized base frequencies,
                     optQ = TRUE, # optimized rate matrix.
                     optGamma = TRUE)


    anova(fit_ini, fit) # ***：Extremely significant
    # Calculating AIC, the lower the value the better.
    if (AIC(fit_ini) > AIC(fit)){
      ltree <- fit
      attr(ltree, "distance") <- "hamming"
      attr(ltree, "mtd") <- NULL
    } else {
      ltree <- fit_ini
      attr(ltree, "distance") <- "hamming"
      attr(ltree, "mtd") <- NULL
    }
    res <- c(res, ltree$tree)
    res_name <- c(res_name, "ltree")
  }

  res <- res[-1]
  names(res) <- res_name
  return(res)
}
