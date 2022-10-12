


## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@

#' @title Extracting tree bipartitions
#' @description  This function is to extracting the information of bipartitions found in a series of trees.
#'
#' @param tree an object of class \code{"multiPhylo"} or \code{"phylo"}. Note that trees in \code{"multiPhylo"} should have the same tip labels.
#'
#' @return an object of class \code{data.frame}. Its columns represent samples and rows represent different bipartitions.
#' @export
#' @importFrom ape .uncompressTipLabel di2multi prop.part



tree2Bidata <- function(tree) {

  # some minor error checking
  if (!inherits(tree, c("multiPhylo", "phylo")))
    stop("tree must be object of class 'multiPhylo' or 'phylo'.")
  if (class(tree) == "phylo")
    tree <- c(tree)

  tree <- .uncompressTipLabel(tree)
  tree <- di2multi(tree)

  labels <- lapply(tree, function(x) sort(x$tip.label)) # Sorting from smallest to largest.
  ulabels <- unique(labels)
  lul <- length(ulabels)
  # compute matrix representation phylogenies
  X <- vector("list", lul) # list of bipartitions
  characters <- 0 # number of characters
  weights <- NULL
  species <- tree[[1]]$tip.label


  for (i in 1:lul) {
    pos <- match(labels, ulabels[i])
    ind <- which(!is.na(pos))
    temp <- prop.part(tree[ind]) # find all bipartitions
    # create matrix representation of tree[[i]] in X[[i]]
    TMP <- matrix(0L, nrow = length(temp) - 1,
                  ncol = length(tree[[ind[1]]]$tip.label))
    for (j in seq_len(nrow(TMP)))  TMP[j, c(temp[[j + 1]])] <- 1L
    colnames(TMP) <- attr(temp, "labels") # label rows
    X[[i]] <- TMP
    species <- union(species, tree[[ind[1]]]$tip.label) # accumulate labels
    characters <- characters + nrow(TMP) # count characters
    weights <- c(weights, attr(temp, "number")[-1])
  }

  data <- matrix(data = NA, nrow = characters, ncol = length(species),
                 dimnames = list(NULL, species))
  j <- 1
  for (i in seq_along(X)) {
    # copy each of X into supermatrix data
    data[c(j:((j - 1) + nrow(X[[i]]))), colnames(X[[i]])] <- X[[i]]
    # [1:nrow(X[[i]]),1:ncol(X[[i]])]
    j <- j + nrow(X[[i]])
  }
  data <- as.data.frame(data)
  # compute contrast matrix
  contrast <- matrix(data = c(1, 0, 0, 1, 1, 1), nrow = 3, ncol = 2,
                     dimnames = list(NULL, c("0", "1")), byrow = TRUE)

  # attr(data, "row.names") <- NULL
  # class(data) <- "phyDat"
  attr(data, "weight") <- weights
  attr(data, "nr") <- length(weights)
  attr(data, "nc") <- 2L
  attr(data, "levels") <- c("0", "1")
  attr(data, "allLevels") <- c("0", "1")
  attr(data, "type") <- "USER"
  attr(data, "contrast") <- contrast
  return(data)
}



# data(GSE45719_268_count)
# processed_data <- getPPdata(GSE45719_268_count)
# d <- getDlist(x = t(processed_data), mtd = c("maximum", "euclidean", "manhattan",
#                                            "minkowski", "chebyshev", "sorensen",
#                                            "gower", "soergel", "kulczynski_d",
#                                            "canberra", "lorentzian", "intersection",
#                                            "non-intersection", "wavehedges", "czekanowski"))
# b <- getBasicPartitions(d, method = "all")
# tree <- b$partition
# bidat <- tree2Bidata(tree)
# str(bidat)

