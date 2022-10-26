

## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@



#' Find the popular trees
#'
#' This function is to find the popular trees with the frequency of tree types.
#'
#' @param treelist an object of class \code{multiPhylo}.
#'
#' @return an object of class \code{list}.
#' @export
#' @import igraph
#' @importFrom stringr str_split
#'


findPOPtrees <- function(treelist) {

  namelist <- ifelse(is.null(names(treelist)),
                     paste("tree", 1:length(treelist), sep = ""),
                     names(treelist))

  # ---- Comparing tip labels.
  eq <- compareTiplabs(treelist)
  # table(eq)
  eq[eq] <- 1
  eq[upper.tri(x = eq, diag = TRUE)] <- 0
  # table(eq)
  c1 <- rep(namelist, times = ncol(eq))
  c2 <- rep(namelist, each = nrow(eq))
  c3 <- unlist(as.data.frame(eq))
  eq_3 <- data.frame(c1, c2, c3)


  if (sum(eq_3$c3) > 1) {
    eq_3 <- eq_3[eq_3$c3 == 1, ]
    t <- igraph::graph_from_data_frame(eq_3[, 1:2], directed = FALSE)
    # class(t)
    # plot(t, vertex.size = 2, edge.width = 2, edge.lty = c("solid"))
    community <- split(names(V(t)), components(t)$membership)
    tree_types <- length(community)
    types_num <- lengths(community, use.names = TRUE)
    num_max <- max(types_num)
    if (sum(types_num == num_max) == 1) {
      POP_tree <- unlist(community[types_num == num_max])
      mtree <- treelist[POP_tree]
      nm_s <- stringr::str_split(names(mtree), " ", simplify = TRUE)
      attr(mtree, "distance") <- nm_s[1]
      attr(mtree, "mtd") <- nm_s[3]
      attr(mtree, "number") <- 1
      return(mtree)
    } else {
      mtrees <- list()
      for(n in 1:sum(types_num == num_max)) {
        POP_tree_tmp <- unlist(community[types_num == num_max][n])
        mtrees_tmp <- treelist[POP_tree_tmp][[1]]

        nm_s <- stringr::str_split(names(mtrees_tmp), " ", simplify = TRUE)
        attr(mtrees_tmp, "distance") <- nm_s[1]
        attr(mtrees_tmp, "mtd") <- nm_s[3]
        attr(mtrees_tmp, "number") <- num_max
        mtrees[[n]] <- mtrees_tmp
      }
      return(mtrees)
    }
  } else {
    stop("No identical tree exists! Please consider another function findBitrees().")
  }
}


# data(GSE45719_268_count)
# processed_data <- getPPdata(GSE45719_268_count)
# d <- getDlist(x = t(processed_data), mtd = c("maximum", "euclidean", "manhattan",
#                                            "minkowski", "chebyshev", "sorensen",
#                                            "gower", "soergel", "kulczynski_d",
#                                            "canberra", "lorentzian", "intersection",
#                                            "non-intersection", "wavehedges", "czekanowski"))
#
# b <- getBasicPartitions(d, method = "all")
#
# m <- findPOPtrees(b$partition)
# attr(m[[1]], "number")




