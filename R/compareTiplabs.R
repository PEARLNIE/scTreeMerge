
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@


#' @title Tip labels comparisons
#' @description This function is to compare tip labels of each trees.
#'
#' @param treelist an object of class \code{multiPhylo}.
#'
#' @return an object of class \code{matrix}. \code{TRUE} indicates that two trees are of the same type.because of the identical tip labels.
#' @export
#' @importFrom ape all.equal.phylo
#'




compareTiplabs <- function(treelist) {

  message(paste("###------------", "*START*", "------------###", sep = " "))
  eq <- matrix(nrow = length(treelist), ncol = length(treelist))
  colnames(eq) <- names(treelist)
  rownames(eq) <- names(treelist)

  for(i in 1:length(treelist)) {
    print(paste("----", names(treelist[i]), "----", sep = " "))
    tree1 <- treelist[[i]]
    for(n in 1:length(treelist)) {
      tree2 <- treelist[[n]]
      tmp <- all.equal.phylo(target = tree1, current = tree2, use.edge.length = FALSE,
                             use.tip.label = TRUE, index.return = FALSE,
                             tolerance = .Machine$double.eps ^ 0.5, scale = NULL)
      eq[i, n] <- tmp
    }
  }
  message(paste("###------------", "*END*", "------------###", sep = " "))

  return(eq)
}


# data(GSE45719_268_count)
# processed_data <- getPPdata(GSE45719_268_count)
# d <- getDlist(x = t(processed_data), mtd = "euclidean")
# b <- getBasicPartitions(d, method = "all")
#
# eq <- compareTiplabs(b$partition)
# head(eq)

