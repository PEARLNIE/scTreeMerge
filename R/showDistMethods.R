

## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@



#' @title disatnce measures used in scTreeMerge
#' @description It returns the names of the measures that can be applied to compute distances using the \code{getDlist} function.
#' @export
#' @examples
#' showDistMethods()


showDistMethods <- function() {

  n <- c("maximum",
         "euclidean",
         "manhattan",
         "minkowski",
         "chebyshev",
         "sorensen",
         "gower",
         "soergel",
         "kulczynski_d",
         "canberra",
         "lorentzian",
         "intersection",
         "non-intersection",
         "wavehedges",
         "czekanowski",
         "motyka",
         "kulczynski_s",
         "tanimoto",
         "ruzicka",
         "inner_product",
         "harmonic_mean",
         "cosine",
         "hassebrook",
         "jaccard",
         "dice",
         "fidelity",
         "bhattacharyya",
         "hellinger",
         "matusita",
         "squared_chord",
         "squared_euclidean",
         "pearson",
         "neyman",
         "squared_chi",
         "prob_symm",
         "divergence",
         "clark",
         "additive_symm",
         "kullback-leibler",
         "jeffreys",
         "k_divergence",
         "topsoe",
         "jensen-shannon",
         "jensen_difference",
         "taneja",
         "kumar-johnson",
         "avg")
  return(n)
}

