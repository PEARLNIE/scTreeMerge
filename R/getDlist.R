


## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@



#' @title Distance computation with lots of disatnce measures.
#' @description  This is a distance matrix computation function which containins 47 measures. It computes between the rows of a matrix and return the distance matrix computed by choosing the specific measure.
#' @param x an object of class \code{matrix}. Each column corresponds to a sample and each row to a variable.
#' @param mtd a character string that represents the distance measure to be used. It could be one or some of the 47 measures(\code{showDistMethods()}). If users want to use all the measures, \code{all} could be selected in this function. For more information, please see Details.
#' @param p The power of Minkowski distance. Note that setting p = 1 is equivalent to calculating the Manhattan distance and setting p = 2 is equivalent to calculating the Euclidean distance. Default: \code{p = 2}.
#' @details This function implements 47 measures to quantify the disatnce between two objects:
#' \itemize{
#' \item Forty-six measures from Philentropy package
#' \itemize{
#' \item euclidean
#' \item manhattan
#' \item minkowski
#' \item chebyshev
#' \item sorensen
#' \item gower
#' \item soergel
#' \item kulczynski_d
#' \item canberra
#' \item lorentzian
#' \item intersection
#' \item non-intersection
#' \item wavehedges
#' \item czekanowski
#' \item motyka
#' \item kulczynski_s
#' \item tanimoto
#' \item ruzicka
#' \item inner_product
#' \item harmonic_mean
#' \item cosine
#' \item hassebrook
#' \item jaccard
#' \item dice
#' \item fidelity
#' \item bhattacharyya
#' \item hellinger
#' \item matusita
#' \item squared_chord
#' \item squared_euclidean
#' \item pearson
#' \item neyman
#' \item squared_chi
#' \item prob_symm
#' \item divergence
#' \item clark
#' \item additive_symm
#' \item kullback-leibler
#' \item jeffreys
#' \item k_divergence
#' \item topsoe
#' \item jensen-shannon
#' \item jensen_difference
#' \item taneja
#' \item kumar-johnson
#' \item avg
#' }
#'
#' \item One measure from stats package
#' \itemize{
#' \item maximum
#' }
#' }
#' @return an object of class \code{"dist"} or an object if class \code{"list"}.
#' @export
#' @importFrom philentropy distance
#' @importFrom stats dist
#' @examples
#' data(GSE45719_268_count)
#' processed_data <- getPPdata(GSE45719_268_count)
#'
#' d <- getDlist(x = t(processed_data), mtd = "euclidean")
#' class(d)
#' d <- getDlist(x = t(processed_data), mtd = "minkowski", p = 2)
#' class(d)
#' # make sure that all vectors sum up to 1.0 ...
#' d <- getDlist(x = t(processed_data), mtd = "matusita")
#'
#' d <- getDlist(x = t(processed_data), mtd = c("euclidean", "manhattan"))
#' class(d)
#' d[[1]]
#' class(d[[1]])
#'
#' d <- getDlist(x = t(processed_data), mtd = "all")
#' class(d)


# Calculating with many distance metrics.
getDlist <- function(x, mtd = "euclidean", p = 2) {

  # Error checking.
  # if (!inherits(x, c("data.frame", "matrix")))
  #   stop("x must be object of class 'data.frame' or 'matrix'.")
  if (!inherits(x, "matrix"))

    x <- as.matrix(x)

  if (all(length(mtd) > 1 & "all" %in% mtd))

    stop("'all' covers all 47 measures in this function. Once 'all' is chosen,
         there's no need to choose any of the 47 measures.")


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

  if ("all" %in% mtd)

    mtd <- n

  # Calculating the specific distance measures.
  res <- list()

  q <- 0

  for (i in mtd) {

    message(paste("[------", i, "-------]", sep = " "))

    if (i == "maximum") {

      q <- q + 1

      res[[q]] <- dist(x = x, method = i, diag = FALSE, upper = FALSE, p = p)

    } else {

      q <- q + 1

      if (i %in% c("kullback-leibler", "k_divergence", "matusita")) {

        x <- t(apply(x, 1, function(i) {i/sum(i)}))

        res[[q]] <- philentropy::distance(x = x, method = i, p = p, test.na = TRUE,

                                          unit = "log", est.prob = NULL,

                                          use.row.names = TRUE, as.dist.obj = TRUE,

                                          diag = FALSE, upper = FALSE,mute.message = TRUE)

      } else {

        res[[q]] <- philentropy::distance(x = x, method = i, p = p, test.na = TRUE,

                                          unit = "log", est.prob = NULL,

                                          use.row.names = TRUE, as.dist.obj = TRUE,

                                          diag = FALSE, upper = FALSE,mute.message = TRUE)
      }

    }
  }

  names(res) <- mtd

  if (q == 1) {

    return(res[[1]])

  } else {

    return(res)
  }

}













