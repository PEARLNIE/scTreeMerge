#' The Stree class
#'
#' S4 class that holds parameters for the supertree.
#'
#' @section Parameters:
#'
#' The Stree class defines the following parameters:
#'
#' \describe{
#'     \item{\code{raw}}{The matrix of P values.}
#'     \item{\code{tree}}{The structure of the supertree.}
#'     \item{\code{k}}{The integer value of the optimal cluster number.}
#'     \item{\code{pvalue}}{The numeric value of individual clsters with the optimal cluster number.}
#' }
#'
#' @name Stree
#' @rdname Stree
#' @aliases Stree-class

setOldClass("phylo")
setClass(Class = "Stree",
         slots = list(raw = "matrix",
                      tree = "phylo",
                      k = "integer",
                      pvalue = "numeric"))

