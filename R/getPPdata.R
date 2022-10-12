

## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@


#' @title Data preprocessing
#' @description This function is used to preprocess input data.
#' @param x an object of class \code{matrix}. Each column corresponds to a sample and each row to a variable.
#' @param nfeatures a numeric value represents the number of the filtered features.
#' @return an object of class \code{matrix}.
#' @export
#' @importFrom Seurat CreateSeuratObject FindVariableFeatures VariableFeatures
#'



getPPdata <- function(x, nfeatures = 2000) {

  # Error checking.
  # if (!inherits(x, c("data.frame", "matrix")))
  #   stop("x must be object of class 'data.frame' or 'matrix'.")
  if (!inherits(x, "matrix"))

    x <- as.matrix(x)

  rownames(x) <- sapply(X = rownames(x), function(i) gsub(pattern = "\\|", replacement = "-", x = i))



  if (!is.null(nfeatures)) {

    # if (!require("Seurat")) BiocManager::install("Seurat")
    # suppressPackageStartupMessages(library(Seurat))
    seu <- Seurat::CreateSeuratObject(counts = x,

                                      min.cells = 3, # Retain genes that are expressed in at least three cells

                                      min.features = 200,)

    seu2 <- Seurat::FindVariableFeatures(seu,

                                         selection.method = "vst",

                                         nfeatures = nfeatures)

    seu_hvgs <- Seurat::VariableFeatures(seu2)
  }

  res <- x[seu_hvgs, ]
}


# data(GSE45719_268_count)
# processed_data <- getPPdata(GSE45719_268_count)
# dim(processed_data)



