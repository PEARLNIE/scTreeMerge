

#' newCellDataSet from monocle
#'
#' @importFrom VGAM negbinomial.size


newCellDataSet <- function(cellData,
                           phenoData = NULL,
                           featureData = NULL,
                           lowerDetectionLimit = 0.1,
                           expressionFamily = VGAM::negbinomial.size())
{
  if (!("gene_short_name" %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }
  if (class(cellData) != "matrix" && isSparseMatrix(cellData) ==
      FALSE) {
    stop("Error: argument cellData must be a matrix (either sparse from the Matrix package or dense)")
  }
  if (!("gene_short_name" %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }
  sizeFactors <- rep(NA_real_, ncol(cellData))
  if (is.null(phenoData))
    phenoData <- annotatedDataFrameFrom(cellData, byrow = FALSE)
  if (is.null(featureData))
    featureData <- annotatedDataFrameFrom(cellData, byrow = TRUE)
  if (!("gene_short_name" %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }
  phenoData$Size_Factor <- sizeFactors
  cds <- new("CellDataSet", assayData = assayDataNew("environment",
                                                     exprs = cellData), phenoData = phenoData, featureData = featureData,
             lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily,
             dispFitInfo = new.env(hash = TRUE))
  validObject(cds)
  cds
}
