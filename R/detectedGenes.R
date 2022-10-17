
#' detectGenes from monocle
#'
#' @importFrom Matrix rowSums

detectGenes <- function (cds, min_expr = NULL)
{
  if (is.null(min_expr)) {
    min_expr <- cds@lowerDetectionLimit
  }
  fData(cds)$num_cells_expressed <- Matrix::rowSums(exprs(cds) >
                                                      min_expr)
  pData(cds)$num_genes_expressed <- Matrix::colSums(exprs(cds) >
                                                      min_expr)
  cds
}
