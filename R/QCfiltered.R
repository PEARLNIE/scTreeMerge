





#' Quality control in scTreeMerge
#'
#' @param data an object of class \code{matrix} with cells on columns and genes on rows.
#' @param meta an object of class \code{matrix}. It contains at least two columns: `cell_id` and `cell_type`.
#' @param species a character string.It could be one of \code{"huaman"} or \code{"mouse"}
#'
#' @return an object of class \code{list}
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment counts reducedDim
#' @importFrom scater addPerCellQC isOutlier runPCA
#' @importFrom DropletUtils emptyDrops
#' @importFrom scDblFinder computeDoubletDensity
#' @importFrom SummarizedExperiment rowData colData rowRanges
#' @importFrom GenomeInfoDb seqnames


QCfiltered <- function(data, meta, species) {

  # ----------------------------------
  # *********** cell-level ***********
  # ----------------------------------

  cat("******* cell-level *******", fill = TRUE)

  # create a SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(data),

                                                                  logcounts = log2(as.matrix(data) + 1)),

                                                    colData = meta)

  # define feature names in feature_symbol column
  SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)

  # Identifying the mitochondrial transcripts in our SingleCellExperiment.
  location <- SummarizedExperiment::rowRanges(sce)

  if(species == "huamn") {

    is_mito <- any(GenomeInfoDb::seqnames(location) == "MT")

  } else {

    is_mito <- any(GenomeInfoDb::seqnames(location) == "Mt")

  }

  sce <- scater::addPerCellQC(sce, subsets = list(Mito = is_mito), detection_limit = 0, flatten = FALSE)



  # Step-1 empty droplets

  cat("Step-1: empty droplets", fill = TRUE)

  x <- try({set.seed(2022)

  e_out <- DropletUtils::emptyDrops(counts(sce))

  sce <- sce[, which(e_out$FDR <= 0.001)]}, silent = TRUE)

  if ("try-error" %in% class(x)) {

    message("No empty droplets in matrix!")

  } else {

    message(length(e_out$FDR <= 0.001), " empty droplets have been excluded!")

  }

  # A rough filtration
  sce <- sce[, colData(sce)$sum >= 200]


  # Step-2 doublets

  cat("Step-2: doublets", fill = TRUE)

  # with simulation
  set.seed(2022)

  sce <- scater::runPCA(sce)

  set.seed(2023)

  dbl_dens <- scDblFinder::computeDoubletDensity(sce,

                                                 d = ncol(reducedDim(sce)))

  sce$DoubletScore <- dbl_dens

  dbl_out <- scater::isOutlier(dbl_dens)

  # which(dbl.dens > 3)

  sce <- sce[, !dbl_out]

  message(sum(dbl_out), " doublets have been excluded!")


  # Step-3 Strict filtering

  # Over 3, the minimum number of MADs away from median required for a value to be called an outlier.
  reads_drop <- scater::isOutlier(as.numeric(sce$sum), nmads = 3, type = "lower", log = TRUE)

  feature_drop <- scater::isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)

  mito_drop <- colData(sce)$subsets$Mito$percent > 10

  discard <- reads_drop | feature_drop | mito_drop

  # Summarize the number of cells removed for each reason.
  # DataFrame(LibSize = sum(reads_drop), NExprs = sum(feature_drop),
  #
  #           MitoProp = sum(mito_drop), Total = sum(discard))

  message(sum(discard), " low-quality cells have been excluded!")

  sce <- sce[ , !discard]


  # ----------------------------------
  # *********** gene-level ***********
  # ----------------------------------

  cat("******* gene-level *******", fill = TRUE)

  num_cells_expressed <- scater::nexprs(sce, byrow = TRUE)

  # keep genes expressed more than five percent of cells
  genes_drop <- num_cells_expressed < dim(sce)[2]*0.05

  message(sum(genes_drop), " genes have been excluded!")

  sce <- sce[!genes_drop, ]



  res <- list(data_qc = counts(sce),

              meta_qc = as.data.frame(colData(sce)[, 1:ncol(meta)]))

  return(res)

}
