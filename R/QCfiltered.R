





#' Quality control in scTreeMerge
#'
#' @param data an object of class \code{matrix} with cells on columns and genes on rows.
#' @param meta an object of class \code{matrix}. It contains at least two columns: `cell_id` and `cell_type`.
#'
#' @return an object of class \code{list}
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment counts reducedDim colData
#' @importFrom scater addPerCellQC isOutlier
#' @importFrom DropletUtils emptyDrops
#' @importFrom BiocSingular runPCA
#' @importFrom scDblFinder computeDoubletDensity

QCfiltered <- function(data, meta) {

  # ----------------------------------
  # *********** cell-level ***********
  # ----------------------------------

  cat("******* cell-level *******")

  # create a SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(data),

                                                                  logcounts = log2(as.matrix(data) + 1)),

                                                    colData = meta)

  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)

  sce <- scater::addPerCellQC(sce, detection_limit = 0, flatten = FALSE)


  # Step-1 empty droplets

  cat("Step-1: empty droplets")

  x <- try({set.seed(2022)

  e_out <- DropletUtils::emptyDrops(counts(sce))

  sce <- sce[, which(e_out$FDR <= 0.001)]}, silent = TRUE)

  if ("try-error" %in% class(x)) {

    message("No empty droplets in matrix!")

  } else {

    message(length(e_out$FDR <= 0.001), " empty droplets have been excluded!")

  }

  # A rough filtration
  sce <- sce[, colData(sce)$sum < 200]


  # Step-2 doublets

  cat("Step-2: doublets")

  # with simulation
  set.seed(2022)

  sce <- BiocSingular::runPCA(sce)

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

  # discard <- qc.lib | qc.nexprs | qc.spike | qc.mito
  discard <- reads_drop | feature_drop | mito_drop

  # Summarize the number of cells removed for each reason.
  DataFrame(LibSize = sum(reads_drop), NExprs = sum(feature_drop),

            MitoProp = sum(qc.mito), Total = sum(discard))

  message(sum(discard), " low-quality cells have been excluded!")

  sce <- sce[ , !discard]


  # ----------------------------------
  # *********** gene-level ***********
  # ----------------------------------

  cat("******* gene-level *******")

  names(rowData(sce)) <- "gene_short_name"

  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))

  # featureNames(pd) <- rownames(pd)

  fd <- new("AnnotatedDataFrame", data = as.data.frame(rowData(sce)))

  # featureNames(fd) <- rowData(sce)$gene_short_name

  cds <- newCellDataSet(counts(sce), featureData = fd, phenoData = pd,

                        lowerDetectionLimit = 0, expressionFamily = negbinomial.size())

  # how many cells each feature in a CellDataSet object that are detectably expressed.
  cds <- detectGenes(cds, min_expr = 0)

  # keep genes expressed more than five percent of cells
  genes_drop <- row.names(subset(fData(cds), num_cells_expressed < dim(cds)[2]*0.05))

  message(sum(genes_drop), " genes have been excluded!")

  cds <- cds[!expr_genes, ]



  res <- list(data_qc = exprs(cds),

              meta_qc = pData(cds)[, 1:length(meta)])

  return(res)

}
