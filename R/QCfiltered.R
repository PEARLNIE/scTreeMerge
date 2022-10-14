





#' Quality control in scTreeMerge
#'
#' @param data an object of class \code{matrix} with cells on columns and genes on rows.
#' @param meta an object of class \code{matrix}. It contains at least two columns: `cell_id` and `cell_type`.
#'
#' @return an object of class \code{list}
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scater addPerCellQC isOutlier
#' @importFrom monocle newCellDataSet detectGenes

QCfiltered <- function(data, meta) {

  message("cell-level")

  # create a SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(data),

                                                                  logcounts = log2(as.matrix(data) + 1)),

                                                    colData = meta)

  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)

  sce <- scater::addPerCellQC(sce, detection_limit = 0, flatten = FALSE)


  # ----------------------------------
  # *********** cell-level ***********
  # ----------------------------------


  # Over 4, the minimum number of MADs away from median required for a value to be called an outlier.
  reads_drop <- scater::isOutlier(as.numeric(sce$sum), nmads = 4, type = "lower", log = TRUE)

  feature_drop <- scater::isOutlier(sce$detected, nmads = 4, type = "lower", log = TRUE)


  keep_samples <- !(reads_drop | feature_drop)

  keep_samples[which(is.na(keep_samples))] <- FALSE

  sce <- sce[ , keep_samples]


  # ----------------------------------
  # *********** gene-level ***********
  # ----------------------------------

  message("gene-level")

  library(monocle)

  names(rowData(sce)) <- "gene_short_name"

  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(sce)))

  # featureNames(pd) <- rownames(pd)

  fd <- new("AnnotatedDataFrame", data = as.data.frame(rowData(sce)))

  # featureNames(fd) <- rowData(sce)$gene_short_name





  cds <- monocle::newCellDataSet(counts(sce), featureData = fd, phenoData = pd,

                                 lowerDetectionLimit = 0, expressionFamily = negbinomial.size())

  # how many cells each feature in a CellDataSet object that are detectably expressed.
  cds <- monocle::detectGenes(cds, min_expr = 0)

  # keep genes expressed more than five percent of cells
  expr_genes <- row.names(subset(fData(cds), num_cells_expressed >= dim(cds)[2]*0.05))

  cds <- cds[expr_genes, ]


  # In each types, keep genes expressed more than five percent of cells
  remove_genes <- list()

  for(i in unique(pData(cds)$cell_type)) {

    cds_tmp <- cds[, rownames(subset(pData(cds), cell_type == i))] %>%

      monocle::detectGenes(min_expr = 0)

    remove_genes[[i]] <- row.names(subset(fData(cds_tmp), num_cells_expressed < dim(cds_tmp)[2]*0.05))

  }

  remove_genes_all <- unique(unlist(remove_genes))

  keep_genes <- setdiff(rownames(fData(cds)), remove_genes_all)

  cds <- cds[keep_genes, ]

  dim(cds)




  res <- list(data_qc = exprs(cds),

              meta_qc = pData(cds)[, 1:length(meta)])

  return(res)

}
