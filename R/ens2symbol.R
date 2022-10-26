
#' Ensemble to gene symbol
#'
#' @import EnsDb.Mmusculus.v79
#' @importFrom AnnotationDbi select
#' @importFrom dplyr add_count
#' @importFrom stats aggregate
#'
#' @export

ens2symbol <- function(data, species = "mouse") {

  if(species == "mouse") {

    # columns(EnsDb.Mmusculus.v79)
    #
    # keytypes(EnsDb.Mmusculus.v79)

    # run the query
    annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                                   keys = rownames(data),
                                   columns = "SYMBOL",
                                   keytype = "GENEID")


    annot <- dplyr::add_count(annot, SYMBOL)
    # annot_filter <- dplyr::filter(annot, n>1)

    zero_drop <- sapply(annot$SYMBOL, nchar) == 0

    annot <- annot[!zero_drop, ]

    annot_keep <- match(annot$GENEID, rownames(data))

    data <- data[annot_keep, ]

    rownames(data) <- annot$SYMBOL

    data <- aggregate(data, by = list(rownames(data)), mean)

    rownames(data) <- data[, 1]

    data <- as.matrix(data[, -1])

    return(data)

  }

}







