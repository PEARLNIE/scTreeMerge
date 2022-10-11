




#' Normalization method for scRNA-seq
#'
#' @param data an object of class \code{matrix} with cells on columns and genes on rows.
#' @param norm_method a character string. This could be one of \code{"CountToCPM"}, \code{"CountToFPKM"}, \code{"CountToCPM"}, \code{"CPMToFPKM"} or \code{"FPKMToCPM"}.
#' @param gene_length an numberic vetcor of gene length.
#' @param species a character string.It could be one of \code{"huaman"} or \code{"mouse"}
#'
#' @return a normalized matrix.
#'
#' @export

Norm_sc <- function(data, norm_method = "CountToCPM",

                          gene_length = NULL, species = "human") {

  ## gene length data
  if(is.null(gene_length)) {

    if(species == "human") {

      gene <- scTreeMerge::hs_gene_length

    }

    if(species == "mouse"){

      gene <- scTreeMerge::mm_gene_length

    }
  }

  ## gene form
  gene_name <- rownames(data)[1:10]

  if(all(startsWith(gene_name, "ENS"))) {

    gene_formal <- "ensembl_id"

  } else {

    if(all(grepl(pattern = "^[1-9]", gene_name))) {

      gene_formal <- "gene_id"

    } else {

      gene_formal <- "symbol"

    }
  }

  ## filter miss-match gene
  intersect_gene <- intersect(rownames(data), gene[, gene_formal])

  message(paste0(nrow(data)-length(intersect_gene), " genes are removed due to the miss-match."))

  data <- data[intersect_gene, ]

  index <- match(intersect_gene, gene[, gene_formal])

  gene_length <- gene$length[index]




  # -------------- CountToCPM --------------

  CountToCPM <- function(counts) {

    N <- sum(counts)

    CPM <- counts / N * 10^6

    return(CPM)

  }



  # -------------- CountToFPKM --------------

  CountToFPKM <- function(counts, gene_length) {

    N <- sum(counts)

    CPM <- counts / N * 10^6

    FPKM <- CPM / gene_length * 10^3

    return(FPKM)
  }




  # -------------- CountToTPM --------------

  CountToTPM <- function(counts, gene_length) {

    counts_by_length <- counts / gene_length * 10^3

    N <- sum(counts_by_length)

    TPM <- counts_by_length / N * 10^6

    return(TPM)

  }



  # -------------- CPMToFPKM --------------

  CPMToFPKM <- function(CPM, gene_length) {

    FPKM <- CPM / gene_length * 10^3

    return(FPKM)

  }



  # -------------- FPKMToCPM --------------

  FPKMToCPM <- function(FPKM, gene_length) {

    CPM <- FPKM * gene_length / 10^3

    return(CPM)

  }






  if(norm_method == "CountToCPM") {

    norm_data <- apply(data, 2, CountToCPM)

  }

  if(norm_method == "CountToFPKM") {

    norm_data <- apply(data, 2, CountToFPKM, gene_length = gene_length)

  }

  if(norm_method == "CountToCPM") {

    norm_data <- apply(data, 2, CountToTPM, gene_length = gene_length)

  }

  if(norm_method == "CPMToFPKM") {

    norm_data <- apply(CPM, 2, CPMToFPKM, gene_length = gene_length)

  }

  if(norm_method == "FPKMToCPM") {

    norm_data <- apply(FPKM, 2, FPKMToCPM, gene_length = gene_length)

  }

  return(norm_data)

}
