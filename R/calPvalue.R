


## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@



#' Calculating P values of individual clsters
#'
#' This function is to calculate the average P values of individual clusters.
#'
#' @param x an object of class \code{'matrix'}. Each column corresponds to a sample and each row to a variable.
#' @param clusters an integer vector of clustering labels.
#' @param nsim integer value. the number of simulation. Default is 10.
#'
#' @return the average of P values.
#' @export
#' @importFrom stats shapiro.test wilcox.test t.test
#'
#' @examples
#' require(ape)
#' data(GSE45719_268_count)
#' processed_data <- getPPdata(GSE45719_268_count)
#' d1 <- getDlist(x = t(processed_data))
#' b1 <- getBasicPartitions(d1)
#' cl <- cutree(as.hclust.phylo(b1$partition$`euclidean + complete`), 3)
#' meanp <- calPvalue(x = processed_data, clusters = cl, nsim = 10)



calPvalue <- function(x, clusters, nsim = 1000) {

  # Error checking.
  # if (!inherits(x, c("data.frame", "matrix")))
  #   stop("x must be object of class 'data.frame' or 'matrix'.")
  if (!inherits(x, c("matrix")))
    stop("x must be object of class 'matrix'.")

  ClustAssign <- as.character(clusters)
  n <- length(unique(ClustAssign))
  pval <- matrix(NA, ncol = n, nrow = n)
  rownames(pval) <- c(paste("cluster", c(1:n), sep = ""))
  colnames(pval) <- c(paste("cluster", c(1:n), sep = ""))

  # set up data to run.
  x_list <- list()
  l_list <- list()
  for (i in 1: n) {
    x_list[[paste0("cluster", i)]] <- as.matrix(x[, ClustAssign == i])
    colnames(x_list[[paste0("cluster", i)]]) <- colnames(x)[ClustAssign == i]
    l_list[[i]] <- ClustAssign[ClustAssign == i]
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      x1 <- x_list[[i]]
      x2 <- x_list[[j]]
      l1 <- l_list[[i]]
      l2 <- l_list[[j]]
      # Run sigclust
      sig_dat <- as.matrix(t(cbind(x1, x2)))
      l1 <- rep(1, length(l1))
      l2 <- rep(2, length(l2))
      sig_label <- c(l1, l2)
      if (length(sig_label) > 2) {
        sig <- sigclust_n(x = sig_dat,
                          nsim = nsim,
                          label = sig_label)$pval
      } else {
        ntest1 <- shapiro.test(sig_dat[1,])$p.value
        ntest2 <- shapiro.test(sig_dat[2,])$p.value
        if (ntest1 < 0.05|ntest1 < 0.05) {
          sig <- wilcox.test(x = sig_dat[1,], y = sig_dat[2,])$p.value
        } else {
          # res.ftest <- var.test(x = sig_dat[1,], y = sig_dat[2,])$statistic
          sig <- t.test(x = sig_dat[1,], y = sig_dat[2,])$p.value
        }
      }
      pval[j, i] <- sig
    }
  }
  pval <- as.vector(pval)
  pval <- pval[!is.na(pval)]
  meanp <- sum(pval)/length(pval)
  meanp
}
