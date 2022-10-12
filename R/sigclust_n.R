

## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@
## -----------------------------------------------------------------
##     ***       ** ***      *** **********   / DATE ：2021/02/27/
##    ** **     **    **  **    ***          /
##   **   **   **      **      **********   / AUTHOR ：Xiner Nie
##  **     ** **    **  **    ***          /
## **       *** ***      *** ***********  / AUTHOR ：Bo Li
## -----------------------------------------------------------------
## @.@-@.@ @.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@.@



#' Modified sigclust
#'
#' This is a modified function.
#'
#' @param x an object of class \code{matrix}. Each column corresponds to a variable and each row to a sample.
#' @param nsim integer value. The number of simulation.
#' @param label a numeric, integer vector of 1s and 2s with length \code{nrow(x)} which indicates given cluster labels.
#' @return an object of class \code{list}.
#'
#' @export
#' @importFrom sigclust .vareigen .simnull
#' @importFrom stats sd pnorm
#'



sigclust_n <- function(x, nsim, label)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (n > 1) {
    x <- as.matrix(x)
    # if (labflag == 0) {
    #   xclust <- .cluster(x, n, p)
    #   for (i in 1:nrep) {
    #     clust.temp <- .cluster(x, n, p)
    #     if (clust.temp$cindex < xclust$cindex)
    #       xclust <- clust.temp
    #     xcindex <- xclust$cindex
    #   }
    # }
    xvareigen <- sigclust::.vareigen(x, n, p, icovest = 3)
    simcindex <- rep(0, nsim)
    for (i in 1:nsim) {
      xsim <- sigclust::.simnull(xvareigen$vsimeigval, n, p)
      simcindex[i] <- xsim$cindex
    }
    # if (labflag == 0) {
    #   index <- (simcindex <= xclust$cindex)
    #   mindex <- mean(simcindex)
    #   sindex <- sd(simcindex)
    #   pval <- sum(index)/nsim
    #   pvalnorm <- pnorm(xclust$cindex, mindex, sindex)
    # }
    # if (labflag == 1) {
    if (length(label[label == 1]) > 1) {
      meanp1 <- colMeans(x[label == 1, ])
      txdiff1 <- t(x[label == 1, ]) - meanp1
    } else {
      meanp1 <- x[label == 1, ]
      txdiff1 <- t(x[label == 1, ]) - meanp1
    }
    if (length(label[label == 2]) > 1) {
      meanp2 <- colMeans(x[label == 2, ])
      txdiff2 <- t(x[label == 2, ]) - meanp2
    } else {
      meanp2 <- x[label == 2, ]
      txdiff2 <- t(x[label == 2, ]) - meanp2
    }

    withinsum <- sum(txdiff1^2) + sum(txdiff2^2)
    meanp <- colMeans(x)
    tx <- t(x)
    txdiff <- tx - meanp
    totalsum <- sum(txdiff^2)
    cindexlab <- withinsum/totalsum
    index <- (simcindex <= cindexlab)
    mindex <- mean(simcindex)
    sindex <- sd(simcindex)
    pval <- sum(index)/nsim
    pvalnorm <- pnorm(cindexlab, mindex, sindex)
    xcindex <- cindexlab
    # }

    return(list(raw.data = x,
                veigval = xvareigen$veigval,
                vsimeigval = xvareigen$vsimeigval,
                simbackvar = xvareigen$simbackvar,
                nsim = nsim,
                simcindex = simcindex,
                pval = pval,
                pvalnorm = pvalnorm,
                xcindex = xcindex))
  } else {
    print("Only one sample left, no need for clustering!")
    return(0)
  }
}


# data <- matrix(rnorm(1500*150),1500,150, dimnames = list(1:1500, 1:150))
# label <- setNames(object = rep(1:2, each = 75), colnames(data))
# p <- sigclust_n(data, 10, label)
# p$pval





