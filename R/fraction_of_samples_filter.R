#' Fraction of samples filter
#'
#' Each circRNA must be detected in at least the specified fraction of samples
#'
#' @param x the matrix to be filtered. It must have a \code{circ_id} column
#' @param rthr the minimum BJR count for a circRNA to be kept
#' @param frac the fraction of the samples in which each circRNA must be
#'   detected. Default is that each circRNA must be detected in al least half
#'   the samples
#'
#' @return the matrix bearing circRNAs detected in at least a fraction of the
#'   samples
#' @export
#'
#' @examples
#' @import data.table
half_the_samples_filter <-
  function(x, rthr = 0, frac = .5) {
    xmat <- as.matrix(data.frame(x, row.names = "circ_id"))
    x[circ_id %in%
        rownames(xmat)[rowSums(xmat >= rthr) >=
                         ceiling(length(colnames(xmat)) * frac)]]
  }
