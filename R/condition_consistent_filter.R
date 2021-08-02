#' Condition consistency filter
#'
#' CircRNAs must be detected in all the samples of at least one condition
#'
#' @param x the matrix to be filtered. It must have a \code{circ_id} column
#' @param rthr the minimum BJR count for a circRNA to be kept
#' @param cond the sample metadata table. Must have \code{sample_id} and
#'   \code{condition} columns
#'
#' @return the matrix bearing circRNAs expressed in all samples of at least one
#'   condition
#' @export
#'
#' @examples
#' @import data.table
condition_consistent_filter <-
  function(x, rthr, cond) {

    circ_ids_to_keep <-
      merge(cond[, .(samples_per_cond = .N),
                 by = condition],
            merge(melt(x,
                       id.vars = "circ_id",
                       variable.name = "sample_id",
                       value.name = "BJR"),
                  cond,
                  by = "sample_id")[BJR >= rthr
                                    ][, .N,
                                      by = .(condition,
                                             circ_id)],
            by = "condition")[samples_per_cond == N,
                              unique(circ_id)]

    x[circ_id %in% circ_ids_to_keep]
  }
