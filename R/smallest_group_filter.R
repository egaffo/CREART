#' Smallest condition size filter
#'
#' Each circRNA must be detected in at least as many samples as in the smallest
#' condition group
#'
#' @param x the matrix to be filtered. It must have a \code{circ_id} column
#' @param rthr the minimum BJR count for a circRNA to be kept
#' @param cond the sample metadata table. Must have \code{sample_id} and
#'   \code{condition} columns
#'
#' @return a matrix bearing circRNAs expressed in as many samples as the size of
#'   the smallest condition group
#' @export
#'
#' @examples
#' @import data.table
smallest_group_filter <-
  function(x, rthr, cond) {

    smallest_group_n <- cond[, .N, by = condition][, min(N)]

    circ_ids_to_keep <-
      melt(x,
           id.vars = "circ_id",
           variable.name = "sample_id",
           value.name = "BJR")[BJR >= rthr,
                               .N,
                               by = .(circ_id)][N >= smallest_group_n,
                                                unique(circ_id)]
    x[circ_id %in% circ_ids_to_keep]
  }
