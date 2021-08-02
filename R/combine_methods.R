#' Combine different methods' expression estimates into one expression matrix
#'
#' @param x the list of the methods' output to be combined or the path of the
#'   CirComPara2's results
#' @param filtfun a function to filter the rows in the combined matrix
#' @param name_sep the character(s) used to concatenate SAMPLE and METHOD names
#' @param as_data_table if TRUE, return a data.table instead of a matrix
#'
#' @return a matrix of the expression values with samples in columns and
#'   features in the rows. Column names have the form SAMPLE.METHOD.
#' @export
#'
#' @examples
#' mat <- get_combined_matrix(x = "/home/user/CirComPara2_analysis")
#' @import data.table
get_combined_matrix <-
  function (x,
            filtfun = NULL,
            name_sep = ".",
            as_data_table = F) {

    gtf <-
      fread(file = file.path(x,
                      "circular_expression", "circrna_collection",
                      "circrnas.gtf"),
            header = F,
            select = c(1, 2, 4, 5, 6, 7,
                       9))[, .(Method = V2,
                               BJR = sum(V6)),
                           by = .(sample_id = sub(".*sample_id \"([^\"]+)\";.*",
                                                  "\\1", V9),
                                  circ_id = paste0(V1, ":", V4, "-", V5))]

    mat_dt <-
      dcast(data = gtf,
            formula = circ_id ~ sample_id + Method,
            sep = name_sep,
            fill = 0,
            value.var = "BJR")

    if(as_data_table) {
      mat_dt
    } else {
      as.matrix(data.frame(mat_dt, row.names = "circ_id"))
    }
  }
