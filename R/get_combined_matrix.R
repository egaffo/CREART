#' Combine different methods' expression estimates into one expression matrix
#'
#' @param x the list of the methods' output to be combined, or the path of the
#'   CirComPara2's results, or the full path to the circrnas.gtf file from the
#'   CirComPara's output
#' @param filtfun a function to filter the rows in the combined matrix. No
#'   filter is applied by default. Set "ccp2" to apply a CirComPara2-like
#'   filter, that is keep only circRNAs detected with at least N backspliced
#'   reads and by at least M methods, where N and M are defined by the
#'   \code{min_bjr} and \code{min_methods} parameters.
#' @param name_sep the character(s) used to concatenate SAMPLE and METHOD names
#' @param as_data_table if TRUE, return a data.table instead of a matrix
#' @param min_bjr minimum amount of backspliced reads for a circRNA to be kept
#'   when the ccp2 filter is enabled
#' @param min_methods minimum number of methods that commonly detected a circRNA
#'   to be kept. It applies to the ccp2 filter.
#' @param hard_threshold FALSE by default, meaning that if a circRNA passes the
#'   filtering in one sample, then its values in other samples are preserved
#'   even when they are below the specified threshold. For instance, if
#'   circRNA_x has 2 BJRs according to Method_M1 and 2 BJRs according to
#'   Method_M2 in Sample_A, then it passes the filter. If, in Sample_B,
#'   circRNA_x has 1 BJRs according to Method_M1 and 1 BJRs according to
#'   Method_M2, it would not pass the filter, but the 1 values are kept since
#'   there is some stronger evidence of expression in another sample. Set TRUE
#'   to apply an hard threshold, meaning that the 1 values will be set to zeros
#'   following the reasoning that if there is not enough BJRs, the circRNA is
#'   not considered as being detected. Below an example of the output matrix
#'   with:
#'
#'   \code{hard_threshold = FALSE} (the default)
#'
#'   \tabular{ccccc}{ circ_id \tab S_A.M_M1 \tab S_A.M_M2 \tab S_B.M_M1 \tab
#'   S_B.M_M2 \cr circRNA_x \tab 2 \tab 2 \tab 1 \tab 1 }
#'
#'
#'   or \code{hard_threshold = TRUE}
#'
#'   \tabular{ccccc}{ circ_id \tab S_A.M_M1 \tab S_A.M_M2 \tab S_B.M_M1 \tab
#'   S_B.M_M2 \cr circRNA_x \tab 2 \tab 2 \tab 0 \tab 0 }
#'
#' @param select_methods a character vector with the names of methods to keep.
#'   NULL to keep all methods
#'
#' @return a matrix of the expression values with samples in columns and
#'   features in the rows. Column names have the form SAMPLE.METHOD.
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- get_combined_matrix(x = "/home/user/CirComPara2_analysis")
#' }
#' @import data.table
get_combined_matrix <-
  function(x,
           filtfun = NULL,
           name_sep = ".",
           as_data_table = F,
           min_bjr = 2,
           min_methods = 2,
           hard_threshold = F,
           select_methods = NULL) {

    ## if a directory is given, then compose the path to the file according to
    ## the CirComPara's directory structure
    if (dir.exists(x)) {

      input_file_path <-
        file.path(x,
                  "circular_expression",
                  "circrna_collection",
                  "circrnas.gtf")

      ## try old CirComPara path
      if (!file.exists(input_file_path)) {

        input_file_path <-
          file.path(x,
                    "circular_expression",
                    "circRNA_collection",
                    "circrnas.gtf")
      }
    } else {

      ## the circrnas.gtf full path was given
      if (file.exists(x)) {

        input_file_path <- x
      } else {
        stop(paste0(x, " does not exists or cannot be read."))
      }
    }

    ## read the file
    gtf <-
      fread(file = input_file_path,
            header = F,
            select = c(1, 2, 4, 5, 6, 7,
                       9))[, .(Method = V2,
                               BJR = sum(V6)),
                           by = .(sample_id = sub(".*sample_id \"([^\"]+)\";.*",
                                                  "\\1", V9),
                                  circ_id = paste0(V1, ":", V4, "-", V5))]

    ## select wanted methods
    if (!is.null(select_methods)) {

      gtf <- gtf[Method %in% select_methods]

    }

    ## apply filter (if any)
    if (!is.null(filtfun)) {

      ## apply the CirComPara2-like filter
      if (filtfun == "ccp2") {

        ## select circRNAs with enough BJRs in at least one samples
        ## and detected by enough methods
        keep <-
          gtf[BJR >= min_bjr, .N, by = .(sample_id, circ_id)
              ][N >= min_methods, unique(circ_id)]

          if (hard_threshold) {

            ## filter out circRNAs with not enough BJRs, as if they were
            ## not detected
            gtf <- gtf[BJR >= min_bjr]

          }

        gtf <- gtf[circ_id %in% keep]

      } else {
        message("Only 'ccp2' currently allowed. No filter will be applied")
      }
    }

    ## put into wide format
    mat_dt <-
      dcast(data = gtf,
            formula = circ_id ~ sample_id + Method,
            sep = name_sep,
            fill = 0,
            value.var = "BJR")

    ## return the matrix
    if (as_data_table) {

      mat_dt

    } else {

      as.matrix(data.frame(mat_dt, row.names = "circ_id"))

    }
  }
