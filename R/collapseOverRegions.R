#' alias to collapseAtDMRs
#'
#' @param x         the SummarizedExperiment-like thing to collapse
#' @param regions   the regions to collapse over, as a GRanges
#' @param how       the function used to collapse x over the regions (median)
#' @param impute    should missing values be imputed first? (TRUE) 
#' 
#' @return          an object not unlike x, but collapsed over the regions.
#'
#' @importFrom matrixStats colMedians 
#'
#' @export
collapseOverRegions <- function(x, regions, 
                                how=c("median", "mean", "sum", "max", "min"),
                                impute=TRUE) {
  how <- match.arg(how)
  collapseAtDMRs(x, regions, how, impute)
}
