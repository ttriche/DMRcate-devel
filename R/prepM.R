#' this is the core of prepMvals but omits sanity checking; useful if collapsed
#'
#' @import impute 
#' 
#' @param M       M-value matrix
#' @param cutoff  extremes at which to truncate M-values 
#' 
#' @return        a matrix 
#'
#' @export
prepM <- function(M, cutoff=10) {
 
  if (is(M, "RangedSummarizedExperiment")) {
    if (is(M, "GenomicRatioSet") | is(M, "GenomicMethylSet")) {
      M <- getM(M) 
    } else {
      M <- logit2(assays(M)$Beta) # just in case
    }
  } else if (is(M, "SummarizedExperiment")) {
    M <- -1 * assays(M)$exprs ## HELP data 
  }

  ## get rid of +/-Inf cells
  M[ which(M > cutoff) ] <- cutoff 
  M[ which(M < -1 * cutoff) ] <- -1 * cutoff

  ## impute NAs via k-NN
  if (any(is.na(M))) {
    message("Imputing NAs...")
    set.seed(1234)
    M <- impute.knn(M)$data
  }
  return(M)

}
