#' get "safe" probes based on EPIC remappings
#' 
#' see Zhou et al., https://doi.org/10.1093/nar/gkw967
#' 
#' @param object  a matrix of M or Beta values
#' @param rmXY    remove chrX and chrY? (FALSE) 
#
#' @return        the same matrix minus "unsafe" probes 
#'
#' @export 
getSafeProbes <- function(object, rmXY=FALSE, ...) {

  data(safeProbes)
  stopifnot(is.matrix(object))
  keep <- intersect(rownames(object), names(safeProbes))
  if(rmXY) {
    keep <- setdiff(keep, 
                    names(subset(safeProbes, seqnames %in% c("chrX","chrY"))))
  }
  return(object[keep, ])

}
