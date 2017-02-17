#' this is the recommended user-facing probe-name version
#' 
#' @import minfi
#' 
#' @param grset         A GenomicRatioSet or something very much like it 
#' @param cutoff        Maximum absolute M-value before truncation (10)
#' @param returnBetas   Return beta values instead of M-values (FALSE) 
#' 
#' @return              A matrix of tidied and truncated M- or beta-values
#'
#' @export
prepMvals <- function(grset, cutoff=10, returnBetas=FALSE) { 
  xx <- prepM(getSafeProbes(getM(grset)))
  if (returnBetas) return(ilogit2(xx))
  else return(xx)
}
