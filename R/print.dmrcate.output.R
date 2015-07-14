#' make it harder for me to shoot myself in the foot calling `some_DMR_results`
#'  
#' @param an object of class dmrcate.output
#'
#' @import scales
#' 
#' @export
print.dmrcate.output <- function(object) { 
  cat("Object of class", class(object), "\n\n")
  cat("  Input data:", nrow(object$input), "features\n")
  cat("     Results:", nrow(object$results), "DMRs\n")
  cat("     Cutoffs: p <", paste(scientific(object$cutoff), collapse=", "),"\n")
}
