#' extract GRanges corresponding to DMRs from the results of a dmrcate run
#' (no longer uses extractCoords since we can do as(character, "GRanges"))
#' 
#' @param   dmrcated  an object of class dmrcate.output
#' @param   bySign    boolean, shall the DMRs be extracted by sign?  (FALSE)
#' @param   atCutoff  a p-value cutoff desired for extraction, or NULL (NULL)
#' 
extractRanges <- function(dmrcated, bySign=FALSE, atCutoff=NULL, ...) {
  if (is(dmrcated, "GenomicRangesORGRangesList")) {
    return(dmrcated)
  } else if (!is(dmrcated, "dmrcate.output")) {
    stop("Argument is not of class dmrcate.output; cannot process.")
  } else {
    gr <- as(dmrcated$results$hg19coord, "GRanges") 
    mcols(gr) <- dmrcated$results
    gr$score <- gr$maxbetafc 
    gr$name <- names(gr) <- as(gr, "character") 
    if (!is.null(atCutoff)) gr <- subset(gr, gr$pcutoff == atCutoff)
    genome(gr) <- "hg19" ## nothing else is supported, yet :-/
    if (bySign == FALSE) {
      # GRanges
      return(gr)
    } else if (bySign == TRUE) {
      # GRangesList
      split(gr, ifelse(sign(score(gr)) > 0, "hyper", "hypo"))
    }
  }
}
