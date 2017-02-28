#' extract disjoint GRanges with methylation fraction (beta) diffs as scores
#' the result can be usefully exported as a bed, bedGraph, or bigWig file 
#' 
#' @param   dmrcated  an object of class dmrcate.output
#' @param   minDiff   get rid of DMRs with less than this beta difference (0.1)
#' @param   bySign    split the results into a GRangesList by sign? (FALSE)
#' @param   withDMLs  return both DMRs and DMLs? (Forces bySign to FALSE)
#' @return  a GRanges or GRangesList, depending values of bySign and withDMLs
#' @export 
getDMRbetas <- function(dmrcated, minDiff=.1, bySign=FALSE, withDMLs=FALSE,...){
  if (is(dmrcated, "GenomicRangesORGRangesList")) {
    return(dmrcated)
  } else if (!is(dmrcated, "dmrcate.output")) {
    stop("Argument is not of class dmrcate.output; cannot process.")
  } else { 
    data(seqinfo.hg19)
    gr <- as(dmrcated$results$hg19coord, "GRanges") 
    DMRs <- disjoin(gr)
    names(DMRs) <- as(DMRs, "character")
    probes <- with(dmrcated$input,
                   GRanges(seqnames=CHR,
                           IRanges(start=pos, 
                                   end=pos),
                           name=ID,
                           score=betafc))
    names(probes) <- probes$name
    probes <- subsetByOverlaps(probes, DMRs)
    probes$DMR <- names(DMRs)[queryHits(findOverlaps(DMRs, probes))]
    DMRs$name <- names(DMRs) 
    DMRs$score <- -1*sapply(split(probes$score,probes$DMR), median)[names(DMRs)]
    if (!is.null(minDiff)) DMRs <- subset(DMRs, abs(score) >= minDiff)
    seqinfo(DMRs) <- seqinfo.hg19[seqlevels(DMRs)] 
    DMRs <- sort(DMRs) # proper order 
    if (withDMLs == TRUE) { 
      seqinfo(probes) <- seqinfo.hg19[seqlevels(probes)] 
      probes <- sort(probes) # proper order 
      return(GRangesList(DMRs=DMRs, DMLs=probes))
    } else if (bySign == FALSE) {
      return(DMRs)
    } else if (bySign == TRUE) {
      return(split(gr, ifelse(sign(score(gr)) > 0, "hyper", "hypo")))
    }
  }
}
