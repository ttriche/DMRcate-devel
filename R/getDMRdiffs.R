#' Create GRanges of imputed beta FC between samples to dump to a BigWig 
#'
#' @param output  DMRcate output (of class "dmrcate.output", as it happens)
#' 
#' @return summarized unscaled effect size estimates for the output 
#'
#' @import GenomicRanges
#' 
#' @export 
getDMRdiffs <- function(output, ...) { 

  data(seqinfo.hg19)
  stopifnot(class(output) == "dmrcate.output")
  betaFcGr <- makeGRangesFromDataFrame(output$input, 
                                       start.field="pos", 
                                       end.field="pos", 
                                       seqnames.field="CHR", 
                                       keep.extra.columns=TRUE)
  betaFcGr$score <- betaFcGr$betafc
  seqinfo(betaFcGr) <- seqinfo.hg19[seqlevels(betaFcGr)]
  return(betaFcGr)

}
