asGR <- function(x) { 
  stopifnot(class(x) == "dmrcate.output")
  gr <- as(x$results$hg19coord, "GRanges")
  score(gr) <- x$results$maxbetafc
  gr$name <- x$results$gene_assoc
  return(gr)
}
