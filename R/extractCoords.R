#' Extract coordinates (for transformation into GRanges; somewhat superfluous)
#' Thank you to Xavier Pastor from Bioconductor mailing list for this patch
#'
#' @param  xx    a set of strings in the format chr:start-end[:strand]
#' @return       a data.frame with columns chrom, chromStart, chromEnd
#' @export
extractCoords <- function(xx)
{
    coords <- sapply(xx, strsplit, '[:-]')
    coords <- as.data.frame(do.call(rbind, coords), stringsAsFactors=F)
    colnames(coords) <- c('chrom', 'chromStart', 'chromEnd')
    return(coords)
}
