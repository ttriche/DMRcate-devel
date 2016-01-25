#' function to allow "piling up" differentially methylated region significance
#'
#' @param output        DMRcate output (of class "dmrcate.output" as it happens)
#' @param stepSize      multiplier for scores. Typically if using 10**-y, it's 1
#' @param bridgeSize    bridge regions between significant DMRs up to this size
#' @param bridgeP       don't bridge gaps with adjusted p-values less than this
#'
#' @export
getDMRdepth <- function(output, stepSize=1, bridgeSize=2000, bridgeP=1e-2, ...){

  stopifnot(class(output) == "dmrcate.output")

  bySign <- extractRanges(output, bySign=TRUE)

  binDepths <- function(x) {
    bins <- disjoin(x)
    bins$score <- countOverlaps(bins, x)
    if(median(sign(score(x)) < 0)) bins$score <- -1 * bins$score
    return(bins)
  }

  ## merge hyper and hypo, should never overlap (!)
  data(seqinfo.hg19)
  depths <- GRangesList(lapply(bySign, binDepths))
  ## add seqinfo so it can be dumped to a bigWig file
  seqinfo(depths) <- seqinfo.hg19[seqlevels(depths)] 

  fixOverlaps <- function(overlaps, segs) { # {{{
    segs$score <- rep(0, length(segs))
    for (i in seq_along(segs)) {
      segs[i]$score <- sum(subsetByOverlaps(overlaps, segs[i])$score)
    }
    return(segs[score(segs) != 0])
  } # }}}

  ## patch up any overlapping DMRs 
  overlapping <- NULL
  if(length(depths) > 1) {
    overlapping <- GRangesList(hyper=subsetByOverlaps(depths$hyper,depths$hypo),
                               hypo=subsetByOverlaps(depths$hypo,depths$hyper))
  } ## only makes sense if we have both kinds!
  if (length(unlist(overlapping)) > 0) {
    message("Overlapping DMRs with opposite signs found, fixing...")
    excludeDisputed <- function(x, y) x[-queryHits(findOverlaps(x, y))]
    hyper <- excludeDisputed(depths$hyper, depths$hypo)
    hypo <- excludeDisputed(depths$hypo, depths$hyper)
    depths$hyper <- hyper
    depths$hypo <- hypo 
    overlaps <- sort(unlist(overlapping))
    overlaps <- split(overlaps, seqnames(overlaps))
    overlaps <- overlaps[lapply(overlaps, length) > 0]
    segs <- lapply(overlaps, disjoin)
    fixed <- unlist(GRangesList(mapply(fixOverlaps, overlaps, segs)))
    depths <- sort(c(unlist(depths), fixed))
  } else { 
    depths <- sort(unlist(depths))
  }
  depths$score <- depths$score * stepSize

  if (bridgeSize > 0) {
    message("Need more testing of gap-bridging code...")
    consider <- which(abs(depths$score) >= abs(log10(bridgeP)))
    d2n <- distanceToNearest(depths)
    keepQ <- which(queryHits(d2n) %in% consider)
    keepS <- which(subjectHits(d2n) %in% consider) 
    d2n <- d2n[ intersect(keepQ, keepS) ]
    d2n <- d2n[ which(mcols(d2n)$distance <= bridgeSize) ]
    mcols(d2n)$sign1 <- sign(depths$score[queryHits(d2n)])
    mcols(d2n)$sign2 <- sign(depths$score[subjectHits(d2n)])
    d2n <- d2n[which(mcols(d2n)$sign1 == mcols(d2n)$sign2)]
    qh <- queryHits(d2n)
    sh <- subjectHits(d2n)
    dupe <- rep(FALSE, length(qh))
    for (i in seq_along(sh)) if (i > 1 && qh[i - 1] == sh[i]) dupe[i] <- TRUE
    d2n <- d2n[!dupe]
    if (length(d2n) > 0) {
      message("Bridging ", length(d2n), " eligible sets of DMRs...")
      bridgeMe <- data.frame(qh=queryHits(d2n), sh=subjectHits(d2n))
      anchors <- apply(bridgeMe, 1, function(x) depths[c(x[["qh"]], x[["sh"]])])
      bridges <- do.call(c,
                       lapply(anchors, function(x) 
                              GRanges(unique(seqnames(x)),
                                      IRanges(min(end(x))+ 1,
                                              max(start(x)) - 1),
                                      score=median(score(x)))))
      depths <- sort(c(depths, bridges))
      message("...done.")
    }
  }

  return(depths)

}
