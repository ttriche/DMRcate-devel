#' pretty much what it says, except that (ironically) it calls ComplexHeatmap
#'
#' @param x   a matrix, GenomicRatioSet, or other rectangular object
#' @param k   choose the top [this many] features by SD
#'
#' @return invisibly, the same thing as Heatmap(...)
#'
#' @import ComplexHeatmap
#'
#' @export
getSimpleMethHeatMap <- function(x, k=100, asSNPs=F, binary=F, rotate=F, 
                                 clustering_method_rows="ward",
                                 clustering_method_columns="ward",
                                 clustering_distance_rows="binary",
                                 clustering_distance_columns="binary",
                                 ...) {

  if(is(x, "SummarizedExperiment") || is(x, "RangedSummarizedExperiment")){ #{{{
    x <- keepSeqlevels(x, paste0("chr", 1:22))
    x <- x[grep("^rs", rownames(x), invert=TRUE), ] 
    x <- x[grep("^ch[1234567890]", rownames(x), invert=TRUE), ] 
    if(all(grepl("^cg", rownames(x)))) {
      xxx <- rmSNPandCH(getBeta(x), mafcut=0.01)
      sds <- rowSds(xxx, na.rm=TRUE)
      names(sds) <- rownames(xxx)
      rowordering <- names(sds)[order(sds, decreasing=TRUE)]
      xx <- xxx[head(rowordering, k), ] 
    } else { 
      sds <- rowSds(getBeta(x), na.rm=TRUE)
      names(sds) <- rownames(x)
      rowordering <- names(sds)[order(sds, decreasing=TRUE)]
      xx <- getBeta(x)[head(rowordering, k), ]
    }
    if("indicator" %in% names(colData(x)) && is.null(ColSideColors)) {
      ColSideColors <- x$indicator
    } # }}}
  } else { # {{{
    sds <- rowSds(data.matrix(x))
    rowordering <- order(sds, decreasing=TRUE)
    xx <- x[head(rowordering, k), ]
  } # }}}
  if(asSNPs != TRUE) { # {{{
    cfun <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                               "yellow", "#FF7F00", "red", "#7F0000")) # }}}
  } else { # {{{ as SNPs 
    xx <- round(xx * 2)
    cfun <- colorRampPalette(c("blue","yellow","red"))
  } # }}}

  if (rotate == TRUE && !is.null(ColSideColors)) { # {{{
    hfun <- function(X) {
      Heatmap(X, col=cfun(255),
              clustering_method_rows=clustering_method_rows,
              clustering_method_columns=clustering_method_columns,
              clustering_distance_rows=clustering_distance_rows,
              clustering_distance_columns=clustering_distance_columns,
              ...) 
    } # }}}
  } else if (rotate != TRUE && !is.null(ColSideColors)) { # {{{
    hfun <- function(X) {
      Heatmap(X, col=cfun(255),
              clustering_method_rows=clustering_method_rows,
              clustering_method_columns=clustering_method_columns,
              clustering_distance_rows=clustering_distance_rows,
              clustering_distance_columns=clustering_distance_columns,
              ...) 
    } # }}}
  } else { # {{{
    hfun <- function(X, ...) {
      Heatmap(X, col=cfun(255),
              clustering_method_rows=clustering_method_rows,
              clustering_method_columns=clustering_method_columns,
              clustering_distance_rows=clustering_distance_rows,
              clustering_distance_columns=clustering_distance_columns,
              ...)
    } # }}}
  } # }}}

  X <- xx 
  if (rotate) X <- t(X)
  hfun(X, ...)

}
