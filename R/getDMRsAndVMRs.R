#' omnibus DMR finder starting from a m-values and a design/Surv matrix
#' 
#' @param x           matrix of M-values or SummarizedExperiment-derived object
#' @param design      design matrix (for limma) or Surv object/matrix (for Cox)
#' @param contrasts   apply contrasts to e.g. paired samples? (only for limma)
#' @param cont.matrix if contrasts=T, supply a contrast matrix here
#' @param coef        column of the design matrix to fit DMRs (or "all")?
#' @param pcutoff     p-value cutoff; if not specified, step from 10**-1:10**-8
#' @param betacutoff  DMRs must have at least this great a maxbetaFC difference
#' @param p.adjust.method   how to control the FDR, default is limma-style
#' @param what        (not used yet) whether to use limma or Cox PH to tag DMRs
#'
#' @return            an object of type dmrcate.output 
#' 
#' @export
getDMRs <- function(x, design=NULL, contrasts=FALSE, cont.matrix=NULL, 
                    coef=2, pcutoff=NULL, betacutoff=0.1, keep.raw.DMLs=TRUE, 
                    p.adjust.method="limma", what=c("limma", "cox"), ...) {

  ## support for Cox PH is coming:
  what <- tolower(match.arg(what))

  if (what == "cox") {
    stop("Cox regression DMRs are not (yet) supported.")
  }
  if (what == "cox" && p.adjust.method == "limma") {
    stop("Cannot use limma-style p-value correction for Cox regression DRMs.")
  }

  if (is.null(pcutoff)) pcutoff <- 10 ** (-1 * seq(1, 8)) ## nice for bigWigs

  if (is(x, "SummarizedExperiment") || is(x, "RangedSummarizedExperiment")) {
    x <- prepMvals(x) ## would rather internalize its rowRanges within dmrcate
  }

  if (is.null(design)) {
    if ((is(x,"SummarizedExperiment") || is(x,"RangedSummarizedExperiment")) &&
        ("design" %in% names(exptData(x)) || "design" %in% names(metadata(x)))){
      if (is(x, "SummarizedExperiment")) design <- exptData(x)$design
      if (is(x, "RangedSummarizedExperiment")) design <- metadata(x)$design
    } else {
      stop("You need a design matrix (perhaps metadata(x)$design) to call DMRs")
    }
  }

  # recurse, if requested to fit all coefs... FIXME: write fast cpg.annotate()
  if (coef == "all") { # hack 
    # remember to add: cpgannot <- annotateForMultiDMRs(x, design=design)
    lapply(seq_len(ncol(design))[-1],
           function(y) getDMRs(x, design=design, contrasts=contrasts,
                               cont.matrix=cont.matrix, coef=y, pcutoff=pcutoff,
                               betacutoff=betacutoff, # cpgannot=cpgannot,
                               p.adjust.method=p.adjust.method, what=what, 
                               ...))
    # return DMRs for all columns
  } else { 
    message("Annotating individual CpGs...")
    DMRannot <- switch(what, 
                       limma=cpg.annotate(x, design=design, coef=coef),
                       cox=cox.annotate(x, Surv(x$OS, x$OSevent)))

    message("Demarcating significant regions...")
    res <- dmrcate(DMRannot, pcutoff=pcutoff, betacutoff=betacutoff, 
                   p.adjust.method=p.adjust.method, ...)
    if (keep.raw.DMLs) res$DMLs <- DMRannot
    return(res) 
  }
}

#' get variably methylated regions 
#' 
#' @param x           matrix of M-values or SummarizedExperiment-derived object
#' @param pcutoff     p-value cutoff; if not specified, step from 10**-1:10**-8
#'
#' @return            an object of type dmrcate.output 
#' 
#' @export 
getVMRs <- function(x, pcutoff=0.1, ...) 
{
  if (is(x, "SummarizedExperiment") || is(x, "RangedSummarizedExperiment")) {
    x <- prepMvals(x)
  }
  VMRannot <- cpg.annotate(x, analysis.type="variability", pcutoff=pcutoff)
  dmrcate(VMRannot, ...)
}

#' get both DMRs and VMRs from the same dataset 
#' 
#' @param x           matrix of M-values or SummarizedExperiment-derived object
#' @param design      design matrix (for limma) or Surv object/matrix (for Cox)
#' @param contrasts   apply contrasts to e.g. paired samples? (only for limma)
#' @param cont.matrix if contrasts=T, supply a contrast matrix here
#' @param coef        column of the design matrix to fit DMRs (or "all")?
#'
#' @return            an object of type dmrcate.output 
#' 
#' @export 
getDMRsAndVMRs <- function(x, design=NULL, contrasts=F, 
                           cont.matrix=NULL, coef=2, ...){
  mvals.x <- prepMvals(x) # don't impute twice (or more)
  list(DMRs=getDMRs(mvals.x, design, coef, ...), VMRs=getVMRs(mvals.x, ...))
}
