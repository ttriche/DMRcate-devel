bridgeP <- 1e-2
bridgeSize <- 2000

depths <- GRanges(c("chr1","chr1","chr1","chr1","chr1","chr2","chr3","chr3"),
                  IRanges(start=c(1, 1999,7000, 9500, 12000, 3000, 4000, 6000),
                          end=c(500,2499,8000, 10500, 14000, 3300, 4500, 7000)),
                  score=c(-2, -3, 4, 4, -3, 2, 5, 1))

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
bridgeMe <- data.frame(qh=queryHits(d2n), sh=subjectHits(d2n))
anchors <- apply(bridgeMe, 1, function(x) depths[c(x[["qh"]], x[["sh"]])])
bridges <- do.call(c,
                   lapply(anchors, function(x) GRanges(unique(seqnames(x)),
                                                       IRanges(min(end(x)),
                                                               max(start(x))),
                                                       score=median(score(x)))))
sort(c(depths, bridges))

