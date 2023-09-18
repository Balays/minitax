require(GenomicAlignments)

cigar.sum <- function(cigar) {
  cigar.df <- t(as.data.frame(explodeCigarOpLengths(cigar, ops=CIGAR_OPS)))
  colnames(cigar.df) <- unlist(explodeCigarOps(cigar, ops=CIGAR_OPS))
  rownames(cigar.df) <- NULL
  cigar.df <- as.data.frame(cigar.df)
}
