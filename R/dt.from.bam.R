
### use ov.from.bam3 function to import a bamfile into a data.table

dt.from.bam <- function(bamfile, pattern='.bam', ...) {
  bam <- data.table(NULL)
  try({
    #bamfile <- bamfiles[i]
    bamname <- gsub('.*\\/', '', bamfile)
    bamname <- gsub(pattern, '', bamname)
    message('Start analyzing ', bamfile, '...')
    bam <- ov.from.bam3(bamfile, ...)
    bam[,sample := ..bamname]
  })
  
  bam <- bam[order(sample, seqnames, qname, -mapq, flag, cigar, strand)]
  bam[ , aln_ID  := .GRP, by = .(sample, seqnames, qname, cigar, flag, mapq, strand)]
  
  bam[ , aln_nr  := rank(aln_ID, ties.method="min"), by = qname]
  
  if(length(unique(bam[,aln_ID])) != length(unique(bam[,aln_ID]))) { message('something is wrong with the IDs !') } 
  
  return(bam)
}
