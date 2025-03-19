
### use ov.from.bam3 function to import a bamfile into a data.table

dt.from.bam <- function(bamfile, pattern='.bam', bamname=NULL, ...) {
  
  message('Start analyzing ', bamfile, '...')
  
  bam <- data.table(NULL)
  
  ## Sample name
  if (is.null(bamname)) {
    bamname <- gsub('.*\\/', '', bamfile)
    bamname <- gsub(pattern, '', bamname)
  }
  
  ## Import and add sample column
  try({
    #bamfile <- bamfiles[i]
    bam <- ov.from.bam3(bamfile, ...)
    bam[,sample := ..bamname]
  })
  
  try({
    ## Order
    bam <- bam[order(sample, seqnames, qname, -mapq, flag, cigar, strand)]
    
    ## Add IDs to each alignment
    bam[ , aln_ID  := .GRP, by = .(sample, seqnames, qname, cigar, flag, mapq, strand)]
    
    ## Add alignment number
    bam[ , aln_nr  := rank(aln_ID, ties.method="min"), by = qname]
    
    ## Add exon number
    bam[ , ex_nr   := fifelse(strand == '+', 1:n_distinct(start), 
                                  n_distinct(end):1), by = aln_ID]
    
    ## check
    if(length(unique(bam[,aln_ID])) != length(unique(bam[,aln_ID]))) { message('something is wrong with the IDs !') }
    
  })
  
  return(bam)
}
