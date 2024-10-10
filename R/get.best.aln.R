#' Filter reads for the highest mapq, and/or filter out supplementary and/or secondary alignments.
#' Requires data.frame containing bam flag information
#' Can remove whole the read that has a supplementary and/or secondary alignment or that particular alignment  only.
#' @export
#'
get.best.aln <- function(bam, bam.flags, best.mapq=T, rm.supplementary=T, rm.secondary=T, rm.reads.or.alns='reads', keep.chim=F,
                         give.qname.ID.to.secondary=F) {

  require(tidyr)
  require(dplyr)
  require(stringr)
  require(data.table)
  
  bam.seq <- bam
  luniq.qname <- length(unique(bam.seq$qname))
  ### get best mapq
  if (best.mapq) {
    bam.mapq <- bam.seq %>% group_by(qname) %>% reframe(qname, mapq.max=max(mapq))
    #bam.mapq <- unique.data.frame(bam.mapq)
    bam.mapq <- bam.mapq %>% distinct(qname, mapq.max)
    bam.seq  <- as.data.frame(merge(
      data.table(bam.seq), 
      data.table(bam.mapq), 
      by.x = c('qname', 'mapq'), 
      by.y = c('qname', 'mapq.max')) )

    if (luniq.qname != length(unique(bam.seq$qname))) {
      message('When filtering for best mapq, ', luniq.qname - length(unique(bam.seq$qname)), ' reads were dropped!')
    } else {
      message('Reads were filtered, based on mapq for the best alignment only!')
    }

  }

  if (give.qname.ID.to.secondary) {
    sec.flag   <- 256
    sec.flags  <- c(sec.flag, bam.flags$Decimal + sec.flag)
    #bam.seq$aln_ID <- bam.seq$qname
    qname.sec <- bam.seq$qname[is.element(bam.seq$flag, c(sec.flags)) ]
    i <- '1:697|081205e4-8beb-4130-9730-8e51e239aa9c' #qname.sec[1]

    bam.sec <- bam.seq[is.element(bam.seq$qname, qname.sec), ]

    aln.count <- bam.sec %>%
      group_by(sample, qname) %>%
      summarise(cigar, start, end, strand, aln_count=seq_along(cigar))

    aln.count$aln_ID <- paste0(aln.count$qname, '_', aln.count$aln_count)

    bam.sec <- merge(bam.sec, aln.count, by=c('sample', 'qname', 'cigar', 'start', 'end', 'strand'))

    bam.nosec <- bam.seq[!is.element(bam.seq$qname, qname.sec), ]
    bam.nosec$aln_ID <- bam.nosec$qname

    bam.seq   <- plyr::rbind.fill(bam.nosec, bam.sec)
  }

  noprim <- NULL
  if (rm.secondary) {
    sec.flag   <- 256
    sec.flags  <- c(sec.flag, bam.flags$Decimal + sec.flag)

    if(rm.reads.or.alns =='reads') {
      qname.sec <- bam.seq$qname[is.element(bam.seq$flag, c(sec.flags)) ]
      bam.nosec   <- bam.seq[!is.element(bam.seq$qname, qname.sec), ]
    } else if(rm.reads.or.alns =='alns') {
      bam.nosec  <- bam.seq[!is.element(bam.seq$flag, c(sec.flags)), ]
    } else {stop('error rm.reads.or.alns must be either "reads" or "alns" ')}
    bam.seq    <- bam.nosec
    if (luniq.qname != length(unique(bam.seq$qname))) {
      noprim <- setdiff(bam$qname, bam.nosupp$qname)
      message('When filtering out secondary alignments, ', luniq.qname - length(unique(bam.seq$qname)), ' reads were dropped!')
    } else {

      message('Secondary alignments were fitered out!')
    }

  }

  nolinear <- NULL
  if (rm.supplementary) {
    supp.flag   <- 2048
    supp.flags  <- c(supp.flag, bam.flags$Decimal + supp.flag)

    if(rm.reads.or.alns =='reads') {
      qname.supp  <- bam.seq$qname[is.element(bam.seq$flag, c(supp.flags)) ]
      bam.nosupp  <- bam.seq[!is.element(bam.seq$qname, qname.supp), ]
    } else if(rm.reads.or.alns =='alns') {
      bam.nosupp  <- bam.seq[!is.element(bam.seq$flag, c(supp.flags)), ]
    } else {stop('error rm.reads.or.alns must be either "reads" or "alns" ')}
    bam.seq     <- bam.nosupp
    if (luniq.qname != length(unique(bam.seq$qname))) {
      nolinear <- setdiff(bam$qname, bam.nosupp$qname)
      message('When filtering out supplementary alignments, ', luniq.qname - length(unique(bam.seq$qname)), ' reads were dropped!')
    } else {

      message('Supplementary alignments were fitered out!')
    }
  }

  chim <- c(noprim, nolinear)
  if(!is.null(chim)) {
    bam.chim <- bam[is.element(bam$qname, chim), ]

    if (keep.chim) {
      message('Keeping the alignment, which has highest mapq and and greatest width')
      #qname <- '25bb9865-f2e5-44f8-b88e-8956af95107f'
      best.reads <- data.frame(NULL)
      for (qname in unique(bam.chim$qname)) {
        read.df <- bam[bam$qname == qname, ]
        read.df <- read.df[order(read.df$mapq, read.df$width, decreasing = T), ]
        read.df <- read.df[1,]
        best.reads <- plyr::rbind.fill(read.df, best.reads)
      }
      bam.chim <- best.reads
    } else {
      message('Keeping primary alignments in those cases where the secondary or supplementary had a higher mapq.')
      bam.chim <- bam.chim[bam.chim$flag == 0, ]
    }
    bam.seq <- plyr::rbind.fill(bam.chim, bam.seq)
  }
  return(bam.seq)
}


