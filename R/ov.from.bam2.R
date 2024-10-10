#' Import alignments from a .BAM file.
#' Optional filtering can be applied for the alignments, based on bam flags, mapping quality, and reference sequences.
#' LoRTIA ouptuts can be imported, along with LoRTIA tags for further filtering.
#' Gaps (Ns) in the alignments can be removed into distinct alignments for the same read using the CIGAR values.
#'
#' @export

ov.from.bam2 <- function (bamfile, flag.tokeep = c(0, 16), flag.tocrop = c(4), crop.na.cigar = T,
          seqnames.tofilt = NA, mapq.filt = NA, rm.gaps.in.aln = F, add.primes = T,
          what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand"),
          flag = scanBamFlag(),
          lortia.tags = c("l3", "l5", "r3", "r5"), is.lortia = F, rm.false.exons = F, rm.non.correct = F) {
  require(GenomicAlignments)
  require(Rsamtools)
  #bamfile <- 'I:/data/SARS-CoV2/mapped_v8/.bam.fastq.pych.bam/Hpi10_A.merged_pychopped.bam'

  ## Is this LoRTIA output?
  if (is.lortia) {
    params <- ScanBamParam(what = what, tag = lortia.tags)
    bam <- as.data.frame(scanBam(bamfile, param = params))
    bam$tags <- paste(bam[1, paste0("tag.", lortia.tags)], collapse = ",")
    bam$tags <- apply(as.data.frame(
      apply(bam[, paste0("tag.", lortia.tags)], 2, function(x) gsub(".*,", "", x))),
                      1, function(x) paste0(unique(x), collapse = ","))
    
    ## Remove alignments designated as having no correct adapter by LoRTIA
    if (rm.non.correct) {
      bam <- bam[grepl("correct", bam$tags), ]
    }
    ## Remove alignments designated as having false exons by LoRTIA
    if (rm.false.exons) {
      bam <- bam[!grepl("false", bam$tags), ]
    }
  } else {
    params <- ScanBamParam(what = what)
    bam <- as.data.frame(scanBam(bamfile, param = params))
  }
  
  ## Filter based on seqnames
  if (all(!is.na(seqnames.tofilt))) {
    bam <- bam[!is.na(bam$rname), ]
    bam <- bam[is.element(bam$rname, seqnames.tofilt), ]
  }
  
  ## Filter based on mapping quality
  if (all(!is.na(mapq.filt))) {
    bam <- bam[bam$mapq >= mapq.filt, ]
  }
  
  ## Filter based on flags
  if (all(!is.na(flag.tokeep))) {
    bam <- bam[is.element(bam$flag, flag.tokeep), ]
  }
  if (all(!is.na(flag.tocrop))) {
    bam <- bam[!is.element(bam$flag, flag.tocrop), ]
  }
  
  if (crop.na.cigar) {
    bam <- bam[!is.na(bam$cigar), ]
  }
  
  ## separate unmapped reads
  bam.na <- bam[is.na(bam$rname),  ]
  bam    <- bam[!is.na(bam$rname), ]
  
  ## Separate alignments with gaps into sub_alignments ??
  if (rm.gaps.in.aln) {
    bam$group <- seq(1, nrow(bam))
    aln  <- as.data.frame(extractAlignmentRangesOnReference(bam$cigar, pos = bam$pos, f = NULL))
    stopifnot(luniq(aln$group) == nrow(bam))
  } else {
    bam$group <- seq(1, nrow(bam))
    aln <- as.data.frame(cigarRangesAlongReferenceSpace(bam$cigar, pos = bam$pos, f = NULL, reduce.ranges=T))
    stopifnot(nrow(aln) == nrow(bam))
  }
  aln <- dplyr::select(aln, -group_name)
  aln <- merge(aln, bam, by = "group")[, -1]
  bam <- aln
  
   
  ## put back unmapped reads
  bam <- plyr::rbind.fill(bam, bam.na)
  
  ##
  if (is.lortia) {
    bam <- bam[, c(what, "start", "end"
                 , paste0("tag.", lortia.tags), "tags")]
  } else {
    bam <- bam[, c(what, "start", "end")]
  }
  colnames(bam)[1] <- "seqnames"
  
  ## add 3' and 5' end positions
  if (add.primes) {
    bam$prime5[!is.na(bam$strand)] <- bam$start[!is.na(bam$strand)]
    bam$prime3[!is.na(bam$strand)] <- bam$end[!is.na(bam$strand)]
    bam$prime5[bam$strand == "-" & !is.na(bam$strand)] <- bam$end[bam$strand   == "-" & !is.na(bam$strand)]
    bam$prime3[bam$strand == "-" & !is.na(bam$strand)] <- bam$start[bam$strand == "-" & !is.na(bam$strand)]
  }
 
  return(bam)
}
