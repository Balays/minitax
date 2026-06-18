
minitax2 <- function(minimap2, db='proGcontigs', db.data=prog.db, db.uni.data=NULL, chunk_n=1, calc.cigar.score=T,
                     ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
                     ) {
  
  require(tidyr)
  require(dplyr)
  require(data.table)
  
  dup       <- function(x) x[duplicated(x)]
  luniq     <- function (x) length(unique(x))
  luniq.df  <- function (x) apply(x, 2, luniq)
  try({ minimap2 <- setnames(minimap2, old='seqnames', new='rname')})
  minimap2 <- data.table::as.data.table(minimap2)
  
  validate_seqname_db <- function(db.uni.data, ranks) {
    if (is.null(db.uni.data)) {
      stop("db.uni.data is NULL. Custom and sequence-name based databases require a taxonomy table with seqnames, taxid, and rank columns.", call. = FALSE)
    }
    db.uni.data <- data.table::as.data.table(db.uni.data)
    cols <- c('seqnames', 'taxid', ranks)
    missing.cols <- setdiff(cols, colnames(db.uni.data))
    if (length(missing.cols) > 0) {
      stop("The database table for db='", db, "' is missing required column(s): ",
           paste(missing.cols, collapse = ", "), call. = FALSE)
    }
    db.uni.data
  }
  
  join_seqname_db <- function(minimap2, db.uni.data, ranks) {
    cols <- c('seqnames', 'taxid', ranks)
    db.uni.data <- validate_seqname_db(db.uni.data, ranks)
    
    minimap2.contigs <- minimap2[, c('qname', 'rname', 'flag', 'mapq', 'cigar')]
    minimap2.contigs[, rname := as.character(rname)]
    minimap2.contigs[, rname := sub("\\s.*$", "", rname)]
    
    db.uni.sub <- data.table::copy(db.uni.data[, ..cols])
    db.uni.sub[, seqnames := as.character(seqnames)]
    db.uni.sub[, taxid := as.character(taxid)]
    for (rank in ranks) {
      db.uni.sub[, (rank) := as.character(get(rank))]
    }
    data.table::setkeyv(db.uni.sub, 'seqnames')
    
    minimap2.taxa <- db.uni.sub[minimap2.contigs, on = .(seqnames = rname), nomatch = 0]
    if (nrow(minimap2.taxa) == 0) {
      example.rnames <- paste(utils::head(unique(minimap2.contigs$rname), 5), collapse = ", ")
      example.seqnames <- paste(utils::head(unique(db.uni.sub$seqnames), 5), collapse = ", ")
      stop("No alignments could be joined to the database taxonomy table for db='", db, "'. ",
           "Check that BAM reference names match db_data.tsv seqnames. ",
           "Example BAM rname(s): ", example.rnames, ". ",
           "Example db seqnames: ", example.seqnames, call. = FALSE)
    }
    data.table::setnames(minimap2.taxa, 'seqnames', 'rname')
    unique(minimap2.taxa, by = c("qname", "mapq", "taxid", ranks))
  }
  
  if (db=='proGcontigs_2') {
    if (is.null(db.uni.data)) {
      prog.db.uni <- db.data %>% distinct(across(all_of(c('taxid', ranks))))
    }
    minimap2.contigs <- data.table(taxid=as.numeric(gsub('\\..*', '', minimap2$rname)), minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')])
    minimap2.taxa    <- data.table::merge.data.table(minimap2.contigs, prog.db.uni, by='taxid', all=F)
    taxa.contigs     <- unique(minimap2.taxa, by=c("qname", "mapq", "taxid", ranks))
    
  } else if (db=='proGcontigs_3.host' | db=='proGcontigs_3.repres') {
    if (is.null(db.uni.data)) {
      prog.db.uni <- db.data %>% distinct(across(all_of(c('seqnames', ranks))))
    }
    minimap2.contigs <- data.table(taxid=as.numeric(gsub('\\..*', '', minimap2$rname)), minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')])
    minimap2.taxa    <- data.table::merge.data.table(minimap2.contigs, prog.db.uni, by.x='rname', by.y='seqnames', all=F)
    taxa.contigs     <- unique(minimap2.taxa, by=c("qname", "mapq", "taxid", ranks))
    
  } else if (db == 'rrnDB') {
    minimap2.contigs <- minimap2[,c('taxid', 'qname', 'rname', 'flag', 'mapq')]
    taxa.contigs     <- minimap2[,c("taxid", "qname", "mapq", 'taxid', ranks)]
    
  } else if (db == 'mouse.toy') {
    taxa.contigs <- unique.data.frame(merge(minimap2, db.data[,c('seq_id', 'taxid', ranks)], by.x='rname', by.y='seq_id'))
    
  } else if (db == 'EMUdb') {
    cols <- c('taxid', ranks)
    minimap2.contigs <- data.table(taxid=as.numeric(gsub(':.*', '', minimap2$rname)), minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')])
    taxa.contigs     <- unique(merge(minimap2.contigs[,c("taxid", "qname", "mapq", 'cigar')], db.uni.data[,..cols], by='taxid'))
    
  } else {
    ## Generic sequence-name database path.
    ## This covers all_NCBI_genomes, GTDB_SSU, ncbi_refseq_16S, and any custom db_data.tsv
    ## with columns: seqnames, taxid, and the configured rank columns.
    taxa.contigs <- join_seqname_db(minimap2, db.uni.data, ranks)
  }
  
  rank.df <- data.frame(rank=ranks, rank.level=c(1:length(ranks)))
  tax.df <- data.table(taxa.contigs)
  tax.df[tax.df == ' '] <- NA
  
  rm(taxa.contigs, minimap2.contigs)
  
  rank_cols <- rev(ranks)
  lowest_common_rank <- function(x) {
    for (rank in rank_cols) {
      unique_taxa <- unique(x[[rank]], na.rm = TRUE)
      if (length(unique_taxa) == 1 && !is.na(unique_taxa)) {
        return(data.table(rank = rank, taxon = unique_taxa))
      }
    }
    return(data.table(rank = '', taxon = ''))
  }
  
  taxa <- tax.df
  result <- taxa[, lowest_common_rank(.SD), by = qname]
  result[result == ''] <- NA
  colnames(result)[-1] <- c('tax.identity.level', 'tax.identity')
  taxa <- merge(taxa, result, by='qname')
  return(taxa)
}
