
minitax2 <- function(minimap2, db='proGcontigs',  db.uni.data=NULL, 
                     ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
                     ) {
  
  require(tidyr)
  require(dplyr)
  require(data.table)
  
  dup       <- function (x) x[duplicated(x)]
  luniq     <- function (x) length(unique(x))
  luniq.df  <- function (x) apply(x, 2, luniq)

  try({ minimap2 <- setnames(minimap2, old='seqnames', new='rname')})
  
         if (db=='proGcontigs_2') {
    #db.data <- prog.db
    if (is.null(db.uni.data)) {
      prog.db.uni      <- db.data %>% distinct(across(all_of(c('taxid', ranks)))) } #unique.data.frame(db.data[,c('taxid', ranks)]) 
    minimap2.contigs <- data.table(taxid=as.numeric(gsub('\\..*', '', minimap2$rname)), minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')])
    minimap2.taxa    <- data.table::merge.data.table(minimap2.contigs, prog.db.uni, by='taxid', all=F) #type = "inner")
    taxa.contigs     <- unique(minimap2.taxa, by=c( "qname", "mapq","taxid", ranks))
    
  } else if (db=='proGcontigs_3.host' | db=='proGcontigs_3.repres') {
    #db.data <- prog.db
    if (is.null(db.uni.data)) {
      prog.db.uni      <- db.data %>% distinct(across(all_of(c('seqnames', ranks)))) } #unique.data.frame(db.data[,c('taxid', ranks)]) 
    minimap2.contigs <- data.table(taxid=as.numeric(gsub('\\..*', '', minimap2$rname)), minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')])
    #minimap2.contigs <- minimap2
    minimap2.taxa    <- data.table::merge.data.table(minimap2.contigs, prog.db.uni, by.x='rname', by.y='seqnames', all=F) #type = "inner")
    taxa.contigs     <- unique(minimap2.taxa, by=c( "qname", "mapq","taxid", ranks))
    
  } else if (db == 'rrnDB') {
    
    minimap2.contigs <- minimap2[,c('taxid', 'qname', 'rname', 'flag', 'mapq')]
    taxa.contigs     <- minimap2[,c("taxid", "qname", "mapq", 'taxid', ranks)]
    
  } else if (db == 'mouse.toy') {
    #db.data <- mgtoy.idx
    taxa.contigs <- unique.data.frame(merge(minimap2, db.data[,c('seq_id', 'taxid', ranks)], by.x='rname', by.y='seq_id'))
    
  } else if (db == 'EMUdb') {
    #db.data <- emu.db
    cols <- c('taxid', ranks)
    minimap2.contigs <- data.table(taxid=as.numeric(gsub(':.*', '', minimap2$rname)), minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')])
    taxa.contigs     <- unique(merge(minimap2.contigs[,c("taxid", "qname", "mapq", 'cigar')], db.uni.data[,..cols], by='taxid'))
    
  } else if (db == 'all_NCBI_genomes') {
    
    cols <- c('seqnames', 'taxid', ranks)
    minimap2.contigs <- minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')]
    minimap2.taxa    <- merge(minimap2.contigs, db.uni.data[,..cols], by.x='rname', by.y='seqnames', all=F)
    taxa.contigs     <- unique(minimap2.taxa, by=c( "qname", "mapq","taxid", ranks))
    
  } else {
    #stop('db must be one of the following: ProGenomes; rrnDB; EMUdb or mouse.toy !')
    
    message('using: ', db, ' database...')
    
    minimap2.contigs <- minimap2[,c('qname', 'rname', 'flag', 'mapq', 'cigar')]
    
    cols  <- unique(c('seqnames', ranks, 'taxid'))
    minimap2.taxa    <- merge(minimap2.contigs, db.uni.data[,..cols], by.x='rname', by.y='seqnames', all=F)
    
    cols  <- c( "qname", "mapq", unique(c('taxid', ranks)) )
    taxa.contigs     <- unique(minimap2.taxa, by=cols)
    
  }
  
  message('The following taxon levels will be used in the minitax step: \n', paste0(ranks, collapse='; '))
  
  rank.df <- data.frame(rank=ranks, rank.level=c(1:length(ranks)))
  
  tax.df <- data.table(taxa.contigs)
  tax.df[tax.df == ' '] <- NA
  
  rm(taxa.contigs, minimap2.contigs)
  
  # Reverse the column order since we want to start from the lowest taxonomic rank
  rank_cols <- rev(ranks)
  
  # Function to find the lowest common rank for each qname
  lowest_common_rank <- function(x) {
    for (rank in rank_cols) {
      unique_taxa <- unique(x[[rank]], na.rm = TRUE)
      if (length(unique_taxa) == 1 && !is.na(unique_taxa)) {
        return(data.table(rank = rank, taxon = unique_taxa))
      }
    }
    return(data.table(rank = '', taxon = ''))
  }
  
  taxa   <- tax.df[, lowest_common_rank(.SD), by = qname]
  taxa[taxa == ''] <- NA
  
  colnames(taxa)[-1] <- c('tax.identity.level', 'tax.identity')
  taxa <- merge(tax.df, taxa, by='qname')
    
  try({ taxa <- setnames(taxa, old='rname', new='seqnames')})
  
  #saveRDS(max_cigar_score, paste0(file, ".rds"))
  return(taxa)
}

