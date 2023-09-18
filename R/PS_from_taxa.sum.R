
PS_from_taxa.sum <- function(
    taxa.sum, metadata=NULL,
    count.colname  ='count', ranks=NULL,
    cols_to_unique = c('tax.identity', ranks),
    taxID='species',
    round_counts_decimal=0
    ) {

  require(phyloseq)

  ### Taxonomy table ####
  taxtab   <- unique(taxa.sum[,..cols_to_unique])
  taxtab   <- data.frame(taxtab[,], row.names = unlist(taxtab[,..taxID]))


  ### Count (OTU) table ####

  if (!is.na(round_counts_decimal)) {
    taxa.sum[, (count.colname) := lapply(.SD, round, round_counts_decimal), .SDcols = count.colname]
  }

  cols_to_sum <- c(taxID, 'sample', count.colname)
  otutab    <- data.table(taxa.sum[,..cols_to_sum] %>% spread('sample',count.colname, fill=0))
  cnames    <- colnames(otutab)[-1]
  otutab    <- data.frame(otutab[,-1], row.names = unlist(otutab[,..taxID]))
  colnames(otutab) <- cnames

  if(length(samples) == 1 ) {colnames(otutab) <- samples}

  ### Metadata ####
  if(is.null(metadata)) {
    metadata <- data.frame(sample=unique(taxa.sum$sample), row.names = unique(taxa.sum$sample))
  }

  ### make PS object ####

  ps <- phyloseq(tax_table(as.matrix(taxtab)),
                 otu_table(as.matrix(otutab), taxa_are_rows = T),
                 sample_data(metadata))


  return(ps)
}
