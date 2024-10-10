


# Function to calculate and assign the most probable taxon for a given rank
assign_most_probable <- function(data, rank, thresholds) {
  
  # Define the taxonomic ranks in order
  rank_df <- data.frame(rank=ranks, rank_level = 1:length(ranks))
  
  # Calculate the total number of alignments for each qname
  total_alignments <- data[, .N, by = .(sample, qname, tax.identity, tax.identity.level, tax.identity_rank_level) ]
  setnames(total_alignments, 'N', 'N_alignments')
  
  # Calculate the number of alignments for each taxon at the given rank within each qname
  cols_by <- c('qname', rank)
  taxon_alignments <- data[, .(taxon_count = .N), by = cols_by]
  
  # Merge the total alignments with taxon alignments
  taxon_alignments <- merge(taxon_alignments, total_alignments, by = "qname")
  
  # Calculate the probability for each taxon at the given rank
  taxon_alignments[, taxon_prob := taxon_count / N_alignments]
  
  # Assign the most probable taxon to each qname
  taxon_alignments[, taxon_most_probable := .SD[which.max(taxon_prob)], by = qname]
  
  # Assing rank level
  taxon_alignments[,taxon_most_probable_rank := rank]
  
  taxon_alignments <- merge(taxon_alignments, rank_df, by.x='taxon_most_probable_rank', by.y='rank')
  setnames(taxon_alignments, 'rank_level', "taxon_most_probable_rank_level")
  
  # Check if the most probable taxon is a better hit than the tax identity
  rank_threshold <- thresholds$threshold[thresholds$rank == rank]
  taxon_alignments[,is_taxon_most_prob := fifelse(tax.identity_rank_level < taxon_most_probable_rank_level &
                                                    taxon_prob > rank_threshold,
                                                  T, F)]
  
  # Rename the columns appropriately
  colnames(taxon_alignments)[-1] <- gsub('taxon_', paste0(rank, '_'), colnames(taxon_alignments)[-1])
  
  cols_by <- c('qname', 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level',
               grep(rank, colnames(taxon_alignments), value = T))
  return(taxon_alignments[, ..cols_by])
  #return(taxon_alignments[,])
}


# Function to update the tax.identity based on probability and level
update_tax_identity <- function(data, rank) {
  
  
  prob_col <- paste0(rank, "_prob")
  most_probable_col <- paste0(rank, "_most_probable")
  most_probable_rank_level_col <- paste0(rank, "_most_probable_rank_level")
  is_taxon_most_prob_col <- paste0('is_', rank, "_most_prob")
  
  cols <- c('sample', 'qname', 'aln_nr', 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level', rank, 
            prob_col, most_probable_col, most_probable_rank_level_col, is_taxon_most_prob_col)
  
  data.subs <- data[,..cols]
  
  setnames(data.subs, names(data.subs)[-c(1:6)], gsub(rank, 'taxa', names(data.subs)[-c(1:6)]))
  
  data.subs[,update_tax.ident         := fifelse(any(is_taxa_most_prob), T, F), by=qname]
  data.subs[,taxa_most_probable_rank  := rank]
  
  data.subs[,tax.identity             := fifelse(update_tax.ident, taxa_most_probable,        tax.identity)]
  data.subs[,tax.identity.level       := fifelse(update_tax.ident, taxa_most_probable_rank,     tax.identity.level)]
  data.subs[,tax.identity_rank_level  := fifelse(update_tax.ident, taxa_most_probable_rank_level, tax.identity_rank_level)]
  
  data[,tax.identity             := NULL]
  data[,tax.identity.level       := NULL]
  data[,tax.identity_rank_level  := NULL]
  
  cols <- c('sample', 'qname', 'aln_nr', 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level')
  data.subs <- data.subs[,..cols]
  
  #cols <- c('sample', 'qname', 'aln_nr', ranks)
  #data <- data[,..cols]
  
  data <- merge(data[,], data.subs[,],
                by = c('sample', 'qname', 'aln_nr'))
  
  return(data)
}




BestAln <- function(taxa.dup, ranks=ranks, thresholds=thresholds) {
  
  # Define the taxonomic ranks in order
  rank_df <- data.frame(rank=ranks, rank_level = 1:length(ranks))
  
  taxa.dup     <- merge(taxa.dup, rank_df, by.x='tax.identity.level', by.y='rank', all.x=T)
  setnames(taxa.dup, 'rank_level', 'tax.identity_rank_level')
  
  taxa.dup[is.na(tax.identity_rank_level), tax.identity_rank_level := 0]
  
  cols_by <- c('sample', 'qname', 'aln_nr', ranks, 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level')
  taxa.most_prob <- taxa.dup[, ..cols_by]
  
  ## 1.
  # Apply a Function that calculate and assign the most probable taxon for a given rank, iteratively for each rank 
  for (rank in ranks) {
    message('Finding probability for ', rank, '...')
    most_probable_taxa <- assign_most_probable(data=taxa.most_prob, rank, thresholds)
    taxa.most_prob     <- merge(taxa.most_prob, most_probable_taxa, 
                                by = c("qname", 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level', rank), all.x = TRUE)
    
  }
  
  ## 2.
  # Apply a Function to update the tax.identity based on probability and level iteratively for each level
  taxa.new.taxident <- data.table(as.data.frame(taxa.most_prob)) #data.table(NULL)
  
  for (rank in ranks) {
    message('Refining tax.identity based on ', rank, '...')
    
    taxa.new.taxident <- update_tax_identity(data=taxa.new.taxident, rank)
    #taxa.newtaxident <- data.table(NULL)
    
  }
  
  ##
  return(taxa.new.taxident)
  
}
