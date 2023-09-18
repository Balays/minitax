
rename_duptaxa <- function(taxa.sum, 
                           ranks = NULL,
                           cols  = c(ranks, 'tax.identity', 'tax.identity.level')
) {
  
  taxa.sum <- data.table(taxa.sum)        
  
  unidf    <- unique(taxa.sum[,..cols])
  
  taxa.sum <- as.data.frame(taxa.sum)
  
  #dups <- dup(unidf$tax.identity)
  #taxa.sum.dup <- taxa.sum[is.element(taxa.sum)]
  
  taxa.sum[is.element(taxa.sum$tax.identity, dup(unidf$tax.identity)), 'tax.identity'] <-
    paste(taxa.sum[is.element(taxa.sum$tax.identity, dup(unidf$tax.identity)),'tax.identity'],
          taxa.sum[is.element(taxa.sum$tax.identity, dup(unidf$tax.identity)),'tax.identity.level'],
          sep='_')
  
  
  taxa.sum <- data.table(taxa.sum)
  
  message('The follwing taxa were renamed, as their names were the same in different taxonomic levels: \n',
          paste0(taxa.sum[is.element(taxa.sum$tax.identity, dup(unidf$tax.identity)),'tax.identity.level'], 
                 collapse = '; ')
  )
  
  return(taxa.sum)
  
  
}